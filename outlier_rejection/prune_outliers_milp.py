'''
prune_outliers_milp.py

Use pairwise compatibility test to select the largest set of inliers.

Lorenzo Shaikewitz for SPARK Lab
'''

import numpy as np
import cvxpy as cp
import pickle
import sys
import mosek

import networkx as nx

def prune_outliers_clique(y, cad_dist_min, cad_dist_max, noise_bound, noise_bound_time, prioroutliers):
    '''
    Compute an approximate inlier set using ROBIN
    '''
    N = int(y.shape[0]/3)
    L = y.shape[1]
    y_list = y.T.reshape([N*L,3]).T

    # Graph of single-time compatibility with shape lib
    graphs = {}
    for l in range(L):
        g = nx.Graph()
        yl = y[:,l].reshape([N,3]).T
        if (len(prioroutliers) < l + 1):
            outliers = []
        else:
            outliers = prioroutliers[l]
        shape_consistency_clique(g, yl, cad_dist_min, cad_dist_max, noise_bound, outliers)
        graphs[l] = g

    # Graph of keypoint compatbility across times
    for l1 in range(L-1):
        for l2 in range(l1+1,L):
            for i1 in range(N-1):
                for i2 in range(i1+1,N):
                    if not (graphs[l1].has_edge(i1,i2) and graphs[l2].has_edge(i1,i2)):
                        continue

                    p1 = l1*N + i1
                    p2 = l1*N + i2
                    d1 = np.linalg.norm(y_list[:,p1]-y_list[:,p2])
                    q1 = l2*N + i1
                    q2 = l2*N + i2
                    d2 = np.linalg.norm(y_list[:,q1]-y_list[:,q2])
                    if (abs(d1-d2) < 4*noise_bound_time):
                        graphs[l1].nodes[i1]['weight'] += 1.0
                        graphs[l1].nodes[i2]['weight'] += 1.0
                        graphs[l2].nodes[i1]['weight'] += 1.0
                        graphs[l2].nodes[i2]['weight'] += 1.0

                    graphs[l1].nodes[i1]['count'] += 1
                    graphs[l1].nodes[i2]['count'] += 1
                    graphs[l2].nodes[i1]['count'] += 1
                    graphs[l2].nodes[i2]['count'] += 1

    ## solve
    inlier_indices = []
    for l in range(L):
        g = graphs[l]
        for n in range(len(g.nodes())):
            g.nodes[n]['weight'] = int(g.nodes[n]['weight'] / g.nodes[n]['count'] * 100)

        [clique, _] = nx.max_weight_clique(g)
        for n in clique:
            inlier_indices.append(n + l*N)
        # print(l+1)
        # inlier_indices.sort()
        # print(inlier_indices)
        # draw(g)

    return inlier_indices

def shape_consistency_clique(g, tgt, cad_dist_min, cad_dist_max, noise_bound, prioroutliers):
    '''
    Build graph of keypoint measurements by consistency with shape library
    '''
    N = tgt.shape[1]
    si, sj = np.meshgrid(np.arange(N), np.arange(N))
    mask_uppertri = (sj > si)
    si = si[mask_uppertri]
    sj = sj[mask_uppertri]

    # distances || tgt_j - tgt_i ||
    tgt_dist_ij = np.linalg.norm(
        tgt[:, sj] - tgt[:, si], axis=0)  # shape (n-1)_tri

    allEdges = np.arange(si.shape[0])
    check1 = tgt_dist_ij >= (cad_dist_min - 2 * noise_bound)
    check2 = tgt_dist_ij <= (cad_dist_max + 2 * noise_bound)
    mask_compatible = check1 & check2
    validEdges = allEdges[mask_compatible]
    sdata = np.zeros_like(si)
    sdata[mask_compatible] = 1

    comp_mat = np.zeros((N, N))
    comp_mat[si, sj] = sdata

    nodes = {}
    # creating a Graph in robin
    for i in range(N):
        g.add_node(i,weight=1.0,count=1)

    for edge_idx in validEdges:
        # print(f'Add edge between {si[edge_idx] + lidx} and {sj[edge_idx] + lidx}.')
        if (si[edge_idx] in prioroutliers) or (sj[edge_idx] in prioroutliers):
            continue
        g.add_edge(si[edge_idx], sj[edge_idx])


# native mosek version
def prune_outliers(y, cad_dist_min, cad_dist_max, noise_bound, noise_bound_time, prioroutliers):
    '''
    y: 3N x L matrix of each keypoint location at each time
    '''
    N = int(y.shape[0] / 3)
    L = y.shape[1]
    y_list = y.T.reshape([N*L,3]).T

    # shape constraints
    all_yis = []
    all_yjs = []
    for l in range(L):
        yl = y[:,l].reshape((N,3)).T
        yis, yjs = shape_consistency(yl, cad_dist_min, cad_dist_max, noise_bound)
        yis += N*l
        yjs += N*l
        all_yis = np.append(all_yis,yis)
        all_yjs = np.append(all_yjs,yjs)

    # rigid body constraints
    p1s = []; p2s = []
    q1s = []; q2s = []
    for l1 in range(L-1):
        for l2 in range(l1+1,L):
            for i1 in range(N-1):
                for i2 in range(i1+1,N):
                    p1 = l1*N + i1
                    p2 = l1*N + i2
                    q1 = l2*N + i1
                    q2 = l2*N + i2
                    d1 = np.linalg.norm(y_list[:,p1]-y_list[:,p2])
                    d2 = np.linalg.norm(y_list[:,q1]-y_list[:,q2])
                    if (abs(d1-d2) >= 4*noise_bound_time):
                        p1s = np.append(p1s,p1); p2s = np.append(p2s,p2)
                        q1s = np.append(q1s,q1); q2s = np.append(q2s,q2)

    # warmstart
    warm_indicies = prune_outliers_clique(y, cad_dist_min, cad_dist_max, noise_bound, noise_bound_time, prioroutliers)
    warmstart = np.zeros(N*L)
    warmstart[warm_indicies] = 1.0

    # solve!
    x = solve_mosek_native(N*L, prioroutliers, [all_yis, all_yjs], [p1s, p2s, q1s, q2s], warmstart)
    # prob.solve(verbose=True)

    # pull out inlier indicies
    inliers = np.array(range(N*L))
    inliers = inliers[x == 1]
    return inliers

# Define a stream printer to grab output from MOSEK
def streamprinter(text):
    sys.stdout.write(text)
    sys.stdout.flush()
inf = 0.0 # for symbolic purposes

def solve_mosek_native(numvar, prior_indicies, shape_indicies, time_indicies, warmstart):
    # Make a MOSEK environment
    with mosek.Env() as env:
        # Attach a printer to the environment
        env.set_Stream(mosek.streamtype.log, streamprinter)
        # Create a task
        with env.Task(0, 0) as task:
            # Attach a printer to the task
            task.set_Stream(mosek.streamtype.log, streamprinter)

            # add empty variables/constraints
            numcon = len(shape_indicies[0]) + len(time_indicies[0])
            cur_con = 0
            task.appendcons(numcon) # constraints
            task.appendvars(numvar) # variables

            # define the objective
            [task.putcj(j, 1.0) for j in range(numvar)]
            task.putobjsense(mosek.objsense.maximize)

            # integer variables
            task.putvartypelist(list(range(numvar)), [mosek.variabletype.type_int]*numvar)

            # set up variable bounds (default to [0,1])
            bkx = [mosek.boundkey.ra]*numvar
            blx = [0.0]*numvar
            bux = [1.0]*numvar

            # set up constraint bounds (default to <= 3)
            bkc = [mosek.boundkey.up]*numcon
            blc = [-inf]*numcon
            buc = [1.0]*len(shape_indicies[0]) + [3.0]*len(time_indicies[0])

            # prior outliers: set variable = 0
            for j in prior_indicies:
                bkx[j] = mosek.boundkey.fx
                bux[j] = 0.0
            
            # shape consistency: xi + xj <= 1
            yis = shape_indicies[0]
            yjs = shape_indicies[1]
            for idx in range(len(yis)):
                i = int(yis[idx])
                j = int(yjs[idx])
                # xi + xj
                task.putarow(cur_con, [i,j], [1.0, 1.0])
                # <= 1 already set
                cur_con += 1

            # time consistency: x[p1]+x[p2]+x[q1]+x[q2] <= 3
            p1s = time_indicies[0]
            p2s = time_indicies[1]
            q1s = time_indicies[2]
            q2s = time_indicies[3]
            for idx in range(len(p1s)):
                p1 = int(p1s[idx])
                p2 = int(p2s[idx])
                q1 = int(q1s[idx])
                q2 = int(q2s[idx])
                # x[p1]+x[p2]+x[q1]+x[q2]
                task.putarow(cur_con, [p1,p2,q1,q2], [1.0]*4)
                # <= 3 already set
                cur_con += 1

            # save variable bounds
            task.putvarboundlist(range(numvar), bkx, blx, bux)
            # save constraint bounds
            task.putconboundlist(range(numcon), bkc, blc, buc)

            # Set max solution time
            task.putdouparam(mosek.dparam.mio_max_time, 60.0);

            # warm start
            task.putxx(mosek.soltype.itg, warmstart)

            # Optimize!
            task.optimize()

            xx = task.getxx(mosek.soltype.itg)
            return np.round(xx)

            # Print a summary containing information
            # about the solution for debugging purposes
            # task.solutionsummary(mosek.streamtype.msg)
            # prosta = task.getprosta(mosek.soltype.itg)
            # solsta = task.getsolsta(mosek.soltype.itg)

            # Output a solution         
            # xx = task.getxx(mosek.soltype.itg)
            # if solsta in [mosek.solsta.integer_optimal]:
            #     print("Optimal solution: %s" % xx)
            # elif solsta == mosek.solsta.prim_feas:
            #     print("Feasible solution: %s" % xx)
            # elif mosek.solsta.unknown:
            #     if prosta == mosek.prosta.prim_infeas_or_unbounded:
            #         print("Problem status Infeasible or unbounded.\n")
            #     elif prosta == mosek.prosta.prim_infeas:
            #         print("Problem status Infeasible.\n")
            #     elif prosta == mosek.prosta.unkown:
            #         print("Problem status unkown.\n")
            #     else:
            #         print("Other problem status.\n")
            # else:
            #     print("Other solution status")


def shape_consistency(tgt, cad_dist_min, cad_dist_max, noise_bound):
    N = tgt.shape[1]
    si, sj = np.meshgrid(np.arange(N), np.arange(N))
    mask_uppertri = (sj > si)
    si = si[mask_uppertri]
    sj = sj[mask_uppertri]

    # distances || tgt_j - tgt_i ||
    tgt_dist_ij = np.linalg.norm(
        tgt[:, sj] - tgt[:, si], axis=0)  # shape (n-1)_tri

    allEdges = np.arange(si.shape[0])
    check1 = tgt_dist_ij >= (cad_dist_min - 2 * noise_bound)
    check2 = tgt_dist_ij <= (cad_dist_max + 2 * noise_bound)
    mask_compatible = check1 & check2
    invalidEdges = allEdges[np.bitwise_not(mask_compatible)]
    sdata = np.zeros_like(si)
    sdata[mask_compatible] = 1

    comp_mat = np.zeros((N, N))
    comp_mat[si, sj] = sdata

    yis = si[invalidEdges]
    yjs = sj[invalidEdges]
    return yis, yjs

# cvx version
def prune_outliers_cvx(y, cad_dist_min, cad_dist_max, noise_bound, noise_bound_time, prioroutliers):
    '''
    y: 3N x L matrix of each keypoint location at each time
    '''
    N = int(y.shape[0] / 3)
    L = y.shape[1]
    y_list = y.T.reshape([N*L,3]).T

    # define variables, objective
    x = cp.Variable(N*L, boolean=True)

    objective = cp.Maximize(cp.sum(x))

    constraints_shape = []
    constraints_time = []

    # prior outliers
    if len(prioroutliers) > 0:
        # TODO: check if this is right
        constraints_prioroutliers = [x[prioroutliers] == 0]
    else:
        constraints_prioroutliers = []

    # shape constraints
    for l in range(L):
        yl = y[:,l].reshape((N,3)).T
        yis, yjs = shape_consistency(yl, cad_dist_min, cad_dist_max, noise_bound)
        yis += N*l
        yjs += N*l
        for i in range(len(yis)):
            constraints_shape.append(x[yis[i]]+x[yjs[i]] <= 1)

    # rigid body constraints
    for l1 in range(L-1):
        for l2 in range(l1+1,L):
            for i1 in range(N-1):
                for i2 in range(i1+1,N):
                    p1 = l1*N + i1
                    p2 = l1*N + i2
                    q1 = l2*N + i1
                    q2 = l2*N + i2
                    d1 = np.linalg.norm(y_list[:,p1]-y_list[:,p2])
                    d2 = np.linalg.norm(y_list[:,q1]-y_list[:,q2])
                    if (abs(d1-d2) >= 4*noise_bound_time):
                        constraints_time.append(x[p1]+x[p2]+x[q1]+x[q2] <= 3)

    # solve!
    prob = cp.Problem(objective, 
                      constraints_prioroutliers + 
                      constraints_shape + constraints_time)
    prob.solve(solver='COPT', verbose=True)
    # prob.solve(verbose=True)

    # pull out inlier indicies
    inliers = np.array(range(N*L))
    inliers = inliers[x.value == 1]
    return inliers

def save(y, cad_dist_min, cad_dist_max, noise_bound, noise_bound_time, prioroutliers):
    db = {}
    db['y'] = y
    db['cad_dist_min'] = cad_dist_min
    db['cad_dist_max'] = cad_dist_max
    db['noise_bound'] = noise_bound
    db['noise_bound_time'] = noise_bound_time
    db['prioroutliers'] = prioroutliers
    
    dbfile = open('test_pickle', 'wb')
    # source, destination
    pickle.dump(db, dbfile)
    dbfile.close()
    
import time

if __name__ == '__main__':
    dbfile = open('test_pickle', 'rb')
    db = pickle.load(dbfile)
    dbfile.close()
    
    start = time.time()
    inliers = prune_outliers_cvx(db['y'], db['cad_dist_min'], db['cad_dist_max'], db['noise_bound'], db['noise_bound_time'], db['prioroutliers'])
    end = time.time()
    print(np.sort(inliers))
    print(end - start)