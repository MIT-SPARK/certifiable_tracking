'''
prune_outliers_milp.py

Use pairwise compatibility test to select the largest set of inliers.

Lorenzo Shaikewitz for SPARK Lab
'''

import numpy as np
import cvxpy as cp
import pickle
import coptpy as cpt
from coptpy import COPT

import networkx as nx

# from prune_outliers import prune_outliers as prune_outliers_clique

# native copt version
def prune_outliers(y, cad_dist_min, cad_dist_max, noise_bound, noise_bound_time, prioroutliers, warmstart):
    '''
    y: 3N x L matrix of each keypoint location at each time
    '''
    N = int(y.shape[0] / 3)
    L = y.shape[1]
    y_list = y.T.reshape([N*L,3]).T
    
    # Make a COPT environment
    envconfig = cpt.EnvrConfig()
    envconfig.set('nobanner', '1')
    env = cpt.Envr(envconfig)

    # set up model
    model = env.createModel("prune")

    # turn off prints
    model.setParam(COPT.Param.Logging,False)
    # add variables
    x = model.addMVar(N*L, vtype=COPT.BINARY)
    model.setObjective(x.sum(),sense=COPT.MAXIMIZE)

    # prior outliers
    for i in prioroutliers:
        x[i].setInfo(COPT.Info.UB,0.0)

    # shape constraints
    for l in range(L):
        yl = y[:,l].reshape((N,3)).T
        yis, yjs = shape_consistency(yl, cad_dist_min, cad_dist_max, noise_bound)
        yis += N*l
        yjs += N*l
        if (len(yis) > 0):
            model.addConstrs(x[yis] + x[yjs] <= 1)

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
                        model.addConstrs(x[p1] + x[p2] + x[q1] + x[q2] <= 3)

    # warmstart
    if warmstart:
        warm_indicies = prune_outliers_clique(y, cad_dist_min, cad_dist_max, noise_bound, noise_bound_time, warmstart)
        warmstart = np.zeros(N*L)
        warmstart[warm_indicies] = 1.0
        model.setMipStart(x,warmstart)

    # solve!
    model.setParam(COPT.Param.TimeLimit, 10.0)
    model.solve()
    x = np.round(model.getValues())
    # prob.solve(verbose=True)

    # pull out inlier indicies
    inliers = np.array(range(N*L))
    inliers = inliers[x == 1]
    return inliers


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

def save(y, cad_dist_min, cad_dist_max, noise_bound, noise_bound_time, prioroutliers, warmstart):
    db = {}
    db['y'] = y
    db['cad_dist_min'] = cad_dist_min
    db['cad_dist_max'] = cad_dist_max
    db['noise_bound'] = noise_bound
    db['noise_bound_time'] = noise_bound_time
    db['prioroutliers'] = prioroutliers
    db['warmstart'] = warmstart
    
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
    inliers = prune_outliers(db['y'], db['cad_dist_min'], db['cad_dist_max'], db['noise_bound'], db['noise_bound_time'], db['prioroutliers'], db['warmstart'])
    end = time.time()
    print(np.sort(inliers))
    print(end - start)