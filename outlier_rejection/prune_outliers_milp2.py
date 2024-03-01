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
import mosek

import networkx as nx

# native copt version
def prune_outliers(y, cad_dist_min, cad_dist_max, noise_bound, noise_bound_time, prioroutliers):
    '''
    y: 3N x L matrix of each keypoint location at each time
    '''
    N = int(y.shape[0] / 3)
    L = y.shape[1]
    y_list = y.T.reshape([N*L,3]).T
    
    # Make a COPT environment
    env = cpt.Envr()
    # set up model
    model = env.createModel("prune")
    # add variables
    x = model.addMVar(N*L, vtype=COPT.BINARY)
    model.setObjective(x.sum(),sense=COPT.MAXIMIZE)

    # prior outliers
    # TODO

    # shape constraints
    for l in range(L):
        yl = y[:,l].reshape((N,3)).T
        yis, yjs = shape_consistency(yl, cad_dist_min, cad_dist_max, noise_bound)
        yis += N*l
        yjs += N*l
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

    # warmstart (TODO)
    # warm_indicies = prune_outliers_clique(y, cad_dist_min, cad_dist_max, noise_bound, noise_bound_time, prioroutliers)
    # warmstart = np.zeros(N*L)
    # warmstart[warm_indicies] = 1.0

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
    inliers = prune_outliers(db['y'], db['cad_dist_min'], db['cad_dist_max'], db['noise_bound'], db['noise_bound_time'], db['prioroutliers'])
    end = time.time()
    print(np.sort(inliers))
    print(end - start)