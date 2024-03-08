'''
prune_outliers.py

A graph-theoretic approach to excluding outliers.
Use compatibility test to determine if two points can simultaneously be inliers.
TODO: NO SUPPORT FOR MISSING KEYPOINTS

Lorenzo Shaikewitz for SPARK Lab
Adapted from code written by Jingnan Shi.
'''

import numpy as np
import cvxpy as cp

import robin_py

def robin_prune_outliers(tgt, cad_dist_min, cad_dist_max, noise_bound, prioroutliers, method='maxclique'):
    '''
    First form a compatibility graph and then 
    Use robin to select inliers

    tgt is 3 x N (all keypoints at given time)
    get cad_dist min and max from helper functions
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

    # creating a Graph in robin
    g = robin_py.AdjListGraph()
    for i in range(N):
        g.AddVertex(i)

    for edge_idx in validEdges:
        if (si[edge_idx] in prioroutliers) or (sj[edge_idx] in prioroutliers):
            continue
        # print(f'Add edge between {si[edge_idx]} and {sj[edge_idx]}.')
        g.AddEdge(si[edge_idx], sj[edge_idx])

    if method == "maxclique":
        inlier_indices = robin_py.FindInlierStructure(
            g, robin_py.InlierGraphStructure.MAX_CLIQUE)
    elif method == "maxcore":
        inlier_indices = robin_py.FindInlierStructure(
            g, robin_py.InlierGraphStructure.MAX_CORE)
    else:
        raise RuntimeError('Prune outliers only support maxclique and maxcore')

    # adj_mat = g.GetAdjMat()

    return inlier_indices, comp_mat
    
'''
Helper functions for min and max distances
'''
def compute_min_max_distances(cad_kpts):
    '''
    cad_kpts is K x 3 x N
    '''
    # print('Computing upper and lower bounds in cad pairwise distances...')

    K = cad_kpts.shape[0]
    N = cad_kpts.shape[2]
    si, sj = np.meshgrid(np.arange(N), np.arange(N))
    mask_uppertri = (sj > si)
    si = si[mask_uppertri]
    sj = sj[mask_uppertri]

    cad_TIMs_ij = cad_kpts[:, :, sj] - cad_kpts[:, :, si]  # shape K by 3 by (n-1)_tri

    # compute max distances
    cad_dist_k_ij = np.linalg.norm(cad_TIMs_ij, axis=1)  # shape K by (n-1)_tri
    cad_dist_max_ij = np.max(cad_dist_k_ij, axis=0)

    # compute min distances
    cad_dist_min_ij = []
    num_pairs = cad_TIMs_ij.shape[2]
    one_tenth = num_pairs / 10
    for i in range(num_pairs):
        tmp = cad_TIMs_ij[:, :, i].T
        min_dist = minimum_distance_to_convex_hull(tmp)
        cad_dist_min_ij.append(min_dist)
        # if i % one_tenth == 1:
        #     print(f'{i}/{num_pairs}.')

    cad_dist_min_ij = np.array(cad_dist_min_ij)

    print('Done')

    return cad_dist_min_ij, cad_dist_max_ij

def minimum_distance_to_convex_hull(A):
    '''
    A is shape 3 by K, compute the minimum distance from the origin to the convex hull of A
    '''
    K = A.shape[1]
    P = A.T @ A
    one = np.ones((K, 1))
    # Use CVXPY to solve
    x = cp.Variable(K)
    prob = cp.Problem(cp.Minimize(cp.quad_form(x, P)),
                      [x >= 0,
                       one.T @ x == 1])
    prob.solve(solver='ECOS', verbose=False)
    x_val = x.value
    min_distance = np.linalg.norm(A @ x_val)
    return min_distance


'''
Invariance-based outlier pruning for time series
TODO: FIX
'''
def prune_outliers(y, cad_dist_min, cad_dist_max, noise_bound, noise_bound_time, prioroutliers, method="maxclique"):
    '''
    Compute an approximate inlier set using ROBIN
    '''
    N = int(y.shape[0]/3)
    L = y.shape[1]
    y_list = y.T.reshape([N*L,3]).T

    g = robin_py.AdjListGraph()
    # Graph of single-time compatibility with shape lib
    for l in range(L):
        yl = y[:,l].reshape([N,3]).T
        if (len(prioroutliers) < l + 1):
            outliers = []
        else:
            outliers = prioroutliers[l]
        shape_consistency(g, l*N, yl, cad_dist_min, cad_dist_max, noise_bound)

    test = g.GetAdjMat()
    # Graph of keypoint compatbility across times
    for l1 in range(L-1):
        for l2 in range(l1+1,L):
            for i1 in range(N-1):
                for i2 in range(i1+1,N):
                    p1 = l1*N + i1
                    p2 = l1*N + i2
                    d1 = np.linalg.norm(y_list[:,p1]-y_list[:,p2])
                    q1 = l2*N + i1
                    q2 = l2*N + i2
                    d2 = np.linalg.norm(y_list[:,q1]-y_list[:,q2])
                    if (p1 == 0):
                        print(p1,p2,q1,q2)
                    if (abs(d1-d2) < noise_bound_time):
                        g.AddEdge(p1,q1)
                        g.AddEdge(p1,q2)
                        g.AddEdge(p2,q1)
                        g.AddEdge(p2,q2)
                    else:
                        if (g.HasEdge(p1,p2) and g.HasEdge(q1,q2)):
                            print(p1,p2,q1,q2,abs(d1-d2))
    test2 = g.GetAdjMat()

    ## solve
    if method == "maxclique":
        inlier_indices = robin_py.FindInlierStructure(
            g, robin_py.InlierGraphStructure.MAX_CLIQUE)
    elif method == "maxcore":
        inlier_indices = robin_py.FindInlierStructure(
            g, robin_py.InlierGraphStructure.MAX_CORE)
    else:
        raise RuntimeError('Prune outliers only support maxclique and maxcore')

    return inlier_indices

def shape_consistency(g, lidx, tgt, cad_dist_min, cad_dist_max, noise_bound):
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

    # creating a Graph in robin
    for i in range(N):
        g.AddVertex(i + lidx)

    for edge_idx in validEdges:
        # print(f'Add edge between {si[edge_idx] + lidx} and {sj[edge_idx] + lidx}.')
        g.AddEdge(si[edge_idx] + lidx, sj[edge_idx] + lidx)


import pickle

def save(y, cad_dist_min, cad_dist_max, noise_bound, method="maxclique"):
    db = {}
    db['y'] = y
    db['cad_dist_min'] = cad_dist_min
    db['cad_dist_max'] = cad_dist_max
    db['noise_bound'] = noise_bound
    db['method'] = method
    
    dbfile = open('examplePickle', 'ab')
    # source, destination
    pickle.dump(db, dbfile)                    
    dbfile.close()
    
if __name__ == '__main__':
    dbfile = open('examplePickle', 'rb')    
    db = pickle.load(dbfile)
    dbfile.close()
    
    inliers = prune_outliers(db['y'], db['cad_dist_min'], db['cad_dist_max'], db['noise_bound'], db['noise_bound'], "maxclique")
    print(np.sort(inliers))