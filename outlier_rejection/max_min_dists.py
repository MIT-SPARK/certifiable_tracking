import numpy as np
import cvxpy as cp

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

