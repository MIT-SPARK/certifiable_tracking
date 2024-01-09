import clipperpy
import numpy as np

def run_clipper(M, C):
    invariant = clipperpy.invariants.PairwiseInvariant()
    params = clipperpy.Params()

    clipper = clipperpy.CLIPPER(invariant, params)
    clipper.set_matrix_data(M, C)
    clipper.solve()
    sol = clipper.get_solution()
    return sol

if __name__ == '__main__':
    M = np.array([[1, 0.272727272727273, 0.181818181818182],
              [0.272727272727273, 1, 0.0909090909090909],
              [0.181818181818182, 0.0909090909090909, 1]])
    C = np.ones(3)

    sol = run_clipper(M,C)
    print(sol)
    print(sol.u)
    print(sol.nodes)