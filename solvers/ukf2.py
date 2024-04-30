'''
A UKF for tracking an object using the constant body velocity and constant angular velocity model.

Lorenzo Shaikewitz for SPARK Lab
'''
import ukfm
from ukfm.geometry import SO3
import numpy as np

# MODEL CLASS
class CONSTBODY:
    class STATE:
        """State of the system.
        keypoints, velocity, rotation rate
        """

        def __init__(self, keypoints, v, dR):
            self.N = keypoints.shape[1]
            self.keypoints = keypoints # 3 x N
            self.v = v # 3
            self.dR = dR # 3 x 3

    class INPUT:
        """Input of the propagation model.
        there is no control input
        """
        def __init__(self):
            pass

    def __init__(self, shape, dt):
        # mapping keypoint -> center of body
        self.shape = shape
        # integration step (s)
        self.dt = dt

    @classmethod
    def f(cls, state, omega, noise, dt):
        """ Propagation function.
        """
        # velocities
        v_new = state.v + noise[:3]
        dR_new = state.dR.dot(SO3.from_rpy(noise[3:6]))

        # new keypoints
        keypoints_new = state.dR.T.dot(state.keypoints + state.v*dt)

        new_state = cls.STATE(
            keypoints=keypoints_new,
            v=v_new,
            dR=dR_new,
        )
        return new_state

    @classmethod
    def h(cls, state):
        """ Observation function.

        .. math::

            h\\left(\\boldsymbol{\\chi}\\right)  = 
            \\mathbf{H} \mathbf{x},

        where 

        .. math::

            \mathbf{H}&= \\begin{bmatrix} 0 & 1 & 0 \\\\  0
            & 0 & 1 \end{bmatrix} \\\\
            \mathbf{x} &= L \\mathbf{C} \mathbf{e}^b

        with :math:`\mathbf{x}` the position of the pendulum.

        :var state: state :math:`\\boldsymbol{\\chi}`.
        """
        y = np.squeeze(state.keypoints.reshape([-1,1]))
        return y

    @classmethod
    def phi(cls, state, xi):
        """Retraction.
        """
        new_state = cls.STATE(
            keypoints=state.keypoints + xi[:(3*state.N)].reshape([3,state.N]),
            v=state.v + xi[(3*state.N):(3*state.N+3)],
            dR=state.dR.dot(SO3.from_rpy(xi[(3*state.N+3):])),
        )
        return new_state

    @classmethod
    def phi_inv(cls, state, hat_state):
        """Inverse retraction.
        """
        xi = np.hstack([
                np.squeeze(state.p - hat_state.p),
                np.squeeze(state.v - hat_state.v),
                SO3.to_rpy(hat_state.R @ (state.R.T)),
                ])
        return xi


def run_ukf(est, shape, dt, P, Q, R):
    '''
    L: number of times to reason over (including initial time)
    p_meas: 3 x 1 x L matrix of position measurements
    R_meas: 3 x 3 x L matrix of rotation measurements
    v_init: initial body velocity
    w_init: initial angular velocity
    Q: propagation noise covariance matrix (6 x 6) -- noise on const. vel. model
    R: measurement noise covariance matrix (6 x 6)
    P: initial uncertainty matrix (12 x 12)
    '''
    ## CREATE MODEL
    # create the model
    model = CONSTBODY(shape, dt)
    L = est.p.shape[2]

    # initial state
    state0 = model.STATE(keypoints=est.keypoints[:,:,0], v=est.v[:,:,0], dR=est.dR[:,:,0])

    ## CREATE UKF
    # sigma point parameters
    alpha = np.array([1e-3]*12)

    ukf = ukfm.UKF(state0=state0, P0=P, f=model.f, h=model.h, Q=Q, R=R,
                phi=model.phi, phi_inv=model.phi_inv, alpha=alpha)
    
    # set variables for recording estimates along the full trajectory
    ukf_states = [state0]
    ukf_Ps = np.zeros((L, 12, 12))
    ukf_Ps[0] = P

    ## FILTERING
    for l in range(1, L):
        # propagation
        ukf.propagation([], model.dt) # empty input

        # update only if a measurement is received
        ukf.update(est.keypoints[:,:,l])

        # save estimates
        ukf_states.append(ukf.state)
        ukf_Ps[l] = ukf.P
        est.keypoints[:,:,l+1] = ukf.state.keypoints
        est.v[:,:,l+1] = ukf.state.v
        est.dR[:,:,l+1] = ukf.state.dR

    # Process states into numpy array for passing back to MATLAB
    return est


from scipy.spatial.transform import Rotation as R_geo

class DataHandler:
    def __init__(self,N,L):
        self.keypoints = np.zeros(3,N,L)
        self.v = np.zeros(3,1,L-1)
        self.dR = np.zeros(3,3,L-1)

    def __init__(self,keypoints,v,dR):
        self.keypoints = keypoints
        self.v = v
        self.dR = dR

# GENERATE DATA FOR UKF TEST
def generate_data(L, N):
    dt = 1.0

    covar_velocity = 0.05
    covar_rotrate = 0.05
    covar_measure = 0.05
    R = np.diag([covar_measure]*N) # measurement noise: keypoint noise
    Q = np.diag([covar_velocity]*3 + [covar_rotrate]*3) # propagation noise: noise added to velocity, rotation rate
    P = np.diag([covar_measure]*N + [0]*3 + [0]*3) # (initial) state uncertanty: covariance of initial estimate

    gt = DataHandler(L)
    est = DataHandler(L)
    ## generate initial points
    # starting shape lib

    # starting keypoints, velocity, rotation rate
    gt.p[:,:,0] = np.random.randn([3,1])
    gt.R[:,:,0] = np.random.randn([3,1])
    gt.keypoints[:,:,0] = np.random.randn([3,N]) + gt.p[:,:,0]
    gt.v[:,:,0] = np.random.rand(3)*5.0
    gt.dR[:,:,0] = R_geo.random().as_matrix()

    ## generate trajectory
    for l in range(1,L):
        # Ground truth: these should match the f function
        gt.v[:,:,l] = gt.v[:,:,l-1] + np.random.randn([3,1])*np.sqrt(covar_velocity)
        gt.dR[:,:,l] = gt.dR[:,:,l-1].dot(SO3.from_rpy(np.random.randn([3,1])*np.sqrt(covar_rotrate)))
        gt.keypoints[:,:,l] = gt.dR[:,:,l-1].dot(gt.keypoints[:,:,l-1] + gt.v[:,:,l-1]*dt) # TODO: SHOULD TRANSFORM CENTROID INSTEAD
        
        # Estimated: keypoints only, should match the phi function
        est.keypoints[:,:,l] = gt.keypoints[:,:,l] + np.random.randn([3,N])*np.sqrt(covar_measure)

    # initialize est
    est.v[:,:,0] = gt.v[:,:,0]
    est.dR[:,:,0] = gt.dR[:,:,0]
    est.keypoints[:,:,0] = gt.keypoints[:,:,0] + np.random.randn([3,N])*np.sqrt(covar_measure)

    return gt, est, P, Q, R, dt



import pickle

def save(keypoints, v, dR, dt, P, Q, R):
    db = {}
    db['dt'] = dt
    db['keypoints'] = keypoints
    db['v'] = v
    db['dR'] = dR
    db['Q'] = Q
    db['R'] = R
    db['P'] = P
    
    dbfile = open('examplePickle', 'ab')
    # source, destination
    pickle.dump(db, dbfile)                    
    dbfile.close()

def read_pickle():
    dbfile = open('examplePickle', 'rb')
    db = pickle.load(dbfile)
    dbfile.close()
    dt = db['dt']
    P = db['P']
    Q = db['Q']
    R = db['R']

    est = DataHandler(db['keypoints'], db['v'], db['dR'])

    return est, P, Q, R, dt

    
if __name__ == '__main__':
    gt, est, P, Q, R, dt = generate_data(L=12, N=10)
    shape = gt.keypoints[:,:,0]
    
    est = run_ukf(est, shape, dt, P, Q, R)
    print(p,r)