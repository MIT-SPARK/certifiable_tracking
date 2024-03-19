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
        Position, rotation, velocity, angular velocity
        """

        def __init__(self, p, R, v, w):
            self.p = p
            self.R = R
            self.v = v
            self.w = w

    class INPUT:
        """Input of the propagation model.

        The input is a position and rotation measurement
        """
        def __init__(self):
            pass

    def __init__(self, L, dt):
        # sequence time (s)
        self.T = L*dt
        # number of sequences
        self.L = L
        # integration step (s)
        self.dt = dt

    @classmethod
    def f(cls, state, omega, noise, dt):
        """ Propagation function.
        """
        v_new = state.v + noise[:3]
        w_new = state.w + noise[3:6]

        # rotation
        # R_new = (state.R).dot(SO3.exp((w_new)*dt+noise[3:6]))
        R_new = (state.R).dot(SO3.exp((w_new)*dt))
        # translation
        # p_new = state.p + state.R @ v_new*dt + noise[:3]
        p_new = state.p + state.R @ v_new*dt

        new_state = cls.STATE(
            p = p_new,
            R = R_new,
            v = v_new,
            w = w_new,
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
        y = np.squeeze(state.p)
        # y = np.hstack([y, SO3.to_rpy(state.R)])
        y = np.hstack([y, SO3.log(state.R)])
        return y

    @classmethod
    def phi(cls, state, xi):
        """Retraction.
        """
        new_state = cls.STATE(
            p=state.p + xi[:3],
            R=state.R @ (SO3.exp(xi[3:6])),
            v=state.v + xi[6:9],
            w=state.w + xi[9:12],
        )
        return new_state

    @classmethod
    def phi_inv(cls, state, hat_state):
        """Inverse retraction.
        """
        xi = np.hstack([
                np.squeeze(state.p - hat_state.p),
                SO3.log(hat_state.R @ (state.R.T)),
                np.squeeze(state.v - hat_state.v),
                np.squeeze(state.w - hat_state.w),
                ])
        return xi
    

def run_ukf(dt, L, p_meas, R_meas, v_init, w_init, Q, R, P):
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
    L = int(L)
    model = CONSTBODY(L, dt)

    ## Convert measurements into states (R -> rpy)
    ys = [] # TODO
    for l in range(L):
        # rpy = SO3.to_rpy(R_meas[:,:,l])
        rpy = SO3.log(R_meas[:,:,l])
        ys.append(np.hstack([np.squeeze(p_meas[:,:,l]), rpy]))

    # initial state
    state0 = model.STATE(R=R_meas[:,:,0], p=np.squeeze(p_meas[:,:,0]), v=v_init, w=w_init)

    ## CREATE UKF
    # sigma point parameters
    alpha = np.array([1e-3]*12) # TODO: 12?

    # P = np.diag([1,1,1,0.0,0.,0.0,1,1,1,0.005,0.005,0.005])
    # Q = np.diag([1,1,1,0.005,0.005,0.005])
    # R = np.diag([1,1,1,0.005,0.005,0.005])*0

    ukf = ukfm.UKF(state0=state0, P0=P, f=model.f, h=model.h, Q=Q, R=R,
                phi=model.phi, phi_inv=model.phi_inv, alpha=alpha)
    
    # set variables for recording estimates along the full trajectory
    ukf_states = [state0]
    ukf_Ps = np.zeros((model.L, 12, 12))
    ukf_Ps[0] = P

    p_est = np.zeros([3,1,model.L-1])
    R_est = np.zeros([3,3,model.L-1])

    ## FILTERING
    for l in range(1, model.L):
        # propagation
        ukf.propagation([], model.dt) # empty input
        # update only if a measurement is received

        ukf.update(ys[l])

        # save estimates
        ukf_states.append(ukf.state)
        ukf_Ps[l] = ukf.P
        p_est[:,:,l-1] = np.reshape(ukf.state.p,[3,1])
        R_est[:,:,l-1] = ukf.state.R

    # Process states into numpy array for passing back to MATLAB
    return p_est, R_est


import pickle

def save(dt, L, p_meas, R_meas, v_init, w_init, Q, R, P):
    db = {}
    db['dt'] = dt
    db['L'] = L
    db['p_meas'] = p_meas
    db['R_meas'] = R_meas
    db['v_init'] = v_init
    db['w_init'] = w_init
    db['Q'] = Q
    db['R'] = R
    db['P'] = P
    
    dbfile = open('examplePickle', 'ab')
    # source, destination
    pickle.dump(db, dbfile)                    
    dbfile.close()
    
if __name__ == '__main__':
    dbfile = open('examplePickle', 'rb')    
    db = pickle.load(dbfile)
    dbfile.close()
    
    p, r = run_ukf(db['dt'], db['L'], db['p_meas'], db['R_meas'], db['v_init'], db['w_init'], db['Q'], db['R'], db['P'])
    print(p,r)