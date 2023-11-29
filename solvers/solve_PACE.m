function soln = solve_PACE(problem,lambda,l)
% Implements PACE with no outliers a specific time step.
%   Analytically remove position & shape, leaving only
%   rotation in the SDP.
% NOTE: THIS HAS A BUG! DO NOT USE!
%
% INPUTS:
% - problem (struct): populated problem data
% - lambda (float): value of lambda to use
% - l (int): time index to use
%
% RETURNS:
% - soln (struct): solution data
% 
% Lorenzo Shaikewitz for SPARK Lab

%% Process inputs
mosekpath = problem.mosekpath;

N = problem.N;
K = problem.K;

B = problem.B; % 3*N x K matrix of b_i(k)
y = problem.y(:,l); % 3*N x 1 matrix of y_i(t_l)

pBound = problem.translationBound;

if ~isfield(problem,'weights'); problem.weights = ones(N,1); end
weights = problem.weights;
weights = weights / problem.noiseBoundSq;

% check lambda
if ((K >= 3*N) && (lambda == 0.0))
    error("lambda must be positive when there are more shapes (K) " + ...
        "than measurement points (3*N)");
end

%% Define objective
% optimization vector
d = 9; % 1 rotation
x = msspoly('x',d);

% convert to matrix
R  = reshape(x,3,3);

% define eliminated variables
% POSITION COEFFICIENTS
yw = zeros(3,1);
for i = 1:N
    yw = yw + weights(i)*y(ib3(i));
end
yw = yw / sum(weights);

Bw = zeros(3, K);
for i = 1:N
    Bw = Bw + weights(i)*B(ib3(i),:);
end
Bw = repmat(Bw / sum(weights),N,1);

% SHAPE
ybar = y;
for i = 1:N
    ybar(ib3(i)) = R'*(y(ib3(i)) - yw);
end
% sumh = kron(eye(N), R')*(ybar);
Bbar = B - Bw;
H = [2*((Bbar'*Bbar) + lambda*eye(K)), ones(K,1);
     ones(1,K), 0];
b = [2*Bbar'*ybar; 1];
invH = H \ b;

c = invH(1:(end-1),:);

% FINALIZE POSITION
p = yw - R * Bw(1:3,:) * c;

% MAIN OPTIMIZATION
prob_obj = 0;
for i = 1:N
    obj2 = R'*y(ib3(i)) - B(ib3(i),:)*c - R'*p;
    prob_obj = prob_obj + obj2' * obj2;
end
prob_obj = prob_obj + lambda*(c'*c);

%% Define constraints
% EQUALITY
% SO(3) constraints
h = so3_constraints(R);

% INEQUALITY
% p in range (p'*p<=pBoundSq)
pBoundSq = pBound^2;
g_p = pBoundSq - p'*p;

% c bound (0<=c<=1)
cBoundSq = 1.0; % should just be 1
g_c = [cBoundSq - c'*c;c];

% g = [];
% g = [g_c;g_v];
g = [g_p; g_c];

%% Complete problem definition
problem.vars = x;
problem.objective = prob_obj;
problem.equality = h; % equality
problem.inequality = g; % inequsality

%% Relax!
kappa = 1; % relaxation order
[SDP,info] = dense_sdp_relax(problem,kappa);

%% Solve using MOSEK
% load("sdpvars.mat");
prob = convert_sedumi2mosek(SDP.sedumi.At,...
                            SDP.sedumi.b,...
                            SDP.sedumi.c,...
                            SDP.sedumi.K);
addpath(genpath(mosekpath))
[~,res] = mosekopt('minimize info',prob);
[Xopt,yopt,Sopt,obj] = recover_mosek_sol_blk(res,SDP.blk);
rmpath(genpath(mosekpath))

% figure; bar(eig(Xopt{1})); % if rank = 1, then relaxation is exact/tight

%% Compare to ground truth
% clip to first eigenvalue
[eigvecs, ~] = eig(Xopt{1});
vecmax = eigvecs(:,end);
% re-normalize so first element is 1
vecmax_normed = vecmax / vecmax(1);
x_est = vecmax_normed(2:end);

% Project to SO(3) and extract results
R = project2SO3(reshape(x_est,[3,3]));

p_est = full(dmsubs(p,x,x_est));
c_est = full(dmsubs(c,x,x_est));

%% Pack into struct
% raw SDP/MOSEK data
soln.raw.Xopt = Xopt;
soln.raw.yopt = yopt;
soln.raw.Sopt = Sopt;
soln.raw.obj = obj;
soln.raw.relax_info = info;

% save estimates
soln.x_est = x_est;

soln.c_est = c_est;
soln.p_est = p_est;
soln.R_est = R;

end