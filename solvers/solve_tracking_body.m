function soln = solve_tracking_body(problem)
% Solves const. vel. (body frame) optimization exactly via SDP
%   Assume the *body frame* velocity is constant and object is spinning.
%   Result: spiral trajectory.
%   Analytically remove velocity & shape. SDP variables are
%   * rotated position (s)
%   * body velocity (v)
%   * rotation (R)
%   * rotation change (dR)
%   MINIMALLY SYMBOLIC VERSION
%
% INPUTS:
% - problem (struct): populated problem data
%
% RETURNS:
% - soln (struct): solution data
% 
% Lorenzo Shaikewitz for SPARK Lab

%% Process inputs
% mosekpath = problem.mosekpath;

N = problem.N_VAR;
K = problem.K;
L = problem.L;

B = problem.B; % 3*N x K matrix of b_i(k)
y = problem.y; % 3*N x L matrix of y_i(t_l)
dt = problem.dt;

% Weights
% TODO: scale weights/covars by noisebound?
W = problem.covar_measure.^(-1); % N x L matrix of w_il
lambda = problem.lambda; % scalar
wv = problem.covar_velocity.^(-1); % L-2 vector
wd = problem.kappa_rotrate;  % L-1 vector

pBound = problem.translationBound;
vBound = problem.velocityBound;

% check lambda
if ((K >= 3*N) && (lambda == 0.0))
    error("lambda must be positive when there are more shapes (K) " + ...
        "than measurement points (3*N)");
end

%% Define objective
% optimization vector
d = 9*(2*L - 1) + 3*L + 3*(L-1); % 2L - 1 rotations, L rotated positions, L-1 body velocities
% 2L-1 rotations: L rotations, L-1 delta rotations
x = msspoly('x',d);

% pull out individual variables
r  = x(1:(9*L));
dr = x((9*L + 1):(9*L + 9*(L-1)));
s  = x((18*L - 9 + 1):(18*L - 9 + 3*L));
v = x((21*L - 9 + 1):(21*L - 9 + 3*(L-1)));

% convert to useful form
R  = reshape(r ,3,3*L)';
dR = reshape(dr,3,3*(L-1))';
for l = 1:L
    R(ib3(l),:) =  R(ib3(l),:)';
    if (l < L)
       dR(ib3(l),:) = dR(ib3(l),:)';
    end
end

% % define eliminated variables
% % SHAPE
% cbar = ones(K,1) / K;
% sumh = msspoly(zeros(3*N,1));
% sumW = zeros(N,1);
% for l = 1:L
%     for i = 1:N
%         sumh(ib3(i), :) = sumh(ib3(i), :) + W(i,l)*( ...
%                           R(ib3(l),:)'*y(ib3(i),l) ...
%                           - s(ib3(l)));
% 
%         sumW(i) = sumW(i) + W(i,l);
%     end
% end
% sumW = diag(reshape(repmat(sumW,1,3)',3*N,1));
% H = [2*((B'*sumW*B) + lambda*eye(K)), ones(K,1);
%      ones(1,K), 0];
% b = [2*B'*sumh + 2*lambda*cbar; 1];
% invH = H \ b;
% 
% c = invH(1:(end-1),:);
% 
% % MAIN OPTIMIZATION
% prob_obj = 0;
% for l = 1:L
%     for i = 1:N
%         obj2 = R(ib3(l),:)'*y(ib3(i),l) - B(ib3(i),:)*c - s(ib3(l));
%         prob_obj = prob_obj + W(i,l) * (obj2' * obj2);
%     end
% end
% % c regularization
% prob_obj = prob_obj + lambda*((c - cbar)'*(c - cbar));
% prob_obj2 = prob_obj;
% for l = 2:L-1
%     % delta v
%     delv = v(ib3(l)) - v(ib3(l-1));
%     prob_obj = prob_obj + wv(l-1)*(delv'*delv);
%     % dR
%     deldR = reshape(dR(ib3(l),:) - dR(ib3(l-1),:),9,1);
%     prob_obj = prob_obj + wd(l-1)*(deldR'*deldR);
% end

%% Symbol-free version
Y = reshape(y,[3,N,L]);
P = kron(eye(L),vecRt_r);

% SHAPE (VERIFIED)
cbar = ones(K,1) / K;

% schur complement inverse
sumW = zeros(N,1);
for l = 1:L
    for i = 1:N
        sumW(i) = sumW(i) + W(i,l);
    end
end
sumW = diag(reshape(repmat(sumW,1,3)',3*N,1));
H11 = 2*((B'*sumW*B) + lambda*eye(K));
invH11 = inv(H11);
H12 = ones(K,1);
invS = inv(-H12'*invH11*H12);
G = invH11 + invH11*H12*invS*H12'*invH11;
g = -invH11*H12*invS;
% test:
% c2 = G*(2*B'*sumh + 2*lambda*cbar) + g;

% simplify shape
ycr_coeff = zeros(3*N,9*L);
ycs_coeff = kron(W,eye(3));
for l = 1:L
    dw = reshape(repmat(W(:,l)',3,1),3*N,1);
    ycr_coeff(:,(9*(l-1)+1):(9*l)) = diag(dw)*kron(Y(:,:,l)',eye(3));
end

Cr = 2*G*B'*ycr_coeff*P;
Cs = 2*G*B'*ycs_coeff;
gbar = 2*lambda*G*cbar+g;
% test:
% c2 = Cr*r - Cs*s + gbar;

% MEASUREMENT TERM (VERIFIED)
rar_coeff = zeros(3*N*L,3*L);
car_coeff = zeros(3*N*L,K);
sas_coeff = zeros(3*N*L,3*L);
for l = 1:L
    dw = reshape(repmat(W(:,l)',3,1),3*N,1).^(1/2);

    rar_coeff((3*N*(l-1)+1):(3*N*l),(9*(l-1)+1):(9*l)) = ...
           diag(dw)*kron(Y(:,:,l)',eye(3));

    car_coeff((3*N*(l-1)+1):(3*N*l),:) = diag(dw)*B;

    sas_coeff((3*N*(l-1)+1):(3*N*l),(3*(l-1)+1):(3*l)) = ...
        ycs_coeff(:,(3*(l-1)+1):(3*l)).^(1/2);
end

Ar = rar_coeff*P - car_coeff*Cr;
As = -(sas_coeff - car_coeff*Cs);
Ag = -car_coeff*gbar;
% test
% meas_term = (Ar*r + As*s + Ag)'*(Ar*r + As*s + Ag);

% SHAPE REGULARIZATION TERM
if (lambda > 0)
    Ar = [Ar;  sqrt(lambda)*Cr];
    As = [As; -sqrt(lambda)*Cs];
    Ag = [Ag;  sqrt(lambda)*(gbar - cbar)];
end

% VELOCITY TERM
dwv = reshape(repmat(wv',3,1),3*(L-2),1).^(1/2);
eye3LL = [zeros(3*(L-2),3), eye(3*(L-2))];
eye3LR = [eye(3*(L-2)), zeros(3*(L-2),3)];

Av = diag(dwv)*(eye3LL - eye3LR);

% ROTATION RATE TERM
dwd = reshape(repmat(wd',9,1),9*(L-2),1).^(1/2);
eye9LL = [zeros(9*(L-2),9), eye(9*(L-2))];
eye9LR = [eye(9*(L-2)), zeros(9*(L-2),9)];

Ad = diag(dwd)*(eye9LL - eye9LR);

% ALL TOGETHER (VERIFIED)
Q = zeros(size(Ag,1) + size(Ad,1) + size(Av,1),d+1);
% measurements
Q(1:size(Ag,1),1:(1+9*L+9*(L-1)+3*L)) = ...
    [Ag, Ar, zeros(size(Ag,1),9*(L-1)), As];
% rotation rate
Q((size(Ag,1) + 1):(size(Ag,1)+size(Ad,1)),...
  (1+9*L+1):(1+9*L+9*(L-1))) = ...
    Ad;
% velocity
Q((size(Ag,1)+size(Ad,1) + 1):(size(Ag,1)+size(Ad,1) + size(Av,1)),...
  (1+9*L+9*(L-1)+3*L+1):(1+21*L-9+3*(L-1))) = ...
    Av;
prob_obj = [1;x]'*(Q'*Q)*[1;x];
% c = Cr*r - Cs*s + gbar;

%% Define constraints
if problem.regen_sdp
% regenerate constraints

% EQUALITY
h = [];
% SO(3) constraints
for l = 1:L
    c1 = so3_constraints( R(ib3(l),:));
    if (l < L)
        c2 = so3_constraints(dR(ib3(l),:));
        h = [h; c1; c2];
    else
        h = [h; c1];
    end
end

% R(t+1) = R(t) dR(t) constraint
for l = 2:L
    h = [h; reshape(R(ib3(l),:) - R(ib3(l-1),:)*dR(ib3(l-1),:),9,1)];
end

% sh(l) = s(l-1) + v(l-1)*dt constraint
for l = 2:L
    % dR version
    h = [h; dR(ib3(l-1),:)*s(ib3(l)) - s(ib3(l-1)) - v(ib3(l-1))*dt];
    h = [h; s(ib3(l)) - dR(ib3(l-1),:)'*s(ib3(l-1)) - dR(ib3(l-1),:)'*v(ib3(l-1))*dt];
end

% constraint on v(t1) as a function of v(t2)
% TODO: this may help solve time?
h = [h; dR(ib3(1),:)'*v(ib3(1))*dt - dR(ib3(2),:)*s(ib3(3)) + dR(ib3(1),:)'*s(ib3(1)) + v(ib3(2))*dt];

% INEQUALITY
% p,s in range for just first time (p'*p<=pBoundSq)
% TODO: enforce in range for ALL time steps - LC: agreed
% TODO: add in sh
pBoundSq = pBound^2;
g_s_first = pBoundSq*L - s(ib3(1))'*s(ib3(1));

% v bound (v'*v<=vBoundSq)
vBoundSq = vBound^2;
g_v = vBoundSq*L - v'*v;
% g_v = [];
% for l = 1:L-1
%     g_v = [g_v; vBoundSq - v(ib3(l))'*v(ib3(l))];
% end

% c bound (0<=c<=1) (if using, uncomment c def)
% cBoundSq = 1.0; % should just be 1
% g_c = [cBoundSq - c'*c;c]; % includes measurements

% g = [g_c];
% g = [g_c;g_v];
% g = [g_s_first; g_v; g_c]; % only use if regenerating each time.
g = [g_s_first; g_v];

% save("data/constraints.mat","g","h");
% else
    % this runs much faster
    % load data/sdp_data.mat g h
end


%% Complete problem definition
problem.vars = x;
problem.objective = prob_obj;
if ~problem.regen_sdp
    load sdpdata.mat sdpdata
else
    sdpdata = [];
    problem.equality = h; % equality
    problem.inequality = g; % inequality
end

%% Relax!
kappa = 1; % relaxation order
% [SDP,info] = dense_sdp_relax(problem,kappa);
[SDP,info,sdpdata] = fast_sdp_relax(problem,kappa,sdpdata);

%% Solve using MOSEK
tic
prob = convert_sedumi2mosek(SDP.sedumi.At,...
                            SDP.sedumi.b,...
                            SDP.sedumi.c,...
                            SDP.sedumi.K);
% addpath(genpath(mosekpath))
[~,res] = mosekopt('minimize info echo(0)',prob);
[Xopt,yopt,Sopt,obj] = recover_mosek_sol_blk(res,SDP.blk);
% rmpath(genpath(mosekpath))
soln.solvetime = toc;

% figure; bar(eig(Xopt{1})); % if rank = 1, then relaxation is exact/tight

%% Compare to ground truth
% clip to first eigenvalue
[eigvecs, ~] = eig(Xopt{1});
vecmax = eigvecs(:,end);
% re-normalize so first element is 1
vecmax_normed = vecmax / vecmax(1);
x_est = vecmax_normed(2:end);

% Project to SO(3) and extract results
rs = x_est(1:(9*L));
Rs = projectRList(rs);
drs = x_est((9*L+1):(18*L-9));
dRs = projectRList(drs);

% s_est = reshape(full(dmsubs(s,x,x_est)),[3,1,L]);
s_est = reshape(x_est((18*L-9+1):(18*L-9+3*L)),[3,1,L]);
% v_est = reshape(full(dmsubs(v,x,x_est)),[3,1,L-1]);
v_est = reshape(x_est((18*L-9+3*L+1):(21*L-9+3*L-3)),[3,1,L-1]);
% c_est = full(dmsubs(c,x,x_est));
c_est = Cr*reshape(Rs,9*L,1,1) - Cs*reshape(s_est,3*L,1,1) + gbar;
% c = Cr*r - Cs*s + gbar;

% estimate p from s
p_est = zeros(3,1,L);
for l = 1:L
    p_est(:,:,l) = Rs(:,:,l)*s_est(:,:,l);
end

% suboptimality gap
% x_proj = [1.0];
% for l = 1:L
%     r_temp = reshape(Rs(:,:,l),9,1);
%     x_proj = [x_proj; r_temp];
% end
% for l = 1:L-1
%     r_temp = reshape(dRs(:,:,l),9,1);
%     x_proj = [x_proj; r_temp];
% end
% x_proj = [x_proj; reshape(s_est,[3*L,1,1]); reshape(v_est,[3*L-3,1,1])];

x_proj = [1.0; reshape(Rs,[9*L,1])];
% dR from R
for l = 1:L-1
    dr = reshape(Rs(:,:,l)'*Rs(:,:,l+1),[9,1]);
    x_proj = [x_proj; dr];
end
x_proj = [x_proj; reshape(s_est,[3*L,1,1])];
for l = 1:L-1
    v = (1/dt)*(Rs(:,:,l)'*Rs(:,:,l+1) * s_est(:,:,l+1) - s_est(:,:,l));
    x_proj = [x_proj; v];
end

% compute gap
% obj_est = dmsubs(prob_obj,x,x_proj); % slow
obj_est = x_proj'*(Q'*Q)*x_proj;
gap = (obj_est - obj(1)) / obj_est;

% compute residuals
residuals = zeros(N, L);
for i = 1:N
    for l = 1:L
        residue = Rs(:,:,l)'*y(ib3(i),l) - B(ib3(i),:)*c_est - s_est(:,:,l);
        residuals(i,l) = (residue'*residue);
    end
end
% residuals = residuals ./ (problem.noiseBoundSq);

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
soln.v_est = v_est;
soln.s_est = s_est;

soln.R_est = Rs;
soln.dR_est = dRs;

soln.gap = gap;
soln.x_proj = x_proj;
soln.obj_est = obj_est;

soln.residuals = residuals;

%% Save SDP data
if problem.regen_sdp
    save("sdpdata.mat","sdpdata");
end

end