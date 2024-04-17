function soln = solve_tracking_body(problem)
% Solves const. vel. (body frame) optimization exactly via SDP
%   Assume the *body frame* velocity is constant and object is spinning.
%   Result: spiral trajectory.
%   Analytically remove velocity & shape. SDP variables are
%   * rotated position (s)
%   * body velocity (v)
%   * rotation (R)
%   * rotation change (dR)
%
%   Sparse relaxation: velocity has its own semidefinite variable
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
W = problem.covar_measure.^(-1); % N x L matrix of w_il
lambda = problem.lambda; % scalar
wv = problem.covar_velocity.^(-1); % L-2 vector
wd = problem.kappa_rotrate;  % L-1 vector

% pBound = problem.translationBound;
vBound = problem.velocityBound;

% check lambda
if ((K >= 3*N) && (lambda == 0.0))
    error("lambda must be positive when there are more shapes (K) " + ...
        "than measurement points (3*N)");
end

%% Define objective
% optimization vector
d = 9*(2*L - 1) + 3*L; % 2L - 1 rotations, L rotated positions
% 2L-1 rotations: L rotations, L-1 delta rotations
x = msspoly('x',d);

% pull out individual variables
r  = x(1:(9*L));
dr = x((9*L + 1):(9*L + 9*(L-1)));
s  = x((18*L - 9 + 1):(18*L - 9 + 3*L));

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
Q = zeros(size(Ag,1) + size(Ad,1),d+1);
% measurements
Q(1:size(Ag,1),1:(1+9*L+9*(L-1)+3*L)) = ...
    [Ag, Ar, zeros(size(Ag,1),9*(L-1)), As];
% rotation rate
Q((size(Ag,1) + 1):(size(Ag,1)+size(Ad,1)),...
  (1+9*L+1):(1+9*L+9*(L-1))) = ...
    Ad;
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

% INEQUALITY
% p,s in range for just first time (p'*p<=pBoundSq)
% pBoundSq = pBound^2;
% g_s_first = pBoundSq*L - s(ib3(1))'*s(ib3(1));
% g_s = [];
% for l = 1:L
%     g_s = [g_s; pBoundSq - s(ib3(l))'*s(ib3(l))];
% end

if isfield(problem,"usecBound")
    if problem.usecBound
        % c bound (0<=c<=1)
        c = Cr*r - Cs*s + gbar;
        g_c = [c; 1-c; c.^2; 1-c.^2];
        g = g_c;
    else
        g = [];
    end
else
g = []; % no constraints works just as well
end

end


%% Complete problem definition
problem.vars = x;
problem.objective = prob_obj;
if ~problem.regen_sdp
    if isfield(problem,"sdp_filename")
        load(problem.sdp_filename +".mat","sdpdata");
    else
        pid = string(feature("getpid"));
        load("sdpdata"+ pid +".mat","sdpdata");
    end
else
    sdpdata = [];
    problem.equality = h; % equality
    problem.inequality = g; % inequality
end

%% Relax!
kappa = 1; % relaxation order
[SDP,info,sdpdata] = fast_sdp_relax(problem,kappa,sdpdata);

%% Solve using MOSEK
prob = convert_sedumi2mosek(SDP.sedumi.At,...
                            SDP.sedumi.b,...
                            SDP.sedumi.c,...
                            SDP.sedumi.K);
% prob = add_v(prob,Av,L,dt,vBound);
prob = add_quad_v_constraints(prob,Av,L,dt,vBound);
tic
[~,res] = mosekopt('minimize info echo(0)',prob);
soln.solvetime = toc;
% blk = {SDP.blk{1}, SDP.blk{2}; SDP.blk{1}, [1+3*(L-1)]};
blk = SDP.blk;
[Xopt,yopt,Sopt,obj] = recover_mosek_sol_blk(res,blk);

%% Compare to ground truth
% clip to first eigenvalue
[eigvecs, ~] = eig(Xopt{1});
vecmax = eigvecs(:,end);
% re-normalize so first element is 1
vecmax_normed = vecmax / vecmax(1);
x_est = vecmax_normed(2:end);
% x_est = Xopt{1};
% x_est = x_est(2:end,1);

% Project to SO(3) and extract results
rs = x_est(1:(9*L));
Rs = projectRList(rs);
drs = x_est((9*L+1):(18*L-9));
dRs = projectRList(drs);

% s_est = reshape(full(dmsubs(s,x,x_est)),[3,1,L]);
s_est = reshape(x_est((18*L-9+1):(18*L-9+3*L)),[3,1,L]);
% v_est = reshape(full(dmsubs(v,x,x_est)),[3,1,L-1]);
v_est = reshape(res.sol.itr.xx(2:end),[3,1,L-1]);
% c_est = full(dmsubs(c,x,x_est));
c_est = Cr*reshape(Rs,9*L,1,1) - Cs*reshape(s_est,3*L,1,1) + gbar;
% c = Cr*r - Cs*s + gbar;

% estimate p from s
p_est = zeros(3,1,L);
for l = 1:L
    p_est(:,:,l) = Rs(:,:,l)*s_est(:,:,l);
end

x_proj = [1.0; reshape(Rs,[9*L,1])];
% dR from R
for l = 1:L-1
    dr = reshape(Rs(:,:,l)'*Rs(:,:,l+1),[9,1]);
    x_proj = [x_proj; dr];
end
x_proj = [x_proj; reshape(s_est,[3*L,1,1])];
v_est_corrected = zeros(3*(L-1),1);
for l = 1:L-1
    v = (1/dt)*(Rs(:,:,l)'*Rs(:,:,l+1) * s_est(:,:,l+1) - s_est(:,:,l));
    v_est_corrected(ib3(l)) = v;
end

% compute gap
% obj_est = dmsubs(prob_obj,x,x_proj); % slow
obj_est = x_proj'*(Q'*Q)*x_proj + v_est_corrected'*(Av'*Av)*v_est_corrected;
gap = (obj_est - obj(2)) / obj(2);
gap_stable = (obj_est - obj(2)) / (obj(2)+1);

% compute residuals
% residuals = zeros(N, L);
% for i = 1:N
%     for l = 1:L
%         residue = Rs(:,:,l)'*y(ib3(i),l) - B(ib3(i),:)*c_est - s_est(:,:,l);
%         residuals(i,l) = (residue'*residue);
%     end
% end
% corrected c
c_est_corrected = c_est;
c_est_corrected(c_est < 0) = 0;
c_est_corrected(c_est > 1) = 1;
% compute residuals with corrected c
residuals = zeros(N, L);
for i = 1:N
    for l = 1:L
        residue = Rs(:,:,l)'*y(ib3(i),l) - B(ib3(i),:)*c_est_corrected - s_est(:,:,l);
        residuals(i,l) = (residue'*residue);
    end
end

%% Alt reprojection: generative model
% % use v, dR to gen p, R from p0, R0
% s2_est = zeros(3,1,L);
% p2_est = zeros(3,1,L);
% R2_est = zeros(3,3,L);
% 
% % initialize p, s, R
% s2_est(:,:,1) = s_est(:,:,1);
% p2_est(:,:,1) = p_est(:,:,1); % this came from s
% R2_est(:,:,1) = Rs(:,:,1);
% 
% % estimate p, s, R from p0, R0, v, dR
% for l = 1:L-1
%     s2_est(:,:,l+1) = dRs(:,:,l)'*(s2_est(:,:,l) + v_est(:,:,l)*dt);
%     p2_est(:,:,l+1) = p2_est(:,:,l) + R2_est(:,:,l)*v_est(:,:,l)*dt;
%     R2_est(:,:,l+1) = R2_est(:,:,l) * dRs(:,:,l);
% end
% 
% % gap
% x_proj2 = [1.0; reshape(R2_est,[9*L,1]); reshape(dRs,[9*(L-1),1]);
%                 reshape(s2_est,[3*L,1])];
% v2 = reshape(v_est,[3*(L-1),1]);
% obj_est2 = x_proj2'*(Q'*Q)*x_proj2 + v2'*(Av'*Av)*v2;
% gap2 = (obj_est2 - obj(2)) / obj(2);
% 
% % gt gap
% x_gt = [1; problem.x_gt];
% v_gt = reshape(problem.v_gt,[3*L-3,1]);
% obj_gt = x_gt'*(Q'*Q)*x_gt + v_gt'*(Av'*Av)*v_gt;
% gap_gt = (obj_gt - obj(2)) / obj(2);
% 
% % save
% soln.s2_est = s2_est;
% soln.p2_est = p2_est;
% soln.R2_est = R2_est;
% soln.dR2_est = dRs;
% soln.obj_est2 = obj_est2;
% soln.gap2 = gap2;
% soln.gap_gt = gap_gt;

% % Alt reprojection: Lasserre
% 1) chol
% [U,S,~] = svd(Xopt{1}); % U = V
% S(S < 1e-6) = 0;
% V = U*sqrt(S);
% V = V(:,1:rank(S));
% % 2) column reduce
% U = rref(V')';


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
soln.v_est_corrected = reshape(v_est_corrected,[3,1,L-1]);
soln.s_est = s_est;

soln.R_est = Rs;
soln.dR_est = dRs;

soln.gap = gap;
soln.gap_stable = gap_stable;
soln.x_proj = x_proj;
soln.obj_est = obj_est;

soln.residuals = residuals;

%% Save SDP data
if problem.regen_sdp
    if isfield(problem,"sdp_filename")
        save(problem.sdp_filename +".mat","sdpdata");
    else
        pid = string(feature("getpid"));
        save("sdpdata"+ pid +".mat","sdpdata");
    end
end

end

%% Add v constraints
function prob = add_v(prob,Av,L,dt, vBound)
    % adds 3*(L-1) x 3*(L-1) SDP of velocity variables
    % and associated constraints

    % create a new PSD variable for velocity
    prob.bardim = [prob.bardim; 1+3*(L-1)];

    % objective
    [r2,c2,v2] = find(tril(Av'*Av));
    prob.barc.subj = [prob.barc.subj,...
                      repmat(2,1,length(v2))];
    prob.barc.subk = [prob.barc.subk, r2'+1];
    prob.barc.subl = [prob.barc.subl, c2'+1];
    prob.barc.val = [prob.barc.val, v2'];

    % constrain 1 = 1
    prob.blc = [prob.blc; 1];
    prob.buc = [prob.buc; 1];
    prob.bara.subi = [prob.bara.subi, prob.bara.subi(end)+1];
    prob.bara.subj = [prob.bara.subj, 2];
    prob.bara.subk = [prob.bara.subk, 1];
    prob.bara.subl = [prob.bara.subl, 1];
    prob.bara.val = [prob.bara.val, 1];

    % add constraints between SDP vars
    slices_d = 1+((9*L+1):(9*L+9*(L-1)));
    slices_s = 1+((18*L-9+1):(18*L-9+3*L));
    slices_v = 1+(1:3*(L-1));
    for l = 1:L-1
        % A1: dR(ib3(l-1),:)*s(ib3(l)) - s(ib3(l-1)) 
        r1 = [repmat(slices_s(ib3(l+1)),[1,3]),slices_s(ib3(l))];
        c1 = [slices_d([1,4,7]+9*(l-1)),slices_d([2,5,8]+9*(l-1)),slices_d([3,6,9]+9*(l-1)),ones(1,3)];
        v1 = [1*ones(1,9),-1*ones(1,3)];
        % A2: -v(ib3(l-1))*dt
        r2 = [slices_v(ib3(l))];
        c2 = [1*ones(1,3)];
        v2 = [-dt*ones(1,3)];
        
        prob.blc = [prob.blc; zeros(3,1)];
        prob.buc = [prob.buc; zeros(3,1)];
        prob.bara.subi = [prob.bara.subi, ...
                            repelem(prob.bara.subi(end) + (1:3),3), ...
                            repmat(prob.bara.subi(end) + (1:3),[1,2])];
        prob.bara.subj = [prob.bara.subj, 1*ones(1,12), 2*ones(1,3)];
        prob.bara.subk = [prob.bara.subk, r1, r2];
        prob.bara.subl = [prob.bara.subl, c1, c2];
        prob.bara.val  = [prob.bara.val , v1, v2];
    end

    % add that v.^2 < vbound^2
    vBoundSq = vBound^2;
    prob.blc = [prob.blc; 0];
    prob.buc = [prob.buc; vBoundSq];
    prob.bara.subi = [prob.bara.subi, prob.bara.subi(end)+ones(1,3)];
    prob.bara.subj = [prob.bara.subj, 2,2,2];
    prob.bara.subk = [prob.bara.subk, 2,3,4];
    prob.bara.subl = [prob.bara.subl, 2,3,4];
    prob.bara.val = [prob.bara.val, 1,1,1];

    prob.a = sparse([], [], [], length(prob.blc), 0); 
end

%% add v constraints, but quadratic
function prob = add_quad_v_constraints(prob, Av, L, dt, vBound)
    % add new scalar (non-SDP) variables:
    % sh (dim 1) and v (dim 3*L-3)
    % in objective: sh replaces v'*(Av'*Av)*v
    % in constraints: add constraint v'*(Av'*Av)*v <= sh
    % result: minimizing sh -> sh = (v'*(Av'*Av)*v)
    [~, res] = mosekopt('symbcon');
    symbcon = res.symbcon;

    % objective: add sh
    prob.c = [1, zeros(1,3*(L-1))];
    % define 1 + 3*(L-1) scalar variables with 3*(L-1) new constraints
    % prob.a = sparse([], [], [], length(prob.blc) + 3*(L-1), 1 + 3*(L-1));
    prob.a = sparse([], [], [], length(prob.blc), 1 + 3*(L-1));
    % bounds: [-vBound, vBound]
    prob.blx = [0, -vBound*ones(1,3*(L-1))];
    prob.bux = [Inf,  vBound*ones(1,3*(L-1))];
    
    % conic constraint: v'*(Av'*Av)*v <= sh
    % (sh, 1, (sqrt(2)*Av)*v) in Q^(2 + 3*(L-1))
    k = size(Av,1);
    fQ = sparse([1/2, zeros(1,3*L-3);zeros(1,1+3*L-3); zeros(k,1) sparse(Av)]); % equiv to sqrt(2)*Av
    gQ = [0, 1, zeros(1,k)]';
    accQ = [symbcon.MSK_DOMAIN_RQUADRATIC_CONE, 2 + k];

    % vbound constraint: (1/2, vBound, v) in Qr^(1 + 3)
    fV = sparse([zeros(2,1+3*L-3); zeros(3, 1), speye(3), zeros(3,3*L-6)]);
    gV = [1/2, vBound^2, zeros(1,3)]';
    accV = [symbcon.MSK_DOMAIN_RQUADRATIC_CONE, 2 + 3];
    % fV = []; gV = []; accV = [];

    prob.f = [fQ; fV];
    prob.g = [gQ; gV];
    prob.accs = [accQ, accV];

    % add constraints between v & semidefinite variable
    slices_d = 1+((9*L+1):(9*L+9*(L-1)));
    slices_s = 1+((18*L-9+1):(18*L-9+3*L));
    slices_v = 1+(1:3*(L-1));
    for l = 1:L-1
        % A1: dR(ib3(l-1),:)*s(ib3(l)) - s(ib3(l-1))
        r1 = [repmat(slices_s(ib3(l+1)),[1,3]),slices_s(ib3(l))];
        c1 = [slices_d([1,4,7]+9*(l-1)),slices_d([2,5,8]+9*(l-1)),slices_d([3,6,9]+9*(l-1)),ones(1,3)];
        v1 = [1*ones(1,9),-1*ones(1,3)];
        % A2: -v(ib3(l-1))*dt
        c2 = [slices_v(ib3(l))];
        r2 = [prob.bara.subi(end) + (1:3)];
        v2 = [-2*dt*speye(3)];
        prob.a(r2, c2) = v2;

        prob.blc = [prob.blc; zeros(3,1)];
        prob.buc = [prob.buc; zeros(3,1)];
        prob.bara.subi = [prob.bara.subi, ...
                            repelem(prob.bara.subi(end) + (1:3),3), ...
                            prob.bara.subi(end) + (1:3)];
        prob.bara.subj = [prob.bara.subj, 1*ones(1,12)];
        prob.bara.subk = [prob.bara.subk, r1];
        prob.bara.subl = [prob.bara.subl, c1];
        prob.bara.val  = [prob.bara.val , v1];
    end

    % add constraints between v and semidefinite variable (acc version)
    % slices_d = 1+((9*L+1):(9*L+9*(L-1)));
    % slices_s = 1+((18*L-9+1):(18*L-9+3*L));
    % slices_v = 1+(1:3*(L-1));
    % for l = 1:L-1
    %     % A1: dR(ib3(l-1),:)*s(ib3(l)) - s(ib3(l-1))
    %     r1 = [repmat(slices_s(ib3(l+1)),[1,3]),slices_s(ib3(l))];
    %     c1 = [slices_d([1,4,7]+9*(l-1)),slices_d([2,5,8]+9*(l-1)),slices_d([3,6,9]+9*(l-1)),ones(1,3)];
    %     v1 = [1*ones(1,9),-1*ones(1,3)];
    %     % A2: -v(ib3(l-1))*dt
    %     c2 = [slices_v(ib3(l))];
    %     if l == 1
    %         r2 = [size(prob.f, 1) + (1:3)];
    %     else
    %         r2 = [prob.barf.subi(end) + (1:3)];
    %     end
    %     v2 = [-dt*eye(3)];
    % 
    %     prob.g = [prob.g; zeros(3,1)];
    %     prob.accs = [prob.accs, symbcon.MSK_DOMAIN_RZERO, 3];
    % 
    %     if l == 1
    %         prob.barf.subi = [size(prob.f,1) + repelem((1:3),3), size(prob.f,1) + (1:3)];
    %         prob.barf.subj = [ones(1,12)];
    %         prob.barf.subk = [r1];
    %         prob.barf.subl = [c1];
    %         prob.barf.val  = [v1];
    %     else
    %         prob.barf.subi = [prob.barf.subi, ...
    %                             repelem(prob.barf.subi(end) + (1:3),3), ...
    %                             prob.barf.subi(end) + (1:3)];
    %         prob.barf.subj = [prob.barf.subj, 1*ones(1,12)];
    %         prob.barf.subk = [prob.barf.subk, r1];
    %         prob.barf.subl = [prob.barf.subl, c1];
    %         prob.barf.val  = [prob.barf.val , v1];
    %     end
    % 
    %     prob.f(r2, c2) = v2;
    % end

end