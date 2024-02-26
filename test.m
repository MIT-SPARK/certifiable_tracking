%% Test script to check properties of solution
% use solve_tracking_body_dense
%
% Lorenzo Shaikewitz for SPARK Lab

%% 1) The velocity is in the span of the eigenvalues of Xopt
[v, lam] = eig(soln.raw.Xopt{1});
A = v(:,end-3:end);
xv = soln.x_proj;
c = A \ xv;

%% 2) Check if the matricies B_k are linearly independent
[v, lam] = eig(soln.raw.Xopt{1});
m = 3*(L-1);
n = 1 + 9*(2*L - 1) + 3*L + 3*(L-1);
B_k = zeros(n,4,m);

slices_d = 1+((9*L+1):(9*L+9*(L-1)));
slices_s = 1+((18*L-9+1):(18*L-9+3*L));
slices_v = 1+((21*L-9+1):(21*L-9+3*(L-1)));
dt = problem.dt;
for k = 0:m-1
    % only check Aks that involve v AND 1
    % dR(ib3(l-1),:)*s(ib3(l)) - s(ib3(l-1)) - v(ib3(l-1))*dt
    l = floor(k/3)+2;
    idx_l = ib3(l);
    idx_lm1 = ib3(l-1);
    if (mod(k,3) == 1)
        % dR(l-1)_1 * s(l)_1 + dR(l-1)_1  * s(l)_2 + dR(l-1)_1 * s(l)_3 - ...
        % s(l-1)_1 - v(l-1)_1*dt
        Ak = sparse(...
            [slices_d(idx_lm1(1)), slices_d(idx_lm1(2)), slices_d(idx_lm1(3)), slices_s(idx_lm1(1)), slices_v(idx_lm1(1))], ...
            [slices_s(idx_l(1))*ones(1,3), 1, 1],...
            [1, 1, 1, -1, -dt]*0.5,n,n);
        Ak = Ak + Ak';
    elseif (mod(k,3) == 2)
        Ak = sparse(...
            [slices_d(idx_lm1(1)), slices_d(idx_lm1(2)), slices_d(idx_lm1(3)), slices_s(idx_lm1(2)), slices_v(idx_lm1(2))], ...
            [slices_s(idx_l(2))*ones(1,3), 1, 1],...
            [1, 1, 1, -1, -dt]*0.5,n,n);
        Ak = Ak + Ak';
    else
        Ak = sparse(...
            [slices_d(idx_lm1(1)), slices_d(idx_lm1(2)), slices_d(idx_lm1(3)), slices_s(idx_lm1(3)), slices_v(idx_lm1(3))], ...
            [slices_s(idx_l(3))*ones(1,3), 1, 1],...
            [1, 1, 1, -1, -dt]*0.5,n,n);
        Ak = Ak + Ak';
    end
    Q1 = v(:,end:-1:end-3);
    Q2 = v(:,end-4:-1:1);
    B_k(:,:,k+1) = [Q1'*Ak*Q1; Q2'*Ak*Q1];
end
% add the Ak involving 1
Ak = sparse(1,1,1,n,n);
B_k(:,:,end+1) = [Q1'*Ak*Q1; Q2'*Ak*Q1];
m = m+1;

% check if they are linearly independent
b = zeros(n*4,m);
for k = 1:m
    b(:,k) = reshape(B_k(:,:,k),[n*4,1,1]);
end
B = squeeze(b);
s = svd(B);
% rank(B) = # nonzero elements of s