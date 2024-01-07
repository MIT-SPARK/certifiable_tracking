function A_list = time_consistency(problem, C_list, NOISE_BOUND)
% For each edge, compute the percentage of edges
% that are rigid-body consistent.
% For example:
%   Let p, q be keypoints seen at t1, t2, t3, t4
%       with shape-consistent edges @ t1, t3, t4
%   Check rigid consistency (~):
%       t1 ~ t3, t4 -> edge p-q(t1) has score 1
%       t1 ~ none -> edge p-q(t1) has score 0.33
%       t1 ~ t3 -> edge p-q(t1) has score 0.67
% 
% Lorenzo Shaikewitz for SPARK Lab

% setup
L = problem.L;
N = problem.N_VAR;
y = problem.y; % 3*N x L

A_list = zeros(N,N,L); % sym
A_counter = ones(N,N,L); % sym

% for each time compare with every other time
% TODO: FIX LOOP!
for l1 = 1:(L-1)
    for l2 = (l1+1):L
        % work through each edge
        for i1 = 1:(N-1)
            for i2 = (i1+1):N
                % only count consistent edges
                if (C_list(i1,i2,l1) + C_list(i1,i2,l2)) == 0
                    continue; end
                % count if consistent
                A_counter(i1,i2,l1) = A_counter(i1,i2,l1) + 1;
                A_counter(i1,i2,l2) = A_counter(i1,i2,l2) + 1;

                % check rigid consistency
                d1 = norm(y(ib3(i1),l1) - y(ib3(i2),l2));
                d2 = norm(y(ib3(i1),l2) - y(ib3(i2),l2));
                if (abs(d1 - d2) < NOISE_BOUND)
                    A_list(i1,i2,l1)  = A_list(i1,i2,l1) + 1;
                    A_list(i1,i2,l2)  = A_list(i1,i2,l2) + 1;
                end
            end
        end
    end
end

% make symetric
A_list = 0.5*(A_list + permute(A_list,[2,1,3]));
A_counter = 0.5*(A_counter + permute(A_counter,[2,1,3]));

% normalize
A_list = A_list ./ A_counter;

end