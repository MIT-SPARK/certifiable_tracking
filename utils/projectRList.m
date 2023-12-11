function Rs = projectRList(r)
% project list of Rs into SO(3) conveniently

N = length(r)/9;
temp_R  = reshape(r ,3, 3*N)';
Rs = zeros(3,3,N);
for i = 1:N
    Rs(:,:,i) = project2SO3(temp_R(ib3(i),:)');
end
end