%% Test script for ROBIN
% Lorenzo Shaikewitz for SPARK Lab

addpath("../../ROBIN/build2/matlab");

%% Create graph

g = zeros(20);

for i = 1:10
    % add vertex
    g(i,i) = 1;
    g(i+10,i+10) = 1;
    % add edge
    g(i,i+10) = 1;
    g(i+10,i) = 1;
end

% add additional edge
g(1,12) = 1;
g(12,1) = 1;
g(11,12) = 1;
g(12,11) = 1;

%% Find inlier structures
[max_core_indicies,~] = findInlierStructure_mex(g,0);
[max_clique_indicies,~] = findInlierStructure_mex(g,1);

max_core_indicies
