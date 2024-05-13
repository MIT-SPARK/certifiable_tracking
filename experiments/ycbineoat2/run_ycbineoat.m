%% Base parameters
params = struct();
params.videoNumber = 3; %***********
params.skipPruning = false;
params.interp = false;

params.savename = "ycbineoat_" + string(params.videoNumber);
params.maxL = 8;

% for outliers
params.noiseBound_GRAPH = 0.01/2;
params.noiseBound_GNC = 0.01/2;
% for solver
params.velocityBound = 1.5;
params.translationBound = 2.5;
params.covar_measure_base = 0.01;
params.covar_velocity_base = 10*0.01;
params.covar_rotrate_base = 10*0.01;

ParamList = params;

%% Set param sweep
% ratios = [0.1, 0.01];
% for r = ratios
%     params2 = params;
%     params2.savename = strrep(params2.savename + "_r" + string(r),'.','');
%     params2.covar_velocity_base = r*params2.covar_velocity_base;
%     params2.covar_rotrate_base = r*params2.covar_rotrate_base;
%     ParamList = [ParamList; params2];
% end
% 
% noisebounds = [0.02, 0.05, 0.1];
% for n = noisebounds
%     params2 = params;
%     params2.savename = strrep(params2.savename + "_ngra" + string(n),'.','');
%     params.noiseBound_GRAPH = n;
%     ParamList = [ParamList; params2];
% end
% 
% noisebounds = [0.02, 0.05, 0.1];
% for n = noisebounds
%     params2 = params;
%     params2.savename = strrep(params2.savename + "_ngnc" + string(n),'.','');
%     params.noiseBound_GNC = n;
%     ParamList = [ParamList; params2];
% end
% videos = [4,5,6];
% for v = videos
%     params2 = params;
%     params2.videoNumber = v;
%     params2.savename = "ycbineoat_" + string(v);
%     ParamList = [ParamList; params2];
% end

%% Run!
for p = 1:length(ParamList)
    par = ParamList(p);
    ycbineoat(par)
end