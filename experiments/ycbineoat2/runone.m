function runone(params)
%% Run just ONE BATCH

%% Experiment Settings
videos = ["cracker_box_reorient", "cracker_box_yalehand0", ...
          "sugar_box1", "sugar_box_yalehand0", ...
          "mustard0", "mustard_easy_00_02", ...
          "bleach0", "bleach_hard_00_03_chaitanya", ...
          "tomato_soup_can_yalehand0"];

% parameters to change
video = videos(params.videoNumber);
params.video = video;
maxL = params.maxL;

% for pruning
noiseBound_GRAPH = params.noiseBound_GRAPH;
% for GNC
noiseBound_GNC = params.noiseBound_GNC;
% for solver
velocityBound = params.velocityBound;
translationBound = params.translationBound;
covar_measure_base = params.covar_measure_base;
covar_velocity_base = params.covar_velocity_base;
covar_rotrate_base = params.covar_rotrate_base;

savename = params.savename;
jsondir = "../datasets/YCBInEOAT/";

%% Load pruned problems
load(savename,"problems","gt","teaser")
numProblemsToSolve = length(problems);

%% Run single batch
batch = params.batch;
solns = {};

curproblem = problems{batch};
curproblem.sdp_filename = "sdpdata" + curproblem.L;
curproblem.regen_sdp = true;

curproblem.noiseBound_GNC = noiseBound_GNC;
curproblem.velocityBound = velocityBound;
curproblem.translationBound = translationBound;
curproblem.covar_measure_base = covar_measure_base;
curproblem.covar_velocity_base = covar_velocity_base;
curproblem.covar_rotrate_base = covar_rotrate_base;

if (maxL ~= curproblem.L)
    curproblem.L = maxL;
    curproblem.y = curproblem.y(:,(end-maxL+1):end);
    problems{batch} = curproblem;
end

soln = solveBatch(curproblem);
solns{1} = soln;

% report bad results
if (soln.gap > 1e-2)
    fprintf("Batch %d failed: %.4e\n",batch,soln.gap)
elseif (isnan(soln.gap))
    fprintf("Batch %d failed: NaN\n",batch)
end

%% Save to mat
save(savename,"params","problems","solns","gt","teaser")

%% Save to JSON
L = curproblem.L;
T_est = repmat(eye(4),[1,1,L]);
for l = 1:L
    T_est(1:3,1:3,l) = soln.R(:,:,l);
    T_est(1:3,4,l) = soln.p(:,:,l)*1000.0;
end

fid = fopen(curproblem.json); 
raw = fread(fid,inf); 
str = char(raw'); 
fclose(fid);
data = jsondecode(str);
data = data((batch+2-L+1):(batch+2));

for l = 1:L
    data(l).cast_pose = T_est(:,:,l);
end

cocoString = jsonencode(data, "PrettyPrint",true);
fid = fopen("test.json", 'w');
fprintf(fid, '%s', cocoString);
fclose(fid);

end