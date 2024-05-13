%% Step 4: Visualize Results
% Load solved problems and display results
% Run Experiment Settings section before this
%
% Lorenzo Shaikewitz for SPARK Lab

videoNumber = "2_baseline";
savename = "ycbineoat_" + videoNumber;
load(savename);

%% Reformat solutions
est = getEstimates(problems, solns);
% est.p = est.p - [0.0016; 0.0076; -0.0072]; % offset for cracker

% plot x,y,z
plotSolns(gt, est);

% compute scores
adds = false; % takes a while
thresh = 0.1;
figure; hold on;
tab_new = computeScores(problems{1}, gt, teaser, est, params.video, adds, thresh)

% save to json
% save2json(est, problems{1})

% get errors
est = getErrors(gt, est);

%% Save everything together

%% Helper function: reformat solutions
function est = getEstimates(problems, solns)
% convert solns (list of structs) to compact form
est.p = zeros(3,1,length(solns))*NaN;
est.R = zeros(3,3,length(solns))*NaN;
est.c = zeros(length(solns),1)*NaN;

for j = 1:length(solns)
    problem = problems{j};
    soln = solns{j};
    L_cur = problem.L;

    idx = problem.startIdx:(problem.startIdx + L_cur-1);
    if j == 1
        % update estimates for times 1, 2, 3
        est.p(:,:,idx) = soln.p;
        est.R(:,:,idx) = soln.R;
        [~, est.c(idx)] = max(soln.c);
    else
        % update estimate for current time only
        est.p(:,:,idx(end)) = soln.p(:,:,end);
        est.R(:,:,idx(end)) = soln.R(:,:,end);
        [~, est.c(idx(end))] = max(soln.c);
    end
end

end

%% Helper function: plot solutions
function plotSolns(gt, est)
% plot x, y, z
p_gt = reshape(gt.p,[3,size(gt.p,3),1]);
p_est = reshape(est.p,[3,size(est.p,3),1]);

L = length(est.R);
w_est = zeros(3,L);
w_gt = zeros(3,L);
for l = 1:length(est.R)
    axang = rotm2axang(est.R(:,:,l));
    w_est(:,l) = axang(1:3)*axang(4);
    axang = rotm2axang(gt.R(:,:,l));
    w_gt(:,l) = axang(1:3)*axang(4);
end

figure
t=tiledlayout(3,2);
nexttile
plot(p_gt(1,:)); hold on; plot(p_est(1,:))
ylabel("x")
title("Positions")
nexttile
plot(w_est(1,:)); hold on; plot(w_gt(1,:),'k');
title("Rotations")

nexttile
plot(p_gt(2,:)); hold on; plot(p_est(2,:))
ylabel("y")
nexttile
plot(w_est(2,:)); hold on; plot(w_gt(2,:),'k');

nexttile
plot(p_gt(3,:)); hold on; plot(p_est(3,:))
ylabel("z")
nexttile
plot(w_est(3,:)); hold on; plot(w_gt(3,:),'k');
end

%% Helper function: compute scores
function tab = computeScores(problem, gt, teaser, est, video, adds, threshold)
models_dir = "~/research/tracking/datasets/YCBInEOAT/models/";
if (problem.object == "cracker") || (problem.object == "sugar")
    pcfiles = models_dir + ["cracker.ply", "sugar.ply", "jello.ply"];
elseif (problem.object == "mustard") || (problem.object == "bleach")
    pcfiles = models_dir + ["mustard.ply", "bleach.ply"];
elseif (problem.object == "tomato")
    pcfiles = models_dir + ["coffee.ply", "tomato.ply", "tuna.ply"];
end
gt.c = find(contains(pcfiles, problem.object));
teaser.c = gt.c*ones(length(teaser.p),1);

% Compute scores!
[add_ours, adds_ours] = get_adds(gt, est, pcfiles, adds); % takes a while if calc adds
score_add_ours = get_auc(add_ours, threshold);
score_adds_ours = get_auc(adds_ours, threshold);

[add_teaser, adds_teaser] = get_adds(gt, teaser, pcfiles, false); % takes a while
score_add_teaser = get_auc(add_teaser, threshold);
score_adds_teaser = get_auc(adds_teaser, threshold);

% put into table
methods = ["CAST"; "TEASER"];
add_adds = [score_add_ours*100, score_adds_ours*100;
            score_add_teaser*100, score_adds_teaser*100];
vid = [video; video];
tab = table(methods,add_adds,vid);
end

%% Helper function: save to JSON
function save2json(est, problem)
L_big = length(est.p);
T_est = repmat(eye(4),[1,1,L_big]);
for l = 1:L_big
    T_est(1:3,1:3,l) = est.R(:,:,l);
    T_est(1:3,4,l) = est.p(:,:,l)*1000.0;
end

fid = fopen(problem.json); 
raw = fread(fid,inf); 
str = char(raw'); 
fclose(fid);
data = jsondecode(str);

for l = 1:length(T_est)
    data(l).cast_pose = T_est(:,:,l);
end

cocoString = jsonencode(data, "PrettyPrint",true);
fid = fopen(problem.savefile, 'w');
fprintf(fid, '%s', cocoString);
fclose(fid);
end

%% Helper function: get errors
function est = getErrors(gt, est)
    est.p_err = squeeze(vecnorm(est.p - gt.p));
    est.R_err = zeros(length(est.R),1);
    for l = 1:length(est.R)
        est.R_err(l) = getAngularError(gt.R(:,:,l),est.R(:,:,l));
    end
end

%% Helper fuction: save everything together
function saveAll(savename, params)
    
end