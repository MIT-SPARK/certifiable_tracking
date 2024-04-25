%% Step 3: Solve Reduced Version
% Load solved spherical problems and rerun without spherical
% Run Experiment Settings section before this
%
% Lorenzo Shaikewitz for SPARK Lab

load(savename + "_b")
% load(savename + "_c") % problems, solns
load(savename + "_gt");

%% Check solutions
% solns = [solns{:}];
L = maxL;
N = problems{1}.N_VAR;

p_err = zeros(L*N,L)*NaN;
R_err = zeros(L*N,L)*NaN;

est_legit.p = zeros(3,1,length(solns))*NaN;
est_legit.R = zeros(3,3,length(solns))*NaN;
est_ignoringbad.p = zeros(3,1,length(solns))*NaN;
est_ignoringbad.R = zeros(3,3,length(solns))*NaN;
est_bestrun.p = zeros(3,1,length(solns))*NaN;
est_bestrun.R = zeros(3,3,length(solns))*NaN;
gaps = [];

figure(1);
for j = 1:length(solns)

problem = problems{j};
soln = solns{j};
L_cur = problem.L;

% eigenvalue plot
% figure; bar(eig(soln.raw.Xopt{1})); % if rank = 1, then relaxation is exact/tight
% hold on

% Plot trajectory!
figure(1);
axis equal
p_est = reshape(soln.p,[3,L_cur,1]);
plot3(p_est(1,:),p_est(2,:),p_est(3,:),'.k', 'MarkerSize',10);
hold on

R_est = soln.R;
quiver3(p_est(1,:)',p_est(2,:)',p_est(3,:)',squeeze(R_est(1,1,:)),squeeze(R_est(2,1,:)),squeeze(R_est(3,1,:)),'r');
quiver3(p_est(1,:)',p_est(2,:)',p_est(3,:)',squeeze(R_est(1,2,:)),squeeze(R_est(2,2,:)),squeeze(R_est(3,2,:)),'g');
quiver3(p_est(1,:)',p_est(2,:)',p_est(3,:)',squeeze(R_est(1,3,:)),squeeze(R_est(2,3,:)),squeeze(R_est(3,3,:)),'b');

idx = problem.startIdx:(problem.startIdx + L_cur-1);
for l = 1:L_cur
    p_err(idx(l),l) = norm(soln.p(:,:,l) - gt.p(:,:,idx(l)));
    R_err(idx(l),l) = getAngularError(gt.R(:,:,idx(l)),soln.R(:,:,l));
end


% compute errors
% 1) legit: use latest estimate no matter what
if j == 1
    est_legit.p(:,:,idx) = soln.p;
    est_legit.R(:,:,idx) = soln.R;
else
    est_legit.p(:,:,idx(end)) = soln.p(:,:,end);
    est_legit.R(:,:,idx(end)) = soln.R(:,:,end);
end
% 2) ignoringbad: throw away bad
if (soln.gap < 0.3)
    est_ignoringbad.p(:,:,idx(end)) = soln.p(:,:,end);
    est_ignoringbad.R(:,:,idx(end)) = soln.R(:,:,end);
end
% 3) bestrun: use estimate that has best run (by gap)
gaps = [gaps; abs(soln.gap)];
for l = 1:L_cur-1
    % for each pose estimate, excluding the current one
    testrange = length(gaps) - (1:(L_cur - l));
    newgood = true;
    for l_test = testrange
        if (l_test < 1)
            break
        end
        % for each other gap we can compare with
        if (abs(soln.gap) > gaps(l_test))
            % fails if any of them have better gaps
            newgood = false;
            break
        end
    end
    if (newgood)
        % update!
        est_bestrun.p(:,:,idx(l)) = soln.p(:,:,l);
        est_bestrun.R(:,:,idx(l)) = soln.R(:,:,l);
    end
end
est_bestrun.p(:,:,idx(end)) = soln.p(:,:,end);
est_bestrun.R(:,:,idx(end)) = soln.R(:,:,end);


end

est_ignoringbad.p = est_ignoringbad.p(:,:,1:end-1);
est_legit.p = est_legit.p(:,:,1:end-1);
est_bestrun.p = est_bestrun.p(:,:,1:end-1);
est_ignoringbad.R = est_ignoringbad.R(:,:,1:end-1);
est_legit.R = est_legit.R(:,:,1:end-1);
est_bestrun.R = est_bestrun.R(:,:,1:end-1);

%% Plot Ground Truth

figure
p_gt = reshape(gt.p,[3,size(gt.p,3),1]);
plot3(p_gt(1,:),p_gt(2,:),p_gt(3,:),'.k', 'MarkerSize',10);
hold on
axis equal

R_est = soln.R;
quiver3(p_gt(1,:)',p_gt(2,:)',p_gt(3,:)',squeeze(gt.R(1,1,:)),squeeze(gt.R(2,1,:)),squeeze(gt.R(3,1,:)),'r');
quiver3(p_gt(1,:)',p_gt(2,:)',p_gt(3,:)',squeeze(gt.R(1,2,:)),squeeze(gt.R(2,2,:)),squeeze(gt.R(3,2,:)),'g');
quiver3(p_gt(1,:)',p_gt(2,:)',p_gt(3,:)',squeeze(gt.R(1,3,:)),squeeze(gt.R(2,3,:)),squeeze(gt.R(3,3,:)),'b');

%% ADD and ADD-S scores
models_dir = "~/research/tracking/datasets/YCBInEOAT/models/";
pcfile = models_dir + "cracker" + ".ply";
pcfile_gt = pcfile;
pcfile_est = pcfile;

% Compute scores!
threshold = 0.1;
[add_ours, adds_ours] = get_adds(gt, est_legit, pcfile_gt); % takes a while
score_add_ours = get_auc(add_ours, threshold);
score_adds_ours = get_auc(adds_ours, threshold);

[add_teaser, adds_teaser] = get_adds(gt, teaser, pcfile_gt); % takes a while
score_add_teaser = get_auc(add_teaser, threshold);
score_adds_teaser = get_auc(adds_teaser, threshold);

%%
methods = ["CAST"; "TEASER"];
add_adds = [score_add_ours*100, score_adds_ours*100;
            score_add_teaser*100, score_adds_teaser*100];
vid = [video; video];
tab_new = table(methods,add_adds,vid);

%% Save Poses into JSON
est = est_legit;
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