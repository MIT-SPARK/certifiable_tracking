%% RAL Experiment: How do process noise and # shapes affect performance?
% Dataset: pascal+car
% Constants: K, N, measurement noise
% Independent variable: process noise, L
% Dependent variables: runtime, duality gap, accuracy (p, R, c)
%
% Lorenzo Shaikewitz for SPARK Lab

clc; clear; close all

%% Experiment settings
indepVar = "accelerationNoiseBoundSqrt";
savename = "pascalcar_fixed_" + indepVar;
domain = 0:0.05:2;
Ldomain = 3:12;
num_repeats = 100;

%% Loop
results = cell(length(domain),1);
for index = 1:length(domain)
iv = domain(index)
resultsIV = struct();
resultsIV.(indepVar) = iv;
resultsIV.R_err_ours = zeros(num_repeats,length(Ldomain));
resultsIV.R_err_ekf = zeros(num_repeats,1);
resultsIV.R_err_pace = zeros(num_repeats,1);
resultsIV.p_err_ours = zeros(num_repeats,length(Ldomain));
resultsIV.p_err_ekf = zeros(num_repeats,1);
resultsIV.p_err_pace = zeros(num_repeats,1);
resultsIV.c_err_ours = zeros(num_repeats,length(Ldomain));
resultsIV.c_err_pace = zeros(num_repeats,1);
resultsIV.gap_ours = zeros(num_repeats,length(Ldomain));
resultsIV.gap_pace = zeros(num_repeats,1);
resultsIV.time_ours = zeros(num_repeats,length(Ldomain));
resultsIV.time_pace = zeros(num_repeats,1);
disp("Starting " + indepVar + "=" + string(iv));
for j = 1:num_repeats

problem = struct();
problem.category = "car";
problem.L = max(Ldomain); % nr of keyframes in horizon

problem.outlierRatio = 0.0;
problem.noiseSigmaSqrt = 0.5*0.3; % [m] (0.3 is length scale for car)
% 0.1*0.3 in original run
problem.velocity_weight_multiplier = 1;
problem.rotrate_kappa_multiplier = 1;

problem.noiseBound = 0.5;
problem.processNoise = 0.5;

problem.translationBound = 10.0;
problem.velocityBound = 2.0;
problem.dt = 1.0;

problem.velprior = "body";       % constant body frame velocity

problem.accelerationNoiseBoundSqrt = iv*0.3;
problem.rotationNoiseBound = 0; % rad

% add shape, measurements, outliers
problem = gen_pascal_tracking(problem);
lambda = 0;
problem.lambda = lambda;

% Solve!
pace = pace_raw(problem);
paceekf = pace_py_UKF(problem,pace);

% Save solutions: only use last error
% rotation error
R_err_ekf = getAngularError(problem.R_gt(:,:,end), paceekf.R(:,:,end));
R_err_pace = getAngularError(problem.R_gt(:,:,end), pace.R(:,:,end));
% position error
p_err_ekf = norm(problem.p_gt(:,:,end) - paceekf.p(:,:,end));
p_err_pace = norm(problem.p_gt(:,:,end) - pace.p(:,:,end));
% shape error
c_err_pace = norm(problem.c_gt - pace.c(:,:,end));
% time and gap
gap_pace = pace.gaps(end);
time_pace = pace.times(end);

R_err_ours = zeros(1,length(Ldomain));
p_err_ours = zeros(1,length(Ldomain));
c_err_ours = zeros(1,length(Ldomain));
gap_ours = zeros(1,length(Ldomain));
time_ours = zeros(1,length(Ldomain));
for lidx = 1:length(Ldomain)
    L = Ldomain(lidx);
    lproblem = problem;
    % regen only first time
    lproblem.regen_sdp = (j == 1);
    pid = string(feature("getpid"))
    lproblem.sdp_filename = "sdpdata" + pid + L;

    lproblem.L = L;
    lproblem.y = problem.y(:,(end-L+1):end);
    soln = solve_weighted_tracking(lproblem);

    R_err_ours(lidx) = getAngularError(problem.R_gt(:,:,end), soln.R_est(:,:,end));
    p_err_ours(lidx) = norm(problem.p_gt(:,:,end) - soln.p_est(:,:,end))
    c_err_ours(lidx) = norm(problem.c_gt - soln.c_est);
    gap_ours(lidx) = soln.gap;
    time_ours(lidx) = soln.solvetime;
end

% save
resultsIV.R_err_ours(j,:) = R_err_ours;
resultsIV.R_err_ekf(j)  = R_err_ekf;
resultsIV.R_err_pace(j) = R_err_pace;
resultsIV.p_err_ours(j,:) = p_err_ours;
resultsIV.p_err_ekf(j)  = p_err_ekf;
resultsIV.p_err_pace(j) = p_err_pace;
resultsIV.c_err_ours(j,:) = c_err_ours;
resultsIV.c_err_pace(j) = c_err_pace;
resultsIV.gap_ours(j,:) = gap_ours;
resultsIV.gap_pace(j) = gap_pace;
resultsIV.time_ours(j,:) = time_ours;
resultsIV.time_pace(j) = time_pace;
end
results{index} = resultsIV;
end
results = [results{:}];
% save
save("../datasets/results/" + savename + ".mat","results")

%% Display Results
% process into displayable form
% settings.OURS = {'DisplayName', 'OURS','LineWidth',3};
ours_colors = ["#000000", "#002e4c", "#00395f", "#004471", "#005084",...
               "#005b97", "#0067aa","#0072bd", "#1a80c4", "#338eca",...
               "#4d9cd1","#66aad7","#80b9de"];
settings.PACEEKF = {'DisplayName', 'PACE-EKF', 'Color', "#D95319"};
settings.PACERAW = {'DisplayName', 'PACE-RAW', 'Color', "#EDB120"};
figure
tiledlayout(2,2);
set(0,'DefaultLineLineWidth',2)

display_range = 1:length(domain);
Llist = [1,4,7,10];

% Rotation figure
nexttile
hold on
c=plot([results.(indepVar)],median([results.R_err_pace]),'x-',settings.PACERAW{:});
errorshade([results.(indepVar)],[results.R_err_pace],get(c,'Color'));
% b=plot([results.(indepVar)],median([results.R_err_ekf]),'x-',settings.PACEEKF{:});
% errorshade([results.(indepVar)],[results.R_err_ekf],get(b,'Color'));
res = [results.R_err_ours];
for lidx = Llist
L = Ldomain(lidx);
lrange = lidx + length(Ldomain)*(0:length(domain)-1);
plotsettings = {'DisplayName', "OURS-" + string(L),'LineWidth',3,'Color',ours_colors(length(Ldomain)-lidx+1)};
a=plot([results.(indepVar)],median(res(:,lrange)),'x-',plotsettings{:});
errorshade([results.(indepVar)],res(:,lrange),get(a,'Color'));
end

xlabel(indepVar); ylabel("Rotation Error (deg)");
title("Rotation Errors")

% position figure
nexttile
hold on
b=plot([results.(indepVar)],median([results.p_err_ekf])/0.3,'x-',settings.PACEEKF{:});
errorshade([results.(indepVar)],[results.p_err_ekf]/0.3,get(b,'Color'));
c=plot([results.(indepVar)],median([results.p_err_pace])/0.3,'x-',settings.PACERAW{:});
errorshade([results.(indepVar)],[results.p_err_pace]/0.3,get(c,'Color'));
res = [results.p_err_ours];
for lidx = Llist
L = Ldomain(lidx);
lrange = lidx + length(Ldomain)*(0:length(domain)-1);
plotsettings = {'DisplayName', "OURS-" + string(L),'LineWidth',3,'Color',ours_colors(length(Ldomain)-lidx+1)};
a=plot([results.(indepVar)],median(res(:,lrange)),'x-',plotsettings{:});
errorshade([results.(indepVar)],res(:,lrange),get(a,'Color'));
end
xlabel(indepVar); ylabel("Position Error (normalized)");
title("Position Errors")

lg = legend('Orientation','horizontal');
lg.Layout.Tile = 'south';

% shape figure
nexttile
hold on
b=plot([results.(indepVar)],median([results.c_err_pace]),'x-',settings.PACERAW{:});
errorshade([results.(indepVar)],[results.c_err_pace],get(b,'Color'));
res = [results.c_err_ours];
for lidx = Llist
L = Ldomain(lidx);
lrange = lidx + length(Ldomain)*(0:length(domain)-1);
plotsettings = {'DisplayName', "OURS-" + string(L),'LineWidth',3,'Color',ours_colors(length(Ldomain)-lidx+1)};
a=plot([results.(indepVar)],median(res(:,lrange)),'x-',plotsettings{:});
errorshade([results.(indepVar)],res(:,lrange),get(a,'Color'));
end
xlabel(indepVar); ylabel("Shape Error (normalized)");
title("Shape Errors")

% gap figure
nexttile
hold on
res = abs([results.gap_ours]);
for lidx = Llist
L = Ldomain(lidx);
lrange = lidx + length(Ldomain)*(0:length(domain)-1);
plotsettings = {'DisplayName', "OURS-" + string(L),'LineWidth',3,'Color',ours_colors(length(Ldomain)-lidx+1)};
a=semilogy([results.(indepVar)],median(res(:,lrange)),'x-',plotsettings{:});
errorshade([results.(indepVar)],res(:,lrange),get(a,'Color'));
end
xlabel(indepVar); ylabel("Gap");
title("Suboptimality Gaps")

% time figure
% nexttile
% hold on
% res = [results.time_ours];
% for lidx = Llist
% L = Ldomain(lidx);
% lrange = lidx + length(Ldomain)*(0:length(domain)-1);
% plotsettings = {'DisplayName', "OURS-" + string(L),'LineWidth',3,'Color',ours_colors(length(Ldomain)-lidx+1)};
% a=plot([results.(indepVar)],median(res(:,lrange)),'x-',plotsettings{:});
% errorshade([results.(indepVar)],res(:,lrange),get(a,'Color'));
% end
% xlabel(indepVar); ylabel("Time (s)");
% title("Solve Time")

%% Change results format
% for i = 1:length(results)
%     results(i).R_err_ekf = results(i).R_err_ekf';
%     results(i).R_err_pace = results(i).R_err_pace';
%     results(i).p_err_ekf = results(i).p_err_ekf';
%     results(i).p_err_pace = results(i).p_err_pace';
%     results(i).c_err_pace = results(i).c_err_pace';
% end