%% RAL Experiment: How do measurement noise and # shapes affect performance?
% Dataset: pascal+car
% Constants: K, N, spiral trajectory
% Independent variable: measurement noise, L
% Dependent variables: runtime, duality gap, accuracy (p, R, c)
%
% Lorenzo Shaikewitz for SPARK Lab

clc; clear; close all

%% Experiment settings
indepVar = "noiseSigmaSqrt";
savename = "pascalaeroplane_mletest_" + indepVar;
lengthScale = 0.2; % smallest dimension
domain = 0.01:0.01:0.24;  % 0:0.025:1
Ldomain = [4,8,12]; % 2: 3:12;
num_repeats = 50; % 2: 50

%% Loop
results = cell(length(domain),1);
parfor (index = 1:length(domain)) % PAR, 20)
iv = domain(index);
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
problem.category = "aeroplane";
problem.L = max(Ldomain); % nr of keyframes in horizon

problem.outlierRatio = 0.0;
problem.noiseSigmaSqrt = iv*lengthScale; % [m]

% MLE parameters
problem.accelerationNoiseBoundSqrt = 0.05*lengthScale;
problem.rotationKappa = 1/(0.05*lengthScale)*1/2;

problem.covar_measure_base = problem.noiseSigmaSqrt^2;
problem.covar_velocity_base = problem.accelerationNoiseBoundSqrt;%^2;
problem.kappa_rotrate_base = problem.rotationKappa;

% problem.covar_measure_base = 0.0001; % 2: 0.01
% problem.covar_velocity_base = 0.001;
% problem.covar_rotrate_base = 0.001;

problem.noiseBound = 3*iv*lengthScale; %3*iv for 8 reg
problem.processNoise = 5e-2;

problem.translationBound = 10.0;
problem.velocityBound = 5.0;
problem.dt = 1.0;

problem.velprior = "body";       % constant body frame velocity

% problem.accelerationNoiseBoundSqrt = 0;
% problem.rotationNoiseBound = 0; % rad

% add shape, measurements, outliers
problem = gen_pascal_tracking(problem);
lambda = 0;
problem.lambda = lambda;

% Solve!
pace = pace_raw(problem);
paceekf = pace_ekf(problem,pace); % try pace_py_ukf (issues with parfor)

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
    pid = string(feature("getpid"));
    lproblem.sdp_filename = "sdpdata" + pid + L;

    lproblem.L = L;
    lproblem.y = problem.y(:,(end-L+1):end);
    soln = solve_weighted_tracking(lproblem);

    R_err_ours(lidx) = getAngularError(problem.R_gt(:,:,end), soln.R_est(:,:,end));
    p_err_ours(lidx) = norm(problem.p_gt(:,:,end) - soln.p_est(:,:,end));
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

%%
load("../datasets/results/" + savename + ".mat","results")

%% Display Results
% data settings
Llist = [1,2,3];
displayRange = 1:length(results);

% visual settings
tile = true;

settings.PACEEKF = {'x-','DisplayName', 'PACE-EKF', 'Color', "#D95319",'LineWidth',1};
settings.PACERAW = {'x-','DisplayName', 'PACE-RAW', 'Color', "#EDB120",'LineWidth',1};

settings.OURS = {'x-','DisplayName', 'OURS','LineWidth',2,'Color','002e4c'};
settings.ours_colors = ["#338eca","#005b97","#002e4c"];
settings.Llist = Llist;
settings.Ldomain = Ldomain;

% restrict domain and normalize positions
resultsAdj = results(:,displayRange);
for j = 1:length(resultsAdj)
resultsAdj(j).p_err_ekf = resultsAdj(j).p_err_ekf/lengthScale;
resultsAdj(j).p_err_pace = resultsAdj(j).p_err_pace/lengthScale;
resultsAdj(j).p_err_ours = resultsAdj(j).p_err_ours/lengthScale;
end

% created tiled figure
if tile
    figure
    t=tiledlayout(1,5);
    title(t,'Measurement Noise')
end

% Positions
if (tile); nexttile; else; figure; end
plotvariable(resultsAdj, indepVar, "p_err", settings)
yscale log;% xscale log
xlabel(indepVar); ylabel("Position Error (normalized)");
title("Position Errors")

settings=rmfield(settings,"PACEEKF");
% Rotations
if (tile); nexttile; else; figure; end
plotvariable(resultsAdj, indepVar, "R_err", settings)
yscale log;% xscale log
xlabel(indepVar); ylabel("Rotation Error (deg)");
title("Rotation Errors")

% Shape
if (tile); nexttile; else; figure; end
plotvariable(resultsAdj, indepVar, "c_err", settings)
yscale log;% xscale log
xlabel(indepVar); ylabel("Shape Error (normalized)");
title("Shape Errors")

% Gap
if (tile); nexttile; else; figure; end
plotvariable(resultsAdj, indepVar, "gap", settings)
yscale log;% xscale log
xlabel(indepVar); ylabel("Gap");
title("Suboptimality Gaps")

% Time
if (tile); nexttile; else; figure; end
plotvariable(resultsAdj, indepVar, "time", settings)
yscale log;% xscale log
xlabel(indepVar); ylabel("Time (s)");
title("Run Times")

%% Display Results
% process into displayable form
% settings.OURS = {'DisplayName', 'OURS','LineWidth',3};
% ours_colors = ["#000000", "#002e4c", "#00395f", "#004471", "#005084",...
%                "#005b97", "#0067aa","#0072bd", "#1a80c4", "#338eca",...
%                "#4d9cd1","#66aad7","#80b9de"];
ours_colors = ["#002e4c", "#005b97","#338eca","#80b9de"];
settings.PACEEKF = {'DisplayName', 'PACE-EKF', 'Color', "#D95319"};
settings.PACERAW = {'DisplayName', 'PACE-RAW', 'Color', "#EDB120"};
figure
t=tiledlayout(2,2);
title(t,'Measurement Noise')
set(0,'DefaultLineLineWidth',1)

% display_range = 2:length(domain);
% display_range = 1:17;
% results = results(:,display_range);
Llist = [1,2,3];

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
lrange = lidx + length(Ldomain)*(0:length(display_range)-1);
plotsettings = {'DisplayName', "OURS-" + string(L),'LineWidth',3,'Color',ours_colors(length(Ldomain)-lidx+1)};
a=plot([results.(indepVar)],median(res(:,lrange)),'x-',plotsettings{:});
errorshade([results.(indepVar)],res(:,lrange),get(a,'Color'));
end
yscale log;% xscale log
xlabel(indepVar); ylabel("Rotation Error (deg)");
title("Rotation Errors")

% position figure
nexttile
hold on
b=plot([results.(indepVar)],median([results.p_err_ekf])/lengthScale,'x-',settings.PACEEKF{:});
errorshade([results.(indepVar)],[results.p_err_ekf]/lengthScale,get(b,'Color'));
c=plot([results.(indepVar)],median([results.p_err_pace])/lengthScale,'x-',settings.PACERAW{:});
errorshade([results.(indepVar)],[results.p_err_pace]/lengthScale,get(c,'Color'));
res = [results.p_err_ours];
for lidx = Llist
L = Ldomain(lidx);
lrange = lidx + length(Ldomain)*(0:length(display_range)-1);
plotsettings = {'DisplayName', "OURS-" + string(L),'LineWidth',3,'Color',ours_colors(length(Ldomain)-lidx+1)};
a=plot([results.(indepVar)],median(res(:,lrange))/lengthScale,'x-',plotsettings{:});
errorshade([results.(indepVar)],res(:,lrange)/lengthScale,get(a,'Color'));
end
yscale log;% xscale log
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
lrange = lidx + length(Ldomain)*(0:length(display_range)-1);
plotsettings = {'DisplayName', "OURS-" + string(L),'LineWidth',3,'Color',ours_colors(length(Ldomain)-lidx+1)};
a=plot([results.(indepVar)],median(res(:,lrange)),'x-',plotsettings{:});
errorshade([results.(indepVar)],res(:,lrange),get(a,'Color'));
end
yscale log;% xscale log
xlabel(indepVar); ylabel("Shape Error (normalized)");
title("Shape Errors")

% gap figure
nexttile
b=plot([results.(indepVar)],median([results.gap_pace]),'x-',settings.PACERAW{:});
hold on
errorshade([results.(indepVar)],[results.gap_pace],get(b,'Color'));
res = [results.gap_ours];
% res(res < 0) = 0;
% res(res > 1) = 1; % only affects first round
for lidx = Llist
L = Ldomain(lidx);
lrange = lidx + length(Ldomain)*(0:length(display_range)-1);
plotsettings = {'DisplayName', "OURS-" + string(L),'LineWidth',3,'Color',ours_colors(length(Ldomain)-lidx+1)};
a=semilogy([results.(indepVar)],median(res(:,lrange)),'x-',plotsettings{:});
errorshade([results.(indepVar)],res(:,lrange),get(a,'Color'));
end
yscale log;% xscale log
xlabel(indepVar); ylabel("Gap");
title("Suboptimality Gaps")

% time figure
nexttile
hold on
b=plot([results.(indepVar)],median([results.time_pace]),'x-',settings.PACERAW{:});
hold on
errorshade([results.(indepVar)],[results.time_pace],get(b,'Color'));
res = [results.time_ours];
for lidx = Llist
L = Ldomain(lidx);
lrange = lidx + length(Ldomain)*(0:length(display_range)-1);
plotsettings = {'DisplayName', "OURS-" + string(L),'LineWidth',3,'Color',ours_colors(length(Ldomain)-lidx+1)};
a=plot([results.(indepVar)],median(res(:,lrange)),'x-',plotsettings{:});
errorshade([results.(indepVar)],res(:,lrange),get(a,'Color'));
end
yscale log;% xscale log
xlabel(indepVar); ylabel("Time (s)");
title("Solve Time")

%% Box and Whisker Display

% parameters
nBins = length(domain);

% make bins
domain = [results.(indepVar)];
edges = linspace(domain(1),domain(end),nBins+1);
bins = discretize(domain,edges);

% generate figure
figure
% t=tiledlayout(2,2);
% title(t,'Measurement Noise')

% rotation figure
% nexttile


hold on
L = Ldomain(lidx);
lrange = lidx + length(Ldomain)*(0:length(display_range)-1);
res = [results.p_err_ours];
boxplot(res(:,lrange), bins, 'Positions', 1:nBins ,'Symbol','','Colors',hex2rgb('#1a80c4'),'PlotStyle','compact','Widths',0.25)
boxplot([results.p_err_pace], bins, 'Positions', (1:nBins) + 0.3,'Symbol','','Colors',hex2rgb('#D95319'),'PlotStyle','compact','Widths',0.25)
hold off
yscale log

% set(gca,'xticklabel',{'#1','#2','#3','#4','#5','#6','#7','#8'})

%% Second Box and Whisker Option

tbl = table;
domain = [results.(indepVar)];
tbl.method = repelem(["OURS", "PACE"], num_repeats*length(domain))';

nBins = length(domain);
edges = linspace(domain(1),domain(end),nBins+1);
bins = discretize(domain,edges);
tbl.domain = repmat(repelem(bins,num_repeats)',[2,1]);

% positions
res = [results.p_err_ours];
lidx = 1;
lrange = lidx + length(Ldomain)*(0:length(display_range)-1);
p_err_ours = res(:,lrange)/lengthScale;
p_err_pace = [results.p_err_pace]/lengthScale;

tbl.p = [p_err_ours(:); p_err_pace(:)];

figure
b=boxchart(tbl.domain,tbl.p,'GroupByColor',tbl.method, 'JitterOutliers','on','MarkerStyle','.', 'Notch','off');
yscale log