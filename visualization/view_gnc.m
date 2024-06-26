function view_gnc(problem, info)
%% View the results of GNC
% Lorenzo Shaikewitz for SPARK Lab

set(0,'DefaultLineLineWidth',2)

% Create figure
global itr;
global first;
first = true;
itr = 1;
figure('KeyPressFcn',{@keypress, "m", problem, info})
plotItr(problem, info)

% plot the current shape of the object
% figure('KeyPressFcn',{@keypress, "c", problem, info})
% plotc(problem, info)


end

function plotItr(problem, info)
    global itr;
    global first;
    ax = gca;
    v = get(ax, 'View');

    % preliminaries
    L = problem.L;
    N = problem.N_VAR;
    y = reshape(problem.y,[3,N,L]);

    if isfield(problem,'prioroutliers')
        weights = zeros(N*L,size(info.history.weights,2));
        for iteration = 1:size(info.history.weights,2)
            w = info.history.weights(:,iteration)';
            % add priors to weights
            for i = 1:length(problem.prioroutliers)
                o = problem.prioroutliers(i);
                w = [w(1:o-1),0.0,w(o:end)];
            end
            weights(:,iteration) = w';
        end
    else
        weights = info.history.weights;
    end

    for l = 1:L
        % plot weighted keypoint measurements
        tbl = array2table(y(:,:,l)','VariableNames',{'x','y','z'});
        tbl.weight = weights((N*(l-1)+1):(N*l),itr);
        s=scatter3(tbl,'x','y','z','filled','ColorVariable','weight');
        s.SizeData=100;

        hold on
        % plot gt pose
        p = problem.p_gt;
        R = problem.R_gt;
        quiver3(p(1,l)',p(2,l)',p(3,l)', ...
            squeeze(R(1,1,l)),squeeze(R(2,1,l)),squeeze(R(3,1,l)),0.25,'r');
        quiver3(p(1,l)',p(2,l)',p(3,l)', ...
            squeeze(R(1,2,l)),squeeze(R(2,2,l)),squeeze(R(3,2,l)),0.25,'g');
        quiver3(p(1,l)',p(2,l)',p(3,l)', ...
            squeeze(R(1,3,l)),squeeze(R(2,3,l)),squeeze(R(3,3,l)),0.25,'b');

        % plot est pose
        p = info.history.solns{itr}.p_est;
        R = info.history.solns{itr}.R_est;
        quiver3(p(1,l)',p(2,l)',p(3,l)', ...
            squeeze(R(1,1,l)),squeeze(R(2,1,l)),squeeze(R(3,1,l)),0.25,'r');
        quiver3(p(1,l)',p(2,l)',p(3,l)', ...
            squeeze(R(1,2,l)),squeeze(R(2,2,l)),squeeze(R(3,2,l)),0.25,'g');
        quiver3(p(1,l)',p(2,l)',p(3,l)', ...
            squeeze(R(1,3,l)),squeeze(R(2,3,l)),squeeze(R(3,3,l)),0.25,'b');
        plot3(p(1,l),p(2,l),p(3,l),'x','MarkerSize',15,'LineWidth',3)
    end
    hold off
    axis equal tight manual

    load("guppy.mat",'guppy')
    set(ax, 'Colormap',guppy,'CLim',[0,1]);
    cb = colorbar;
    % cb.Layout.Tile = 'south';
    title(sprintf("itr=%d",itr))

    if first
        plot_view = get(ax, 'View');
        first = false;
    else
        view(v)
    end
end

function plotc(problem, info)
    global itr

    % plot estimated shape
    shape = problem.B*info.history.solns{itr}.c_est;
    shape = reshape(shape, [3,problem.N_VAR]);

    s=scatter3(shape(1,:),shape(2,:),shape(3,:),'filled');
    s.SizeData=100;

    title(sprintf("itr=%d",itr))

    % plot gt shape?
    hold on
    shape = problem.B*problem.c_gt;
    shape = reshape(shape, [3,problem.N_VAR]);

    s=scatter3(shape(1,:),shape(2,:),shape(3,:),'x','LineWidth',4);
    s.SizeData=100;
    hold off
    
    axis equal
    axis tight manual
end

function keypress(src, event, key, problem, info)
    max_itr = info.Iterations;
    global itr

    switch string(event.Key)
        case "downarrow"
            itr = itr + 1;
            overwrite_plot = true;
        case "rightarrow"
            itr = itr + 1;
            overwrite_plot = true;
        case "uparrow"
            itr = itr - 1;
            overwrite_plot = true;
        case "leftarrow"
            itr = itr - 1;
            overwrite_plot = true;
        otherwise
            overwrite_plot = false;
    end
    if (itr > max_itr)
        itr = max_itr;
        overwrite_plot = false;
    elseif (itr < 1)
        itr = 1;
        overwrite_plot = false;
    end
    
    if overwrite_plot
        switch key
            case "m"
                plotItr(problem, info)
            case "c"
                plotc(problem,info)
        end
    end
end