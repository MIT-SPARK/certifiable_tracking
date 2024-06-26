function plotvariable(results, indepVar, var, settings)
%% Plot a named struct variable on an existing figure
% line plot with errorshading according to the IQR
    hold on
    if isfield(settings, "PACEEKF")
        errorshade([results.(indepVar)],[results.(var + "_ekf")],hex2rgb(settings.PACEEKF{5}));
    end
    errorshade([results.(indepVar)],[results.(var + "_pace")],hex2rgb(settings.PACERAW{5}));

    if isfield(settings, "CASTP")
        errorshade([results.(indepVar)],[results.(var + "_castp")],hex2rgb(settings.CASTP{5}));
    end
    
    Llist = settings.Llist;
    Ldomain = settings.Ldomain;
    ours_colors = settings.ours_colors;
    for lidx=Llist
        lrange = lidx + length(Ldomain)*(0:length([results.(indepVar)])-1);
        res = [results.(var + "_ours")];
        errorshade([results.(indepVar)],res(:,lrange),hex2rgb(ours_colors(lidx)));
    end

    if isfield(settings, "PACEEKF")
        plot([results.(indepVar)],median([results.(var + "_ekf")]),settings.PACEEKF{:});
    end
    plot([results.(indepVar)],median([results.(var + "_pace")]),settings.PACERAW{:});
    if isfield(settings, "CASTP")
        plot([results.(indepVar)],median([results.(var + "_castp")]),settings.CASTP{:});
    end
    for i=1:length(Llist)
        lidx = Llist(i);
        lrange = lidx + length(Ldomain)*(0:length([results.(indepVar)])-1);
        res = [results.(var + "_ours")];
        plotsettings = settings.OURS{i}; L = Ldomain(lidx);
        plotsettings{3} = plotsettings{3} + "-" + string(L);
        plotsettings{7} = ours_colors(lidx);
        plot([results.(indepVar)],median(res(:,lrange)),plotsettings{:});
    end
    hold off
end