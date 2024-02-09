%% Plotting Helper
function errorshade(x,y,color)

curve1 = prctile(y,75);
curve2 = prctile(y,25);
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];
obj = fill(x2, inBetween,color,'FaceAlpha',0.2,'EdgeColor','none');
obj.Annotation.LegendInformation.IconDisplayStyle = "off";

% TF = isoutlier(y);
% x_rep = repmat(x,[size(y,1),1]);
% plot(x_rep(TF),y(TF),"x",'MarkerFaceColor',colors{i},'MarkerEdgeColor',colors{i},'MarkerSize',10);
end