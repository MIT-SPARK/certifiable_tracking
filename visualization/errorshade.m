%% Plotting Helper
function errorshade(x,y,color,type)

if nargin < 4
    type = "shade";
end

curve1 = prctile(y,75);
curve2 = prctile(y,25);
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];

if (type=="shade")
    % shade
    obj = fill(x2, inBetween,color,'FaceAlpha',0.15,'EdgeColor','none');
    % obj1 = plot(x,curve1,'--','Color',color);
    % obj2 = plot(x,curve2,'--','Color',color);
    % obj1.Annotation.LegendInformation.IconDisplayStyle = "off";
    % obj2.Annotation.LegendInformation.IconDisplayStyle = "off";
elseif (type=="bar")
    % error bar
    med = median(y);
    obj = errorbar(x, med, med - curve2, curve1 - med,"LineStyle","none",'Color',color,'LineWidth',2);
elseif (type=="dotted")
    % dotted boundary lines
    obj = plot(x,curve1,'--','Color',color);
end

obj.Annotation.LegendInformation.IconDisplayStyle = "off";

% TF = isoutlier(y);
% x_rep = repmat(x,[size(y,1),1]);
% plot(x_rep(TF),y(TF),"x",'MarkerFaceColor',colors{i},'MarkerEdgeColor',colors{i},'MarkerSize',10);
end