function score = get_auc(metric, threshold)
%% step 2: compute area under curve (AUC)
% curve in question is accuracy-threshold curve
% see pose-cnn figure 8
% generate curve
thresh = linspace(0,threshold,100); % x-axis
accuracy = zeros(length(thresh),1); % y-axis
for t = 1:length(thresh)
    accuracy(t) = sum(metric < thresh(t))/length(metric);
end
% figure
% plot(thresh,accuracy);

% area under curve!
max_score = threshold*1;
score = trapz(thresh, accuracy) / max_score;

end