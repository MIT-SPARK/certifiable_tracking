function [degcm, p_err, R_err] = compute_degcm(gt, est, varargin)
%% Compute xdeg xcm metric
% This is the percent of poses within x cm and x deg of target
% gt and est are structs with p, R fields
% 
% Lorenzo Shaikewitz for SPARK Lab

%% default to 5deg5cm
params = inputParser;
params.addParameter('cmThreshold', 5, @(x) isscalar(x));
params.addParameter('degThreshold', 5, @(x) isscalar(x));
params.parse(varargin{:});

mThreshold = params.Results.cmThreshold / 100;
degThreshold = params.Results.degThreshold;

%% Compute errors
degcm_count = 0;
p_err = squeeze(vecnorm(gt.p - est.p));
R_err = NaN*ones(length(p_err),1);
for l = 1:length(p_err)
    R_err(l) = getAngularError(gt.R(:,:,l), est.R(:,:,l));

    if (R_err(l) <= degThreshold) && (p_err(l) <= mThreshold)
        degcm_count = degcm_count + 1;
    end
end
degcm = (degcm_count / length(p_err)) * 100; % percentage

end