% setup_func("ycbineoat2")
% domain = [4,5,6];
% failed = [];
% 
% for d = domain
%     try
%         ycbineoat(d, false);
%     catch
%         failed = [failed, d];
%         fprintf("--------")
%     end
% end

setup_func("pascal")
failed = [];
% try
%     EXP_outlier_ablation;
% catch
% end
try
    EXP_Landm_UKF;
catch
    % failed = [failed,-1];
end
try
    EXP_Landp_UKF;
catch
    % failed = [failed,-2];
end