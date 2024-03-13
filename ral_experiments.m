setup_func("pascal")
try
    EXP_outlier_ablation;
catch
end
try
    EXP_Landm; % try changing to ukf
catch
end
try
    EXP_Landp; % try changing to ukf
catch
end