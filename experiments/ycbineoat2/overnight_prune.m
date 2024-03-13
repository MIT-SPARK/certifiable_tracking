clc; close all; clear;

domain = 2:9;

for d = domain
    try
        A_prune(d);
    catch
    end
end