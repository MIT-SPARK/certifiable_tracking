clc; close all; clear;

domain = [2,3,4,5,6];
failed = [];

for d = domain
    try
        ycbineoat(d, true);
    catch
        failed = [failed, d];
        fprintf("--------")
        playSound();
    end
end

function playSound
load("handel.mat")
obj = audioplayer(y,Fs);
obj.TimerFcn = 'showSeconds';
obj.TimerPeriod = 5;

playblocking(obj);
end