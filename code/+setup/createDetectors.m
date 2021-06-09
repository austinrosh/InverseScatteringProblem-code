%% createDetectors.m
% creates detector locations to evaluate the scattered field at specified
% points xd, yd, and zd (an array of such points).
function [xd,yd,zd] = createDetectors(Ndetect_x,Ndetect_y)
    xd = linspace(-1,2,Ndetect_x);
    yd = linspace(-1,2,Ndetect_y);
    zd = 1.5; %hold z position fixed
end

