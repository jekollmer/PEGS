function [ xv ] = syncirdif( dat, map, np1, np2, rkerd1, rkerd2, varargin )
% This function exists to be used in MATLAB's optimization. It takes the
% center locations found from max2D and the real kernel sizes and lines
% them up to create a synthetic image (see snythcircreate.m). This image is
% then subtracted from the original to create "holes." Funcfit should
% intake this function and use it in optimization to find subpixel accurate
% locations by minimizing the sum of hol as the minimal value will be where
% all circles are perfectly alligned and canceling one another.
%
%   Perhaps the dat input should be the original image, however, this gives
%   us that cool problem of noise and the odd reflections of the
%   photoelastic discs.
%
%	Input:  rkerd1 - Real diameter of the first kernel
%
%           rkerd2 - Real diameter of the second kernel
%
%           dat - Data; Image output from imdef preferrably.
%
%           map - Concatenated list of circle locations arranged so the
%           first np1 values are the row location of kernel 1 and the
%           np1+1:2*np1 values are column locations of kernel 1. Same
%           applies for kernel 2 but it's just the values afterwards. This
%           can be created by doing mt1 = reshape(ov1,[1,2*length(ov1)]);
%           mt2 = reshape(ov2,[1,2*length(ov2)]); map = [mt1,mt2];
%
%           np1 - Number of particles of size 1
%
%           np2 - Number of particles of size 2
%
%           varargin - optimoptions. This is for users who want more
%           controls over the optimization options used.
%
%   Output: xv - The vector of information. Same format as map.

if isempty(varargin) == 1 %Assigns default options if none are input
    options = optimoptions('lsqnonlin','Algorithm','trust-region-reflective','Display','off', ...
        'FunctionTolerance',1e-25,'StepTolerance',1e-25,'OptimalityTolerance',1e-25,'MaxFunctionEvaluations',100); %Honestly, these settings are generally okay. If you want "more accuracy" but want to waste more time, just change MaxFunctionEvaluations up to like 300+ and you'll get it
else %Checks the data input
    S = varargin{1}; 
    S = whos('S'); %Creates a definition for what the option input is
    
    if strcmp(S.class,'optim.options.Lsqnonlin') ~= 1 %If the option input is not the proper type, defaults are assigned.
        fprintf('Incorrect custom option input: Using default options.\n')
        options = optimoptions('lsqnonlin','Algorithm','trust-region-reflective','Display','off', ...
            'FunctionTolerance',1e-25,'StepTolerance',1e-25,'OptimalityTolerance',1e-25,'MaxFunctionEvaluations',100);
    end
    
end

dat = dat/(max(dat(:))); %Reduces our data down to values of -1 to 1.

xv1 = []; %Creates an empty vector to store new values
xv2 = []; %Creates a second empty vector
ap = 50;
dat = padarray(dat,[round(rkerd2)+ap,round(rkerd2)+ap],0);

dim1 = 65/100; %%%%% Dimensional constant that decides how much of our image is shown. Must be above .5 to capture the whole circle. .65 to .8 are the best range %%%%%
for l = 1:np1 %Runs a loop for optimization because of how it works
    datt1 = dat(round(map(l)-dim1*rkerd1 + rkerd2 + ap):round(map(l)+rkerd1*dim1  + rkerd2 + ap), ...
        (round(map(np1+l)-dim1*rkerd1 + rkerd2 + ap):round(map(np1+l)+rkerd1*dim1 + rkerd2 + ap))); %Sections our data best on the center location. Goes dim1*radius to the left, right, above, and below the center.
    [rt,ct] = size(datt1); %Finds the size of our data
    x0 = [ceil(rt/2), ceil(ct/2), map(2*np1 + l), map(3*np1 + l)]; %Guesses our center is in the middle of our image and pulls radius and weight guess.
    optim = @(mapv) syncirmin(datt1, mapv); %Creates an optimization problem based on our data section and our map varaible
    [x0,~,~,~,~] = lsqnonlin(optim, x0, [rkerd1/2,rkerd1/2,rkerd1 - 3,0], [rkerd1*3/2,rkerd1*3/2,rkerd1 + 3,rkerd1/5], options); %Optimizes x0
    x0(1) = x0(1) + map(l) - round(rt/2); %Readjusts row and column locations
    x0(2) = x0(2) + map(np1+l) - round(ct/2);
    rectt = synthcirccreate(dat, [x0(1) + rkerd2 + ap, x0(2) + rkerd2 + ap, x0(3:4)], 1, 0); %Creates a single circle at the proper location
    dat = dat - rectt; %%%%% Subtracts the calculated circle from the data to minimize overlaps. If you want overlap, just comment this line out. %%%%%
    x0(3) = x0(3);
    xv1 = [xv1;x0]; %Takes x0 and stores it in the vector
end

dim2 = 65/100;
for l = 1:np2
    datt2 = dat(round(map(4*np1 + l) - ceil(dim2*rkerd2) + rkerd2 + ap): ...
        round(map(4*np1 + l) + ceil(rkerd2*dim2) + rkerd2 + ap), ...
        (round(map(4*np1 + np2 +l) - ceil(dim2*rkerd2) + rkerd2 + ap)): ...
        round(map(4*np1 + np2 +l) + ceil(dim2*rkerd2) + rkerd2 + ap));
    [rt,ct] = size(datt2);
    x0 = [ceil(rt/2), ceil(ct/2), map(4*np1 + 2*np2 + l), map(4*np1 + 3*np2 + l)]; %Pulls out row and column guess
    optim = @(mapv) syncirmin(datt2, mapv);
    [x0,~,~,~,~] = lsqnonlin(optim, x0, [rkerd2/2,rkerd2/2,rkerd2 - 3,0], [rkerd2*3/2,rkerd2*3/2,rkerd2 + 3,rkerd1/5], options);
    x0(1) = x0(1) + map(4*np1 + l) - round(rt/2);
    x0(2) = x0(2) + map(4*np1 + np2 + l) - round(ct/2);
    rectt = synthcirccreate(dat, [x0(1) + rkerd2 + ap, x0(2) + rkerd2 + ap, x0(3:4)], 0, 1);
    dat = dat - rectt; %%%%% See above comment on this command %%%%%
    x0(3) = x0(3);
    xv2 = [xv2;x0];
end 

xv1 = reshape(xv1,[1,4*length(xv1)]); %Reshapes the first vector
xv2 = reshape(xv2,[1,4*length(xv2)]); %Reshapes the second vector

xv = [xv1,xv2]; %Combines vector for one output.


end

