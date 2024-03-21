% Jonathan's Photoelastic Disk Solver
%
% Takes input from preprocessing script joDiskPrep.m as of 2016/09/27
% inspired from peDiskSolve by James Puckett (Phd-Thesis 2012)
% http://nile.physics.ncsu.edu
%
% If you use this sovler please cite the follwoing paper
% K.E. Daniels, J. E. Kollmer & J. G. Puckett, "Photoelastic force measurements in granular materials", Rev. Sci. Inst. (201X)
% DOI: XXXXXX
%
% The generation of the synthetic force images is parallelized 
% in a way that each row of the ouput image can be its own worker/thread
% it is usually limited by your number of processing cores though. 
%
% Running this script may take a long time (up to one minute per particle
% on a 2016 desktop computer) so it makes sense to send this off to a high
% performance computer. 
%
% Depending on your experimental data you can play around with the fit
% options, which might speed up processing considerabely.
%
% Last edit on 2016/09/28 by Jonathan Kollmer (jekollme@ncsu.edu)

close all %housekeeping
clear all %housekeeping

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%User defined values:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%which files are we processing ?
directory = 'DATA/test/';
files = dir([directory, 'Step09*_preprocessing.mat']); 

%how much of the particle diameter is used to fit the synthetic image 
%(1 = use everything). Change this parameter only if the fit doesn't work 
%correctly because of imaging artifacts at the particle boundary.
maskradius = 0.96;% 
scaling = 0.5; %scale the image by this factor before doing the fit

%do we want to see each particle on screen while it is fitted ?
verbose = false; 

%fit options: play around with these, depending on the quality of your
%images and the accuracy you want your result to have. 
%Setting a good TolFun value can considerabely speed up your processing.
fitoptions = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt','MaxIter',100,'MaxFunEvals',400,'TolFun',0.01,'Display','final-detailed');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%There should be no need for user input below this line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
nFrames = length(files); %how many files are we processing ?


for frame = 1:nFrames %loop over these frames 
    fileName = [directory,files(frame).name];
    data = load(fileName);
    particle = data.particle;
    particle = disksolve(particle, maskradius, scaling, fileName, verbose, fitoptions);
end
