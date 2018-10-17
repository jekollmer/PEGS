% Particle Detector, Neighbour finder and Contact validator that creates input files for my Photoelastic Disk Solver
% Particle Detection and Neigbour Finding Adapted from my Earlier Script (joCentersNonMonodisperse.m as of 2016/05/03)
% Photoelastic Disk Solver inspired from peDiskSolve by James Puckett (Phd-Thesis 2012) http://nile.physics.ncsu.edu

% If you use this please cite the follwoing paper
% K.E. Daniels, J. E. Kollmer & J. G. Puckett, "Photoelastic force measurements in granular materials", Rev. Sci. Inst. (201X)
% DOI: XXXXXX

% last edit on 2018/08/09 by Joshua Miller (jsmille9@ncsu.edu)

close all % Housekeeping
clear all % Housekeeping

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           User defined values                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pxPerMeter = 0.0077 / 74;
verbose = false; %Generates lots of plots showing results

directory = 'DATA/test/';
files = dir([directory, 'Step09*.jpg']); %Which files are we processing?
%files = dir('Centers*0001.txt'); %Alternatively, centers files can be loaded. This requires that both particle detections be flagged false however.
nFrames = length(files); %How many files are we processing ?

% Hough Transform Values

doParticleDetectionH = false; %Detect particles using Hough Transform?
HoughDebug = false; %Debugs Hough Sensitivities so particles are found "better"

DS = 0.0025; % How much should we adjust sensitivity if wrong number of particles are found
RlargeH = [70 80]; %What radius (in pixels) range do we expect for the large discs?
RsmallH = [45 55]; %What radius (in pixels) range do we expect for the small discs?
SL = 0.95; %Sensitivity of the Hough Transform disc detetcor, exact value is Voodo magic...
SS = 0.89; %Sensitivity of the Hough Transform disc detetcor, exact value is Voodo magic...

NsmallH = 15; %Number of small discs. Only used in Hough Debug.
NlargeH = 14; %Number of large discs. Only used in Hough Debug.

% Convolution Method Values

doParticleDetectionC = true; %Detect particles using convolution method?
ConvDebug = false;

RlargeC = 31; %What radius (in pixels) do we expect for the large discs?
RsmallC = 21; %What radius (in pixels) do we expect for the small discs? 
%Note: The above can be input in a range ONLY if the ConvDebug is set to
%Note: true. Otherwise, the program needs the radius that works.
NsmallC = 507; %Number of small discs. Needed for Convolution.
NlargeC = 273; %Number of large discs. Needed for Convolution.

% Neighbour Finding Values

findNeighbours = true;

fsigma = 100; %photoelastic stress coefficient
g2cal = 100; %Calibration Value for the g^2 method, can be computed by joG2cal.m
dtol = 10; % How far away can the outlines of 2 particles be to still be considered Neighbours

contactG2Threshold = 0.5; %sum of g2 in a contact area larger than this determines a valid contact
CR = 10; %radius around a contactact point that is checked for contact validation


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 User Input Not Required Below This Line                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if HoughDebug
    
    imageFile = [directory,files(1).name]; %input filename
    img = imread(imageFile);
    Rimg = img(:,:,1);
    [SL, SS] = PeGSHoughDebug(Rimg, RlargeH, SL, RsmallH, SS, DS, NlargeH, NsmallH);
    
elseif ConvDebug
    
    imageFile = [directory,files(1).name]; %input filename
    img = imread(imageFile);
    [RsmallC, RlargeC, NlargeC, NsmallC] = PeGSConvDebug(img, RsmallC, NsmallC, RlargeC, NlargeC);
    
end


for frame = 1:nFrames %Loops for total number of images
    
    if doParticleDetectionH || doParticleDetectionC
        
        imageFile = [directory,files(frame).name]; %input filename
        img = imread(imageFile); %read a color image that has particles in red and forces in green channel
        Rimg = img(:,:,1); %particle image
        Gimg = img(:,:,2); %force image
        
    else
        
        centersfile = [directory, files(frame).name]; %input filename
        gImgFile = [directory, 'Img',centersfile(8:end-3),'jpg'];  %adjusted force image filename
        rImgFile = [directory, 'Img',centersfile(8:end-3),'jpg'];  %adjusted force image filename
        Rimg = imread(rImgFile); %particle image
        Gimg = imread(gImgFile); %force image
        
    end
    
    Gimg = im2double(Gimg);
    Rimg = im2double(Rimg);
    Gimg = Gimg-0.5*Rimg;
    Gimg = Gimg.*(Gimg > 0);
    Gimg = imadjust(Gimg,stretchlim(Gimg));
    
    if (verbose)
        
        figure(1); %Draw the particle Image
        imshow(Rimg);
        
        figure(2); %Draw the Force Image
        imshow(Gimg);
        
    end
    
    if doParticleDetectionH
        
        particle = PeGSDiskFindH(Rimg, RlargeH, SL, RsmallH, SS, pxPerMeter, fsigma);
        
    elseif doParticleDetectionC
        
        particle = PeGSDiskFindC(img, RsmallC, NsmallC, RlargeC, NlargeC);
        
    else
        
        pData = dlmread(centersfile); %Read Position data from centers file
        N = size(pData,1);
        particle(1:N) = struct('id',0,'x',0,'y',0,'r',0,'rm',0,'color','','fsigma',0,'z',0,'f',0,'g2',0,'forces',[],'betas',[],'alphas',[],'neighbours',[],'contactG2s',[],'forceImage',[]);
        for n=1:N %Bookkeeping
            particle(n).id= n;
            particle(n).x = pData(n,2); %-xoffset;
            particle(n).y = pData(n,3); %-yoffset;
            particle(n).r = pData(n,4);
        end
        
    end
    
    N = length(particle);
    
    if(verbose)
        %add some information about the particles to the plots
        figure(1)
        for n=1:N
            viscircles([particle(n).x; particle(n).y]', particle(n).r,'EdgeColor',particle(n).color); %draw particle outline
            hold on
            plot(particle(n).x,particle(n).y,'rx'); %Mark particle centers
            text(particle(n).x,particle(n).y,num2str(particle(n).id),'Color','w');
        end
        figure(2)
        for n=1:N
            viscircles([particle(n).x; particle(n).y]', particle(n).r,'EdgeColor',particle(n).color); %draw particle outline
            hold on
            plot(particle(n).x,particle(n).y,'rx'); %Mark particle centers
            text(particle(n).x,particle(n).y,num2str(particle(n).id),'Color','w');
        end
        drawnow;
    end
    
    for n=1:N
        %create a circular mask
        % => Find a better way yo do this masking!
        r = particle(n).r;
        mask = abs(-r:r);
        mask = mask.^2 + mask.^2';
        mask1 = double(sqrt(mask) <= r);
        
        %This crops out a particle
        cropXstart = round(particle(n).x-r);
        cropXstop = round(particle(n).x-r)+ size(mask1,1)-1;
        cropYstart = round(particle(n).y-r);
        cropYstop = round(particle(n).y-r)+ size(mask1,2)-1;
        cimg = im2double(Gimg(cropYstart:cropYstop, cropXstart:cropXstop));
        particleImg = cimg.*mask1;
        particle(n).forceImage=particleImg;
        
        %create a circular mask with a radius that is one pixel smaller
        %for cropping out the relevant gradient

        mask2 = double(sqrt(mask) <= r-1);
        
        %Compute G^2 for each particle
        [gx,gy] = gradient(particleImg);
        g2 = (gx.^2 + gy.^2).*mask2;
        particle(n).g2 = sum(sum(g2));
        particle(n).f = particle(n).g2/g2cal;
    end
    
    if findNeighbours
        
        particle = PeGSNeighbourFind(Gimg, contactG2Threshold, dtol, CR, verbose, particle);
        
    end
    
    figure(3);
    imshow(Gimg); hold on
    for n = 1:N
        z = particle(n).z; %get particle coordination number
        if (z>0) %if the particle does have contacts
            for m = 1:z %for each contact
                %draw contact lines
                lineX(1)=particle(n).x;
                lineY(1)=particle(n).y;
                lineX(2) = lineX(1) + particle(n).r * cos(particle(n).betas(m));
                lineY(2) = lineY(1) + particle(n).r * sin(particle(n).betas(m));
                cX = lineX(1) + (particle(n).r-CR) * cos(particle(n).betas(m));
                cY = lineY(1) + (particle(n).r-CR) * sin(particle(n).betas(m));
                hold on; % Don't blow away the image.
                plot(lineX, lineY,'-y','LineWidth',2);hold on;
            end
        end
    end
    
    %Save what we got so far
    save([directory, files(frame).name(1:end-4),'_preprocessing.mat'],'particle');
    
end
