function [ Rsmall, Rlarge, Nlarge, Nsmall ] = PeGSConvDebug( img, Rsmallr, Nsmall, Rlarger, Nlarge )
% This is the wrapper for the discFind function. It impletements each
% function in the correct manner and order while also displaying small
% progress updates. Output is a data matrix that gives the row, column,
% radius, and weight of circles.
%
%   Input:  img - The image
%
%           Rsmall - Radius of the smaller particle
%
%           Nsmall - Number of smaller particles present
%
%           Rlarge - Radius of the larger particle
%
%           Nlarge - Number of larger particles present
%
%   Output: Values to be used for the rest of the data.
%

if length(Rsmallr) == 2
    Rsmallr = Rsmallr(1):Rsmallr(2);
end
if length(Rlarger) == 2
    Rlarger = Rlarger(1):Rlarger(2);
end

lc = false;
while lc ~= true
    fprintf('\nThis procedure will generate %d figures based on your input and may take up to %d seconds (on Eno). ',length(Rsmallr)*length(Rlarger),40+round(2.2*length(Rsmallr)*length(Rlarger)));
    choice = input('Would you like to continue? [Y/N]\n','s');
    choice = upper(choice);
    if strcmp(choice,'Y') || strcmp(choice,'N')
        lc = true;
    else
        fprintf('\nInvalid Input')
    end
end
if choice == 'N'
    return;
end

[imrf, imgf] = fftfilt(img); %Filters the image

lc = 3;

fig = figure(1);
imshow(imrf + imgf/4)
axis equal

fprintf('\nIMPORTANT: Use the data courser to check the value inside a disc. If it is not positive, the image SHOULD be negated.')
while lc == 3
    nc = input('\nShould the image be negated? [Y/N]\n','s');
    nc = upper(nc);
    if strcmp(nc,'Y')
        lc = true;
        nc = true;
    elseif strcmp(nc,'N')
        lc = false;
        nc = false;
    else
        fprintf('\nInvalid Input')
    end
end

close all

for l = 1:5
    imrf = im2double(imrf);
    imrf(imrf > 3/4 * max(imrf)) = 2/3*max(imrf(:));
end
imgf = im2double(imgf);
im = imrf + imgf/(1/5*max(imgf(:)));
im = imadjust((im),stretchlim(im));
im = im + 1/4*imrf;
im(im < -10) = im(im < -10)/min(im(im < -10));
if nc
    im = - im;
end

runagain = true;
while runagain
    for l = 1:length(Rsmallr)
        for i = 1:length(Rlarger)
            kerd1=Rsmallr(l)*2;
            kerd2=Rlarger(i)*2;
            
            ker1=kerncreate(0,kerd1+2,kerd1/2);
            ker2=kerncreate(0,kerd2+2,kerd2/2); %Creates each kernel
            
            rkerd1 = kerd1*1.1;
            rkerd2 = kerd2*1.1; %Estimates the real kernel diameters
            
            imn1=kernfind(im,ker1); %Makes some nice peaks where small circles are
            imn1 = imn1./max(imn1(:));
            
            imn2=kernfind(im,ker2); %Makes some nice peaks where circles are
            imn2 = imn2./max(imn2(:)); %Puts image on a more workable scale
            
            imn1(1:kerd1,:) = 0;
            imn1(:,1:kerd1) = 0;
            imn2(1:kerd1,:) = 0;
            imn2(:,1:kerd1) = 0; %Unfortunately, edge's do NOT play nice with fft. So, we just sorta pad them out :^)
            
            [cv1,cv2] = max2D(imn1, imn2, kerd1, kerd2, Nsmall, Nlarge); %Finds center locations to pixel accurate locations
            
            cv1(:,2) = cv1(:,2);
            cv1(:,1) = cv1(:,1);
            cv2(:,2) = cv2(:,2);
            cv2(:,1) = cv2(:,1);
            
            figure
            imshow(img)
            colormap gray
            title(sprintf('Original Image with Circle Overlays using Rsmall = %d and Rlarge = %d',Rsmallr(l),Rlarger(i)))
            axis equal
            hold on
            
            viscircles([cv1(:,2),cv1(:,1)],rkerd1/2*ones([1,length(cv1)]),'edgecolor','g');
            viscircles([cv2(:,2),cv2(:,1)],rkerd2/2*ones([1,length(cv2)]),'edgecolor','m');

        end
    end
    
    lc = false;
    while lc == false
        runagain = upper(input('Run again with different parameters? (Can edit Nsmall, Nlarge, Rsmallr, and Rlarger) [Y/N]\n','s'));
        if strcmp(runagain,'Y')
            lc = true;
            runagain = true;
            close all
        elseif strcmp(runagain,'N')
            lc = true;
            runagain = false;
        else
            fprintf('\nInvalid Input')
        end
    end
    
    lc = false;
    if runagain
        while lc ~= true
            c = input(sprintf('New input for Nsmall? (If no change is desired, simply hit enter. If not, enter only integer values.)\nPrevious value was %d\n',Nsmall));
            if isempty(c)
                lc = true;
            elseif rem(c,1) == 0
                Nsmall = c;
                lc = true;
            else
                fprintf('Invalid input\n')
            end
        end
        
        lc = false;
        c = [];
        while lc ~= true
            c = input(sprintf('\nNew input for Nlarge? (If no change is desired, simply hit enter. If not, enter only integer values.)\nPrevious value was %d\n',Nlarge));
            if isempty(c)
                lc = true;
            elseif rem(c,1) == 0
                Nlarge = c;
                lc = true;
            else
                fprintf('Invalid input\n')
            end
        end
        
        lc = false;
        c = [];
        while lc ~= true
            c = input(sprintf('\nNew input for Rsmallr? (If no change is desired, simply hit enter. If not, enter lower and upper bounds as [lower, upper].)\nPrevious values where [%d %d]\n',min(Rsmallr),max(Rsmallr)));
            if isempty(c)
                lc = true;
            elseif length(c) == 2
                Rsmallr = c(1):c(2);
                lc = true;
            else
                fprintf('Invalid input\n')
            end
        end
        
        lc = false;
        c = [];
        while lc ~= true
            c = input(sprintf('\nNew input for Rlarger? (If no change is desired, simply hit enter. If not, enter lower and upper bounds as [lower, upper].)\nPrevious values where [%d %d]\n',min(Rlarger),max(Rlarger))');
            if isempty(c)
                lc = true;
            elseif length(c) == 2
                Rlarger = c(1):c(2);
                lc = true;
            else
                fprintf('Invalid input\n')
            end
        end
    end
end

fprintf('If being used with PeGS, choose the smallest values for both radii (i.e if choices of 31 and 32 are both optimal, choose 31).\n')

if length(Rsmallr) > 1
    lc = false;
    while lc == false
        Rsmall = input('Which value of Rsmall produces the best fit image?\n');
        if rem(Rsmall,1) ~= 0
            fprintf('\nInvalid Input.\n')
        else
            lc = true;
        end
    end
else
    Rsmall = Rsmallr;
end

if length(Rlarger) > 1
    lc = false;
    while lc == false
        Rlarge = input('\nWhich value of Rlarge produces the best fit image?\n');
        if rem(Rlarge,1) ~= 0
            fprintf('\nInvalid Input.\n')
        else
            lc = true;
        end
    end
else
    Rlarge = Rlarger;
end
close all
end