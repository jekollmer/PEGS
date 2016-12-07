% Particle Detector, Neighbour finder and Contact validator that creates input files for my Photoelastic Disk Solver
% Particle Detection and Neigbour Finding Adapted from my Earlier Script (joCentersNonMonodisperse.m as of 2016/05/03)
% Photoelastic Disk Solver inspired from peDiskSolve by James Puckett (Phd-Thesis 2012) http://nile.physics.ncsu.edu

% If you use this please cite the follwoing paper
% K.E. Daniels, J. E. Kollmer & J. G. Puckett, "Photoelastic force measurements in granular materials", Rev. Sci. Inst. (201X)
% DOI: XXXXXX

% last edit on 2016/09/28 by Jonathan Kollmer (jekollme@ncsu.edu)

close all % Housekeeping
clear all % Housekeeping

%User defined values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pxPerMeter = 0.0077 / 74 % pixel were measured to be 0.0077 m
%((1700-1659)/2)
%Particle information
Rlarge = [70 80]; %What radius (in pixels) range do we expect for the small discs?
Rsmall = [45 55]; %What radius (in pixels) range do we expect for the small discs?
SL = 0.95; % Sensitivity of the Hough Transform disc detetcor, exact value is Voodo magic...
SS = 0.89; % Sensitivity of the Hough Transform disc detetcor, exact value is Voodo magic...


fsigma = 100; %photoelastic stress coefficient
g2cal = 100; %Calibration Value for the g^2 method, can be computed by joG2cal.m 
dtol = 10; % How far away can the outlines of 2 particles be to still be considered Neighbours

%This is for contact validation:
contactG2Threshold = 0.5; %sum of g2 in a contact area larger than this determines a valid contact
CR = 10; %radius around a contactact point that is checked for contact validation

%Position of the confing walls
%TODO: Add automatic wall finding
leftwall=59; %pixel coordinate of the right wall %Particle ID -1
rightwall=805; %pixel coordinate of the right wall %Particle ID -2
topwall=112; %pixel coordinate of the right wall %Particle ID -3
bottomwall=790; %pixel coordinate of the right wall %Particle ID -4

% xoffset = 1000;
% yoffset = 700;
% xsize = 2470-xoffset;
% ysize = 2260-yoffset;


%Automatic Filename Generation for Batch Processing

%files = dir('Centers*0001.txt'); %which files are we processing ?
directory = 'DATA/test/';
files = dir([directory, 'Step09*.jpg']); %which files are we processing ?
nFrames = length(files); %how many files are we processing ?

doParticleDetection=true; %detect particles or use centers-file ?

verbose = true; %generate lots of plots showing the results

for frame = 1%:nFrames %loop over these cycles 

    close all
    clear particle
    
    if doParticleDetection
        imageFile = [directory,files(frame).name]; %input filename 
        img = imread(imageFile); %read a color image that has particles in red and forces in green channel
        %img = imcrop(img,[xoffset, yoffset, xsize, ysize]); %particle image
        Rimg = img(:,:,1); %particle image
        Gimg = img(:,:,2); %force image
        
        
        
    else
        centersfile = [directory, files(frame).name]; %input filename 
        gImgFile = [directory, 'Img',centersfile(8:end-3),'jpg'];  %adjusted force image filename 
        rImgFile = [directory, 'Img',centersfile(8:end-3),'jpg'];  %adjusted force image filename 
        Rimg = imread(rImgFile); %particle image
        Gimg = imread(gImgFile); %force image
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %There should be no need for user input below this line
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %read image data and crop out region of interest
    %Rimg = imcrop(Rimg,[xoffset, yoffset, xsize, ysize]); %particle image
    %Gimg = imcrop(Gimg,[xoffset, yoffset, xsize, ysize]); %force image

    %This is the Camera Image
    %template = im2double(particle(n).forceImage);
    
    Gimg = im2double(Gimg);
    Rimg = im2double(Rimg);
    
    Gimg = Gimg-0.5*Rimg;
    
    Gimg = Gimg.*(Gimg > 0); %fine tunes the image, this should happen in preprocessing!
    %Gimg = Gimg*3; %fine tunes the image, this should happen in preprocessing!
    Gimg = imadjust(Gimg,stretchlim(Gimg));
    %display the particle image and the force image
    if (verbose)
        figure(1);
        
        %Draw the Particle Image
        imshow(Rimg); hold on
        %Draw the Boundary
        rectangle('Position',[leftwall topwall rightwall-leftwall bottomwall-topwall],'EdgeColor','r');

        figure(2);
        %Draw the Force Image
        imshow(Gimg); hold on 
        %Draw the Boundary
        rectangle('Position',[leftwall topwall rightwall-leftwall bottomwall-topwall],'EdgeColor','r');
    end
    
    if doParticleDetection
        
        %Disc detection for large particles happens here
        [centers, radii, metric] = imfindcircles(Rimg,Rlarge,'ObjectPolarity','bright','Method','TwoStage','Sensitivity',SL);
        Nlarge=length(centers);
        particle(1:Nlarge) = struct('id',0,'x',0,'y',0,'r',0,'rm',0,'color','','fsigma',0,'z',0,'f',0,'g2',0,'forces',[],'betas',[],'alphas',[],'neighbours',[],'contactG2s',[],'forceImage',[]);
        
        %Bookkeeping
        for n=1:Nlarge
            particle(n).id= n;
            particle(n).x = centers(n,1);
            particle(n).y = centers(n,2);
            particle(n).r = radii(n);
        end
        
        %Disc detection for small particles happens here
        [centers, radii, metric] = imfindcircles(Rimg,Rsmall,'ObjectPolarity','bright','Method','TwoStage','Sensitivity',SS);
        Nsmall=length(centers);
        N = Nsmall + Nlarge;
        particle(Nlarge+1:N) = struct('id',0,'x',0,'y',0,'r',0,'rm',0,'color','','fsigma',0,'z',0,'f',0,'g2',0,'forces',[],'betas',[],'alphas',[],'neighbours',[],'contactG2s',[],'forceImage',[]);
        
        %Bookkeeping
        for n=1:Nsmall
            particle(Nlarge+n).id= Nlarge+n;
            particle(Nlarge+n).x = centers(n,1);
            particle(Nlarge+n).y = centers(n,2);
            particle(Nlarge+n).r = radii(n);
        end

        

    else
        pData = dlmread(centersfile); %Read Position data from centers file
        N = size(pData,1);
        particle(1:N) = struct('id',0,'x',0,'y',0,'r',0,'rm',0,'color','','fsigma',0,'z',0,'f',0,'g2',0,'forces',[],'betas',[],'alphas',[],'neighbours',[],'contactG2s',[],'forceImage',[]);
        for n=1:N%Bookkeeping
            particle(n).id= n;
            particle(n).x = pData(n,2);%-xoffset;
            particle(n).y = pData(n,3);%-yoffset;
            particle(n).r = pData(n,4);
        end
    end
    
    %Bookkeeping
    for n=1:N
        particle(n).rm = particle(n).r * pxPerMeter;
        particle(n).fsigma = fsigma;
        if particle(n).r < ( mean(Rlarge)+ mean(Rsmall) ) /2
            particle(n).color = 'b'; 
        else
            particle(n).color = 'r';
        end
    end

    if(verbose)
        %add some information about the particles to the plots
        figure(1)
        for n=1:N
         viscircles([particle(n).x; particle(n).y]', particle(n).r,'EdgeColor',particle(n).color); %draw particle outline
         hold on
         plot(particle(n).x,particle(n).y,'rx'); %mark particle centers
         text(particle(n).x,particle(n).y,num2str(particle(n).id),'Color','w');
        end
        figure(2)
        for n=1:N
         viscircles([particle(n).x; particle(n).y]', particle(n).r,'EdgeColor',particle(n).color); %draw particle outline
         hold on
         plot(particle(n).x,particle(n).y,'rx'); %mark particle centers
         text(particle(n).x,particle(n).y,num2str(particle(n).id),'Color','w');
        end
        drawnow;
    end
    
    %Crop out Particles and store the PE response of each particle in the
    %particle datastructure
    for n=1:N
        %create a circular mask
        % => Find a better way yo do this masking!
        px = particle(n).r*2;
        cx=round(px/2);cy=round(px/2);ix=px;iy=px;r=px/2;
        [x,y]=meshgrid(-(cx):(ix-cx),-(cy):(iy-cy)); %does this alswys guarantee the right size ?? (rounding)
        mask=double(((x.^2+y.^2)<=r^2));

        %This crops out a particle
        cropXstart = round(particle(n).x-r);
        cropXstop = round(particle(n).x-r)+ size(mask,1)-1;
        cropYstart = round(particle(n).y-r);
        cropYstop = round(particle(n).y-r)+ size(mask,2)-1;
        cimg = im2double(Gimg(cropYstart:cropYstop, cropXstart:cropXstop));
        particleImg = cimg.*mask;
        particle(n).forceImage=particleImg;

        %create a circular mask with a radius that is one pixel smaller
        %for cropping out the relevant gradient
        % => Find a better way yo do this masking!
        px = particle(n).r*2;
        cx=px/2;cy=px/2;ix=px;iy=px;r=px/2-1;
        [x,y]=meshgrid(-(cx):(ix-cx),-(cy):(iy-cy)); %does this alswys guarantee the right size ?? (rounding)
        mask=double(((x.^2+y.^2)<=r^2));

        %Compute G^2 for each particle
        [gx,gy] = gradient(particleImg);
        g2 = (gx.^2 + gy.^2).*mask;
        particle(n).g2 = sum(sum(g2));
        particle(n).f = particle(n).g2/g2cal;
    end


    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %HERE STARTS THE NEIGBOUR FINDING PART    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Find Particles that are Neighbours
    %TODO: This becomes somewhat inefficient for large N, there is a couple of
    %TODO: ways to improve on this, mainly by sorting particles onto a grid and then
    %TODO: only looking at neighbouring cells when comparing positions (optimal cell
    %TODO: size is then the radius of the largest particle
    for l = 1:N
        for m = 1:N
            if ( l ~= m ) 
                %make writing the code a little more easy by subsituting names
                x1 = particle(l).x;
                y1 = particle(l).y;
                x2 = particle(m).x;
                y2 = particle(m).y;
                r1 = particle(l).r;
                r2 = particle(m).r;

                surfaceDistance = sqrt((x2-x1)^2+(y2-y1)^2)-(particle(m).r+particle(l).r);
                if(surfaceDistance<dtol)
                    display(['taking a closer look at particles ',num2str(l),' and ', num2str(m)]);
                    %now check if the contact area has some pe signal
                    %first determine contact point
                    contactXp1 = x1 + (r1-CR) * cos(atan2(y2-y1,x2-x1));
                    contactYp1 = y1 + (r1-CR) * sin(atan2(y2-y1,x2-x1));

                    contactXp2 = x1 + (r1+CR+surfaceDistance) * cos(atan2(y2-y1,x2-x1));
                    contactYp2 = y1 + (r1+CR+surfaceDistance) * sin(atan2(y2-y1,x2-x1));

                    %create a circular mask with a radius that is one pixel smaller
                    % => Find a better way yo do this masking!
                    px = CR*2;
                    cx=px/2;cy=px/2;ix=px;iy=px;r=px/2-1;
                    [x,y]=meshgrid(-(cx):(ix-cx),-(cy):(iy-cy)); %does this alswys guarantee the right size ?? (rounding)
                    mask=double(((x.^2+y.^2)<=r^2));

                    %then look at g2 in a (small) region around the contact
                    contactImg = im2double(imcrop(Gimg,[contactXp1-CR contactYp1-CR CR*2 CR*2]));
                    contactImg = contactImg.*mask;

                    [gx,gy] = gradient(contactImg);
                    g2 = (gx.^2 + gy.^2);
                    contactG2p1 = sum(sum(g2));
                    contactIp1 = sum(sum(contactImg));

                    contactImg = im2double(imcrop(Gimg,[contactXp2-CR contactYp2-CR CR*2 CR*2]));
                    contactImg = contactImg.*mask;

                    [gx,gy] = gradient(contactImg);
                    g2 = (gx.^2 + gy.^2);
                    contactG2p2 = sum(sum(g2));
                    contactIp2 = sum(sum(contactImg));

                    %if we declare our contact valid
                    if(contactG2p1 > contactG2Threshold && contactG2p2 > contactG2Threshold)
                    %cI = sum(sum(contactImg)); %Use integrated intensity instead of g2             
                    %if(cI > contactIThreshold) %Use integrated intensity instead of g2      
                        %Plot contact area
                        if (verbose)
                            display(['contact found between particles ',num2str(l),' and ',num2str(m)]);
                            viscircles([contactXp1; contactYp1]', CR,'EdgeColor','w');
                            viscircles([contactXp2; contactYp2]', CR,'EdgeColor','w');
                            text(contactXp1, contactYp1,num2str(contactG2p1),'Color','w');
                            %drawnow;
                        end
                        %this is a valid contact, remember it                    
                        particle(l).z= particle(l).z +1; %increase coordination number
                        particle(l).contactG2s(particle(l).z)=contactG2p1; %remember the g2 value of the current contact area
                        particle(l).contactIs(particle(l).z)=contactIp1;
                        particle(l).neighbours(particle(l).z) = m; %particle m is now noted as a neigbour in the particle l datastructure
                        particle(l).betas(particle(l).z) = atan2(y2-y1,x2-x1); %the contact angle to particle m is now noted in the particle l datastructure
                    end
                end
            end
        end
    end


    %create a circular mask with a radius that is one pixel smaller
    % => Find a better way yo do this masking!
    px = CR*2;
    cx=px/2;cy=px/2;ix=px;iy=px;r=px/2-1;
    [x,y]=meshgrid(-(cx):(ix-cx),-(cy):(iy-cy)); %does this alswys guarantee the right size ?? (rounding)
    mask=double(((x.^2+y.^2)<=r^2));


    %Check if any of the walls is a neighbour as well
    %TODO: Do this in a less verbose manner, I remember there was a function
    %somewhere that tells if somethig is inside a polygon or not
    for l = 1:N
       if particle(l).x-particle(l).r < ( leftwall + dtol )
            %now check if the contact area has some pe signal
            %first determine contact point
            contactX = particle(l).x-(particle(l).r-CR);
            contactY = particle(l).y;
            %then look at g2 in a (small) region around the contact
            %contactImg = im2double(imcrop(Gimg,[contactX-r1/4 contactY-r1/4 r1/2 r1/2]));
            contactImg = im2double(imcrop(Gimg,[contactX-CR contactY-CR CR*2 CR*2]));
            contactImg = contactImg.*mask;

            [gx,gy] = gradient(contactImg);
            g2 = (gx.^2 + gy.^2);
            contactG2 = sum(sum(g2));

            if(contactG2 > contactG2Threshold)
            cI = sum(sum(contactImg));             
            %if(cI > contactIThreshold)
                %this is a valid contact, remember it
                if(verbose)
                    text(contactX,contactY,num2str(contactG2),'Color','w');
                    viscircles([contactX; contactY]', CR,'EdgeColor','w');
                end
                particle(l).z= particle(l).z +1; %increase coordination number
                particle(l).contactG2s(particle(l).z)=contactG2;
                particle(l).contactIs(particle(l).z)=cI;
                particle(l).neighbours(particle(l).z) = -1; %the wall is now noted as a neigbour in the particle l datastructure
                particle(l).betas(particle(l).z) = pi; %the contact angle to the wall is now noted in the particle l datastructure
            end
       end

       if particle(l).x+particle(l).r > ( rightwall - dtol )
            %now check if the contact area has some pe signal
            %first determine contact point
            contactX = particle(l).x+(particle(l).r-CR);
            contactY = particle(l).y;
            %then look at g2 in a (small) region around the contact
            %contactImg = im2double(imcrop(Gimg,[contactX-r1/4 contactY-r1/4 r1/2 r1/2]));
            contactImg = im2double(imcrop(Gimg,[contactX-CR contactY-CR CR*2 CR*2]));
            contactImg = contactImg.*mask;

            [gx,gy] = gradient(contactImg);
            g2 = (gx.^2 + gy.^2);
            contactG2 = sum(sum(g2));
            if(contactG2 > contactG2Threshold)
            cI = sum(sum(contactImg));             
            %if(cI > contactIThreshold)
                %this is a valid contact, remember it
                if(verbose)
                    text(contactX,contactY,num2str(contactG2),'Color','w');
                    viscircles([contactX; contactY]', CR,'EdgeColor','w');
                end
                particle(l).z= particle(l).z +1; %increase coordination number
                particle(l).contactG2s(particle(l).z)=contactG2;
                particle(l).contactIs(particle(l).z)=cI;
                particle(l).neighbours(particle(l).z) = -2; %the wall is now noted as a neigbour in the particle l datastructure
                particle(l).betas(particle(l).z) = 0; %the contact angle to the wall is now noted in the particle l datastructure
            end
       end

       if particle(l).y-particle(l).r < ( topwall + dtol )
            %now check if the contact area has some pe signal
            %first determine contact point
            contactX = particle(l).x;
            contactY = particle(l).y-(particle(l).r-CR);
            %then look at g2 in a (small) region around the contact
            %contactImg = im2double(imcrop(Gimg,[contactX-r1/4 contactY-r1/4 r1/2 r1/2]));
            contactImg = im2double(imcrop(Gimg,[contactX-CR contactY-CR CR*2 CR*2]));
            contactImg = contactImg.*mask;

            [gx,gy] = gradient(contactImg);
            g2 = (gx.^2 + gy.^2);
            contactG2 = sum(sum(g2));
            if(contactG2 > contactG2Threshold)
            cI = sum(sum(contactImg));             
            %if(cI > contactIThreshold)
                %this is a valid contact, remember it
                if(verbose)
                    text(contactX,contactY,num2str(contactG2),'Color','w');
                    viscircles([contactX; contactY]', CR,'EdgeColor','w');
                end
                particle(l).z= particle(l).z +1; %increase coordination number
                particle(l).contactG2s(particle(l).z)=contactG2;
                particle(l).contactIs(particle(l).z)=cI;
                particle(l).neighbours(particle(l).z) = -3; %the wall is now noted as a neigbour in the particle l datastructure
                particle(l).betas(particle(l).z) = -pi/2; %the contact angle to the wall is now noted in the particle l datastructure
            end
       end

       if particle(l).y+particle(l).r > ( bottomwall - dtol )
            %now check if the contact area has some pe signal
            %first determine contact point
            contactX = particle(l).x;
            contactY = particle(l).y+(particle(l).r-CR);
            %then look at g2 in a (small) region around the contact
            %contactImg = im2double(imcrop(Gimg,[contactX-r1/4 contactY-r1/4 r1/2 r1/2]));
            contactImg = im2double(imcrop(Gimg,[contactX-CR contactY-CR CR*2 CR*2]));
            contactImg = contactImg.*mask;

            [gx,gy] = gradient(contactImg);
            g2 = (gx.^2 + gy.^2);
            contactG2 = sum(sum(g2));
            if(contactG2 > contactG2Threshold)
            cI = sum(sum(contactImg));             
            %if(cI > contactIThreshold)
                %this is a valid contact, remember it
                if(verbose)
                    text(contactX,contactY,num2str(contactG2),'Color','w');
                    viscircles([contactX; contactY]', CR,'EdgeColor','w');
                end
                particle(l).z= particle(l).z +1; %increase coordination number
                particle(l).contactG2s(particle(l).z)=contactG2;
                particle(l).contactIs(particle(l).z)=cI;
                particle(l).neighbours(particle(l).z) = -4; %the wall is now noted as a neigbour in the particle l datastructure
                particle(l).betas(particle(l).z) =  pi/2; %the contact angle to the wall is now noted in the particle l datastructure
            end
       end

    end

    %HERE ENDS THE NEIGBOUR FINDING PART
    
    %plot the preliminary contact network
    %if(verbose)
        figure(3);
        imshow(Gimg); hold on 
        for n = 1:N
            z = particle(n).z; %get particle coordination number
            %viscircles([particle(n).x; particle(n).y]', particle(n).r,'EdgeColor',particle(n).color); %draw particle outline
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
                    %viscircles([cX; cY]', CR,'EdgeColor','w');
                end
            end
        end
    %end
    
    %Save what we got so far
    %save([directory, files(frame).name(1:end-4),'_fsigma',num2str(fsigma),'_preprocessing.mat'],'particle');
    save([directory, files(frame).name(1:end-4),'_preprocessing.mat'],'particle');

end