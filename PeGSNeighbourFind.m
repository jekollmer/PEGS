function [particle] = PeGSNeighbourFind(Gimg, contactG2Threshold, dtol, CR, verbose, particle)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

N = length(particle);

xmat = zeros([N,1]);
ymat = zeros([N,1]); % Preallocation
rmat = zeros([N,1]);

for l = 1:N
    xmat(l) = particle(l).x;
    ymat(l) = particle(l).y; %Pulls data from particle structure
    rmat(l) = particle(l).r;
end

rmats = rmat; %Saves our radius matrix for later

dmat = pdist2([xmat,ymat],[xmat,ymat]); %Creates a distance matrix for particle center locations
rmat = rmat + rmat'; %Makes a combination of radii for each particle

friendmat = dmat < (rmat + dtol) & dmat~=0; %Logical "friend" matrix

friendmat = triu(friendmat); %Only examine the upper triangle portion (no repeats)
[f1, f2] = find(friendmat == 1); %Creates an index of particles that are considered touching

mask = abs(-CR:CR);
mask = mask.^2 + mask.^2';
mask = double(sqrt(mask) <= CR-1);

for l = 1:length(f1)
    x1 = particle(f1(l)).x;
    y1 = particle(f1(l)).y;
    r1 = particle(f1(l)).r;
    x2 = particle(f2(l)).x;
    y2 = particle(f2(l)).y;
    r2 = particle(f2(l)).r;
    
    contactXp1 = x1 + (r1 - CR) * cos(atan2(y2-y1,x2-x1));
    contactYp1 = y1 + (r1 - CR) * sin(atan2(y2-y1,x2-x1));
    
    contactXp2 = x1 + (r1 + CR + dmat(f1(l),f2(l)) - rmat(f1(l),f2(l))) * cos(atan2(y2-y1,x2-x1));
    contactYp2 = y1 + (r1 + CR + dmat(f1(l),f2(l)) - rmat(f1(l),f2(l))) * sin(atan2(y2-y1,x2-x1));
    
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
            display(['contact found between particles ',num2str(f1(l)),' and ',num2str(f2(l))]);
            viscircles([contactXp1; contactYp1]', CR,'EdgeColor','w');
            viscircles([contactXp2; contactYp2]', CR,'EdgeColor','w');
            text(contactXp1, contactYp1,num2str(contactG2p1),'Color','w');
            %drawnow;
        end
        %this is a valid contact, remember it
        particle(f1(l)).z= particle(f1(l)).z +1; %increase coordination number
        particle(f1(l)).contactG2s(particle(f1(l)).z)=contactG2p1; %remember the g2 value of the current contact area
        particle(f1(l)).contactIs(particle(f1(l)).z)=contactIp1;
        particle(f1(l)).neighbours(particle(f1(l)).z) = f2(l); %particle m is now noted as a neigbour in the particle l datastructure
        particle(f1(l)).betas(particle(f1(l)).z) = atan2(y2-y1,x2-x1); %the contact angle to particle m is now noted in the particle l datastructure
        particle(f2(l)).z= particle(f2(l)).z +1; %increase coordination number
        particle(f2(l)).contactG2s(particle(f2(l)).z)=contactG2p2; %remember the g2 value of the current contact area
        particle(f2(l)).contactIs(particle(f2(l)).z)=contactIp2;
        particle(f2(l)).neighbours(particle(f2(l)).z) = f1(l); %particle m is now noted as a neigbour in the particle l datastructure
        particle(f2(l)).betas(particle(f2(l)).z) = atan2(y1-y2,x1-x2);
    end
end

%Check if any of the walls is a neighbour as well
%TODO: Do this in a less verbose manner, I remember there was a function
%somewhere that tells if somethig is inside a polygon or not

circs = [ymat, xmat, rmats]; %Makes a circs matrix from old matrices

rightwall = max(circs(:,2) + circs(:,3));
leftwall = min(circs(:,2) - circs(:,3)); %Finds our theorhetical wall locations
topwall = min(circs(:,1) - circs(:,3));
bottomwall = max(circs(:,1) + circs(:,3));

rwi = find(circs(:,2) + circs(:,3) + dtol*1.5 >= rightwall);
lwi = find(circs(:,2) - circs(:,3) - dtol*1.5 <= leftwall); %Indexes based on particles that would be considered to be touching the wall
bwi = find(circs(:,1) + circs(:,3) + dtol*1.5 >= bottomwall);
twi = find(circs(:,1) - circs(:,3) - dtol*1.5 <= topwall);

for l = 1:length(lwi) %Runs through each index to check for contacts via gradients
    x = circs(lwi(l),2);
    y = circs(lwi(l),1);
    r = round(circs(lwi(l),3));
    
    contactX = x-(r-CR);
    contactY = y;
    
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
        particle(lwi(l)).z= particle(lwi(l)).z +1; %increase coordination number
        particle(lwi(l)).contactG2s(particle(lwi(l)).z)=contactG2;
        particle(lwi(l)).contactIs(particle(lwi(l)).z)=cI;
        particle(lwi(l)).neighbours(particle(lwi(l)).z) = -1; %the wall is now noted as a neigbour in the particle l datastructure
        particle(lwi(l)).betas(particle(lwi(l)).z) = pi; %the contact angle to the wall is now noted in the particle l datastructure
    end
end

for l = 1:length(rwi)
    x = circs(rwi(l),2);
    y = circs(rwi(l),1);
    r = round(circs(rwi(l),3));
    
    contactX = x+(r-CR);
    contactY = y;
    
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
        particle(rwi(l)).z= particle(rwi(l)).z +1; %increase coordination number
        particle(rwi(l)).contactG2s(particle(rwi(l)).z)=contactG2;
        particle(rwi(l)).contactIs(particle(rwi(l)).z)=cI;
        particle(rwi(l)).neighbours(particle(rwi(l)).z) = -2; %the wall is now noted as a neigbour in the particle l datastructure
        particle(rwi(l)).betas(particle(rwi(l)).z) = 0; %the contact angle to the wall is now noted in the particle l datastructure
    end
end

for l = 1:length(twi)
    x = circs(twi(l),2);
    y = circs(twi(l),1);
    r = round(circs(twi(l),3));
    
    contactX = x;
    contactY = y-(r-CR);
    
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
        particle(twi(l)).z= particle(twi(l)).z +1; %increase coordination number
        particle(twi(l)).contactG2s(particle(twi(l)).z)=contactG2;
        particle(twi(l)).contactIs(particle(twi(l)).z)=cI;
        particle(twi(l)).neighbours(particle(twi(l)).z) = -3; %the wall is now noted as a neigbour in the particle l datastructure
        particle(twi(l)).betas(particle(twi(l)).z) = -pi/2; %the contact angle to the wall is now noted in the particle l datastructure
    end
end

for l = 1:length(bwi)
    x = circs(bwi(l),2);
    y = circs(bwi(l),1);
    r = round(circs(bwi(l),3));
    
    contactX = x;
    contactY = y+(r-CR);
    
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
        particle(bwi(l)).z= particle(bwi(l)).z +1; %increase coordination number
        particle(bwi(l)).contactG2s(particle(bwi(l)).z)=contactG2;
        particle(bwi(l)).contactIs(particle(bwi(l)).z)=cI;
        particle(bwi(l)).neighbours(particle(bwi(l)).z) = -4; %the wall is now noted as a neigbour in the particle l datastructure
        particle(bwi(l)).betas(particle(bwi(l)).z) = pi/2; %the contact angle to the wall is now noted in the particle l datastructure
    end
end

end

