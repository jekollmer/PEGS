function [ cv1, cv2 ] = max2D( dat1, dat2, kerd1, kerd2, np1, np2 )     %MAY NEED OPTIMIZATION
% max2D will intake a single data matrix (a 2D array) and "scan" it for
% local maxima and attempt to avoid any false maxima by taking away
% partcles too close to one another.
%
%   Input:  dat1 - Convolution of the image and the first kernel
%
%           dat2 - Convolution of the image and the second kernel
%
%           kerd1 - Diameter of the first kernel
%
%           kerd2 - Diameter of the second kernel
%
%           np1 - Number of particles of size 1
%
%           np2 - Number of particles of size 2
%
%   Output: ov1 - Vector of center locations for particles of size 1
%
%           ov2 - Vector of center locations for particles of size 2

edg = 2; %We need to assign an edge length of >1 so we do not attempt to "escape" our image.

datav1 = mean(dat1(:)); %Averages the data
dat1 = padarray(dat1, [edg,edg], 0, 'both'); %Pads the data with zero around the very small edge
datav2 = mean(dat2(:));
dat2 = padarray(dat2, [edg,edg], 0, 'both');

sd1 = size(dat1);
sd2 = size(dat2); %Sizes both data sets

[x1,y1] = find(dat1(edg:sd1(1)-edg,edg:sd1(2)-edg)); %Gives an appropriate x and y data set based on edge size
[x2,y2] = find(dat2(edg:sd2(1)-edg,edg:sd2(2)-edg));

y1 = y1+edg-1;
x1 = x1+edg-1; %Shifts x and y to adjust for the edge
y2 = y2+edg-1;
x2 = x2+edg-1;

ths = 0;
xv1 = [];
yv1 = [];
I1 = [];
xv2 = xv1;
yv2 = yv1;
I2 = I1; %Assigns values that are used to store max locations.

for l=1:length(y2)
    if (dat1(x1(l), y1(l)) >= (1.000+ths)*dat1(x1(l)-1, y1(l)-1 )) && ...
            (dat1(x1(l), y1(l)) > (1.000+ths)*dat1(x1(l)-1, y1(l))) && ...
            (dat1(x1(l), y1(l)) >= (1.000+ths)*dat1(x1(l)-1, y1(l)+1)) && ...
            (dat1(x1(l), y1(l)) > (1.000+ths)*dat1(x1(l), y1(l)-1)) && ...
            (dat1(x1(l), y1(l)) > (1.000+ths)*dat1(x1(l), y1(l)+1)) && ...
            (dat1(x1(l), y1(l)) >= (1.000+ths)*dat1(x1(l)+1, y1(l)-1)) && ...
            (dat1(x1(l), y1(l)) > (1.000+ths)*dat1(x1(l)+1, y1(l))) && ...
            (dat1(x1(l), y1(l)) >= (1.000+ths)*dat1(x1(l)+1, y1(l)+1)) && ...
            (dat1(x1(l), y1(l)) >= .1) %Checks all the sorrouding points to assure the point x1(l) y1(l) is larger than them. Additionally, checks that it is at least .3 times the mean values of the image.
        I1 = [I1; dat1(x1(l), y1(l))];
        xv1 = [xv1; x1(l)];
        yv1 = [yv1; y1(l)]; %Assigns found maximums and intensity of peak
    end
    if (dat2(x2(l), y2(l)) >= (1.000+ths)*dat2(x2(l)-1, y2(l)-1 )) && ...
            (dat2(x2(l), y2(l)) > (1.000+ths)*dat2(x2(l)-1, y2(l))) && ...
            (dat2(x2(l), y2(l)) >= (1.000+ths)*dat2(x2(l)-1, y2(l)+1)) && ...
            (dat2(x2(l), y2(l)) > (1.000+ths)*dat2(x2(l), y2(l)-1)) && ...
            (dat2(x2(l), y2(l)) > (1.000+ths)*dat2(x2(l), y2(l)+1)) && ...
            (dat2(x2(l), y2(l)) >= (1.000+ths)*dat2(x2(l)+1, y2(l)-1)) && ...
            (dat2(x2(l), y2(l)) > (1.000+ths)*dat2(x2(l)+1, y2(l))) && ...
            (dat2(x2(l), y2(l)) >= (1.000+ths)*dat2(x2(l)+1, y2(l)+1)) && ...
            (dat2(x2(l), y2(l)) >= .1)
        I2 = [I2; dat2(x2(l), y2(l))];
        xv2 = [xv2; x2(l)];
        yv2 = [yv2; y2(l)]; %Does the same as above but for image 2
    end
end

cv1 = [I1, xv1, yv1]; %Intensity, Row Location, Column Location
cv2 = [I2, xv2, yv2]; %Creates compound vectors 1 and 2

dmatc = pdist2([cv1(:,2),cv1(:,3)],[cv2(:,2),cv2(:,3)]); %The points going downwards are from cv1. 
inc = dmatc < (kerd2/2+kerd1/4); %Creates a distance array that contains the distance between every particle of differing sizes.
[rc, cc] = find(inc == 1); %Finds where the distance between two particles is less than the largest radius

for l = 1:length(rc)
    Ic1 = cv1((rc(l)),1);
    Ic2 = cv2((cc(l)),1); %Compares the intensity at both locations
    if Ic1 > Ic2
        cv2((cc(l)),1) = 0; %If the first intensity is higher, than it is taken as the true max. By setting the second equal to zero, we remove it later on.
    elseif Ic2 > Ic1
        cv1((rc(l)),1) = 0; %Does the same as above but for if the second intensity is larger
    end
end

dmat1 = pdist2([cv1(:,2),cv1(:,3)],[cv1(:,2),cv1(:,3)]);
in1 = dmat1 < (kerd1/2 ) & dmat1 ~= 0; %Creates a distance matrix for particles of size 1 to assure they do not have an accidental double peak.
in1 = triu(in1);
[rc1, cc1] = find(in1 == 1);

dmat2 = pdist2([cv2(:,2),cv2(:,3)],[cv2(:,2),cv2(:,3)]);
in2 = dmat2 < (kerd2/2) & dmat2 ~= 0; %Creates the same distance matrix as above for size 2 particles.
in2 = triu(in2);
[rc2, cc2] = find(in2 == 1);

for l = 1:length(rc1)
    I11 = cv1((rc1(l)),1);
    I12 = cv1((cc1(l)),1); %Pulls the intensity from every potential double peak location
    if I11 > I12
        cv1((cc1(l)),1) = 0; %Assigns the lower intensity to be zero
    elseif I12 > I11
        cv1((rc1(l)),1) = 0;
    else
        cv1((cc1(l)),1) = 0;
        cv1((rc1(l)),2) = mean(cv1((cc1(l)),2),cv1((rc1(l)),2));
        cv1((rc1(l)),3) = mean(cv1((cc1(l)),3),cv1((rc1(l)),3));
    end %If both intensities are equal, the location between the peaks is used as a best guess.
end

for l = 1:length(rc2)
    I21 = cv2((rc2(l)),1);
    I22 = cv2((cc2(l)),1); %Repeats the above process for size 2 particles.
    if I21 > I22
        cv2((cc2(l)),1) = 0;
    elseif I22 > I21
        cv2((rc2(l)),1) = 0;
    else
        cv2((cc2(l)),1) = 0;
        cv2((rc2(l)),2) = mean(cv2((cc2(l)),2),cv2((rc2(l)),2));
        cv2((rc2(l)),3) = mean(cv2((cc2(l)),3),cv2((rc2(l)),3));
    end
end

cv1 = sortrows(cv1, 1); %Sorts the rows based on intensity (lowest to highest)
cv1 = flipud(cv1); %Flips the vector so the most intense value is at the top
cv1 = cv1(1:np1, :); %Takes the first np1 particles (since they have the brightest peaks and are thus the locations)
cv1 = cv1(:,2:3) - 2*edg - kerd1/2 - .5 - 2; %Adjusts the center locations from where they are reported to where they actually are.

cv2 = sortrows(cv2, 1);
cv2 = flipud(cv2);
cv2 = cv2(1:np2, :);
cv2 = cv2(:,2:3) - 2*edg - kerd2/2 - .5 - 2; % ^

end