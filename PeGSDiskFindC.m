function [particle] = PeGSDiskFindC(img, Rsmall, Nsmall, Rlarge, Nlarge)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

kerd1=Rsmall*2;
kerd2=Rlarge*2;

[imrf, imgf] = fftfilt(img); %Filters the image

for l = 1:5
    imrf = im2double(imrf);
    imrf(imrf > 3/4 * max(imrf)) = 2/3*max(imrf(:));
end
imgf = im2double(imgf);
im = imrf + imgf/(1/5*max(imgf(:)));
im = imadjust((im),stretchlim(im));
im = im + 1/4*imrf;
im(im < -10) = im(im < -10)/min(im(im < -10));

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

dmat1 = pdist2([cv1(:,1),cv1(:,2)],[cv1(:,1),cv1(:,2)]);
in1 = dmat1 < kerd1/2 & dmat1 ~=0;

dmat2 = pdist2([cv2(:,1),cv2(:,2)],[cv2(:,1),cv2(:,2)]);
in2 = dmat2 < kerd2/2 & dmat2 ~=0;

if sum(in1(:)) > Nsmall/10 || sum(in2(:)) > Nlarge/10
    im = -im;
    imn1=kernfind(im,ker1);
    imn1 = imn1./max(imn1(:));
    imn2=kernfind(im,ker2);
    imn2 = imn2./max(imn2(:));
    imn1(1:kerd1,:) = 0;
    imn1(:,1:kerd1) = 0;
    imn2(1:kerd1,:) = 0;
    imn2(:,1:kerd1) = 0;
    [cv1,cv2] = max2D(imn1, imn2, kerd1, kerd2, Nsmall, Nlarge);
end

cv1 = [cv1,ones([length(cv1),1])*kerd1,ones(length(cv1),1)*1]; %Creates vectors such that they are [row,column,diameter,weight]
cv2 = [cv2,ones([length(cv2),1])*kerd2,ones(length(cv2),1)*1];
mt1 = reshape(cv1,[1,4*length(cv1)]); %Reshapes vectors
mt2 = reshape(cv2,[1,4*length(cv2)]);
x0 = [mt1,mt2]; %Creates map

xv = syncirdif(double(img(:,:,1)), x0, Nsmall, Nlarge, rkerd1, rkerd2);

cir1 = [xv(1:Nsmall)',xv((Nsmall+1):(2*Nsmall))',(xv((2*Nsmall+1):(3*Nsmall))/2)',(xv((3*Nsmall+1):(4*Nsmall)))']; %[r,c,r,w]
cir2 = [xv((4*Nsmall+1):(4*Nsmall+Nlarge))',(xv((4*Nsmall+Nlarge+1):(4*Nsmall+2*Nlarge)))',(xv((2*Nlarge+4*Nsmall+1):(3*Nlarge+4*Nsmall))/2)',(xv((3*Nlarge+4*Nsmall+1):(4*Nsmall+4*Nlarge)))'];
circ = [cir1;cir2];
circ(:,2:3) = circ(:,2:3) - 2;

circ = sortrows(circ,3);
circs = flipud(circ);
circ = circs;
N = length(circ);

particle(1:N) = struct('id',0,'x',0,'y',0,'r',0,'rm',0,'color','','fsigma',0,'z',0,'f',0,'g2',0,'forces',[],'betas',[],'alphas',[],'neighbours',[],'contactG2s',[],'forceImage',[]);
for l = 1:N
    particle(l).x = circ(l,2);
    particle(l).y = circ(l,1);
    particle(l).r = circ(l,3);
    particle(l).id = l;
    if l <= Nsmall
        particle(l).color = 'b';
    else
        particle(l).color = 'r';
    end
end

end

