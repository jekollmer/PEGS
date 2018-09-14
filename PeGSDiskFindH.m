function [ particle ] = PeGSDiskFindH( Rimg, Rlarge, SL, Rsmall, SS, pxPerMeter, fsigma)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%Disc detection for large particles happens here
[centers, radii, ~] = imfindcircles(Rimg,Rlarge,'ObjectPolarity','bright','Method','TwoStage','Sensitivity',SL);
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
[centers, radii, ~] = imfindcircles(Rimg,Rsmall,'ObjectPolarity','bright','Method','TwoStage','Sensitivity',SS);
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

for n=1:N
    particle(n).rm = particle(n).r * pxPerMeter;
    particle(n).fsigma = fsigma;
    if particle(n).r < ( mean(Rlarge)+ mean(Rsmall) ) /2
        particle(n).color = 'b';
    else
        particle(n).color = 'r';
    end
end

end

