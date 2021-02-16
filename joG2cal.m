close all
clear all
%G2 Calibration Tool


% %test values
z=2; %Number of contacts this particle has
%f = [0.6 0.6]; %Absolute forces on this particle
alpha = [0 0];  %Alpha contact angles on this particle
beta = [0 -pi];  %Beta contact angles on this particle
fsigma = 100;  %Photoelastic stress coefficient of this particle
rm = 0.00816516; %Particle radius in meters
px = 100; %Return image size in pixels
aa = 1; %anti-aliasing factor

data=[];
n=1;

for n=1:100
    
    f1= 0.01*n ;
    f=[f1 f1];
    data(n).f = f1;
    
    %compute photoelasstic response
    particleImg = joForceFunc (z, f, alpha, beta, fsigma, rm, px, aa);
    
    %Compute G^2 
    [gx,gy] = gradient(double(particleImg));
    g2 = gx.^2 + gy.^2;
    data(n).g2 = sum(sum(g2));
    
end

%Fit a straigt Line trough the data
P = polyfit([data.f],[data.g2],1);
yfit = P(1)*[data.f]+P(2);
    
   
plot([data.f],[data.g2],'ko'); hold on;
plot([data.f],yfit,'r-','LineWidth',2);
text(0.2,200,sprintf('%.3f * f + %0.3f',P(1),P(2))),'FontSize',16';

