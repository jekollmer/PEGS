%This function creates a synthetic photoelastic response image 
%last change on 2016/09/27 by Jonathan Kollmer


% %test values
% z=2; %Number of contacts this particle has
% f = [0.6 0.6]; %Absolute forces on this particle
% alpha = [0 0];  %Alpha contact angles on this particle
% beta = [0 -pi];  %Beta contact angles on this particle
% fsigma = 100;  %Photoelastic stress coefficient of this particle
% rm = 0.00816516; %Particle radius in meters
% px = 100; %Return image size in pixels



function img = joForceImg (z, f, alpha, beta, fsigma, rm, px, verbose)

    %make sure the forces are balanced
    [alpha,f] = forceBalance(f,alpha,beta);

    %create an empty placeholder image for our result/return value
    img = zeros(px);

    %Create a scale that maps the diameter of the particle onto the image size
    xx = linspace(-rm, rm, px); 
    parfor (x=1:px) %loop throgh image width
        xRow=zeros(px,1); %placeholder for the current row beeing prcessed, parallel for makes it necessary to split up the result data and later combine it
        for y=1:px  %loop throgh image height
            if ((xx(x)^2+xx(y)^2)<=rm^2) %check if we are actually inside
                xRow(y) = StressEngine(xx(x), xx(y), z, f, alpha, beta, fsigma, rm);  %call the StressEngine to compute the photoelastic response at each pixel
            end
        end
        img(x,:)=xRow; %consolidate processed row into the output image, necessary data shuffling to use parallel for
    end

    if(verbose)
        %plot the synthetic image each time it is generated (i.e. you will
        %see the fit converge on screen)
        subplot(1,2,2)
        imshow(img);
        drawnow
    end
end
