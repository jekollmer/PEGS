%This function creates a synthetic photoelastic response image 
%last change on 2016/08/05 by Jonathan Kollmer


% %test values
% z=2; %Number of contacts this particle has
% f = [0.6 0.6]; %Absolute forces on this particle
% alpha = [0 0];  %Alpha contact angles on this particle
% beta = [0 -pi];  %Beta contact angles on this particle
% fsigma = 100;  %Photoelastic stress coefficient of this particle
% rm = 0.00816516; %Particle radius in meters
% px = 100; %Return image size in pixels
% aa = 1; %anti-aliasing factor


function img = joForceFunc (z, f, alpha, beta, fsigma, rm, px, aa)

    %in case we use anti-aliasing, blow up image size here (will be
    %downsampled again at the end of this function)
    px = px*aa; 
    
    %create an empty placeholder image for our result/return value
    img = zeros(px);

    %Create a scale that maps the diameter of the particle onto the image size
    xx = linspace(-rm, rm, px); 

    parfor x=1:px %loop throgh image width
        xRow=zeros(px,1); %placeholder for the current row beeing prcessed, parallel for makes it necessary to split up the result data and later combine it
        for y=1:px  %loop throgh image height
            %r=px/2;
            %if ((px^2+px^2)<=r^2) %check if we are actually inside
            %particle, there is a bug in this check
                xRow(y) = StressEngine(xx(x), xx(y), z, f, alpha, beta, fsigma, rm);  %call the StressEngine to compute the photoelastic response at each pixel
            %end
        end
        img(x,:)=xRow; %consolidate processed row into the output image, necessary data shuffling to use parallel for
    end
    
    %create a circular mask
    cx=px/2;cy=px/2;ix=px;iy=px;r=px/2;
    [x,y]=meshgrid(-(cx-1):(ix-cx),-(cy-1):(iy-cy));
    c_mask=((x.^2+y.^2)<=r^2);
    
    %apply the circular mask to the result
    img = img .* c_mask;
    
    %in case we use anti-aliasing, do the downsamlpling here
    if (aa ~= 1) 
        img = imresize(img,1/aa,'bilinear'); %2x2 average downsampling
    end

end