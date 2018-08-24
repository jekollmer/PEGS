function [ img1, img2 ] = fftfilt( img )
% fftfilt uses the fast fourier transform the apply a high pass filter to
% remove the large amount of mean values. This evens out our images
% lighting.
%
%   Input:  img - The image that is to be filtered.
%
%   Output: img - The filtered image

[r,ci,~] = size(img); %Calculates the size of the image

edg = round((r+ci)/8); %%%%% This is the border size that will be smoothed %%%%%

cm = ((1 + tanh((abs(0:(2*edg)) - edg/2)/(2/3*edg)))/2)';
rm = cm';
cm = repmat(cm, [1, (ci+(4*edg))]);
rm = repmat(rm, [(r+(4*edg)), 1]); %Uses hyperbolic tangest to make a nice curve to reduce how sharp the edges are.

imgo = img; %This just allows the usage a saved original image since we now output red and green components.

for l = 1:2
    
    img = double(imgo(:,:,l)); %Converts the image to a double
    
    img = padarray(img,[(2*edg),(2*edg)],mean(img(:)));
    
    img(1:round((2*edg+1)),:) = img(1:round((2*edg+1)),:) .* cm;
    cm = flipud(cm);
    img((end-round(2*edg)):(end),:) = img((end-round(2*edg)):(end),:) .* cm;
    
    img(:,1:round(edg*2+1)) = img(:,1:round(edg*2+1)) .* rm;
    rm = fliplr(rm);
    img(:,(end-round(2*edg)):(end)) = img(:,(end-round(2*edg)):(end)) .* rm; %Applies the edge smoothing to our image
    
    pc = 3; %%%%% This is the padding constant which determines how large pad will be ((2*pc+1)*imagerows will be the total amount of rows for isntance) %%%%%
    
    img = padarray(img,[pc*(r+2*round(edg)),pc*(ci+2*round(edg))],0); %Pads our array with zeros so our end image 7 times as large as the original (size is now [7r,7c]).
    
    fftim = fft2(img); %Performs the fft of our image
    fftim = fftshift(fftim); %Shifts our fft'd image
    [rf,cf] = size(fftim);
    
    c = round((rf+cf)/500); %%%%% This is the radius of the circle used for filtering %%%%%
    w = c/10; %Sets a weight and c for the high pass filter
    
    fftims = fftim((round(end/2)-c):(round(end/2)+c),(round(end/2)-c):(round(end/2)+c)); %The center section of the fft image.
    d = sqrt((-c:c).^2+(-c:c)'.^2); %Takes a section of the fft image and finds a distance matrix and logical index.
    in = d < c; %This is the index of values within the circular range
    
    msk = zeros(size(in));
    msk(in) = (1 + tanh((abs(d(in)) - c)/w))/2; %Creates a mask for the filter
    
    msk(msk == min(msk(:))) = 0; %Sets the minimum to the mask to be zero as opposed to something close to zero
    fftims = fftims .* msk; %Applies the filter
    
    fftim((round(end/2)-c):(round(end/2)+c),(round(end/2)-c):(round(end/2)+c)) = fftims;
    fftim = fftshift(fftim); %Applies the filtered image to the shifted image and unshifts it
    img = ifft2(fftim); %Returns the image to real space
    
    img = (real(img((pc*(r+2*round(edg))+2*edg):((pc+1)*(r+2*round(edg))-1), ...
        (pc*(ci+2*round(edg))+2*edg):((pc+1)*(ci+2*round(edg))-1))));
    
    if l == 1
        img1 = img;
    else
        img2 = img; %Assigns red and green outputs properly.
    end
end

end

