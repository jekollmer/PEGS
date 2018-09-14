function [ conv ] = kernfind( img, ker )
%kernfind takes an image and a kernel and converts them into Fourier space
%in order to more efficient convelute. The two are multiplied in fourier
%space and brought back using fft2 and ifft2 respectively to produce an
%image that gives the highest "matches" found.
%
%   Input:  img - The image that is to be conveluted. Input is expected to
%                 be a 2D array, however, 3D arrays are acceptable but will
%                 be converted using rgb2gray
%           
%           ker - The kernel that is to be used to search the image. As
%                 above, a 2D array is expected but 3D arrays can be used.
%
%   Output: conv - The convolution of the kernel and image.

if ndims(img) > 2
    fprintf('Image found to be 3D array. Converting to 2D.\n')
    img = rgb2gray(img);
end

if ndims(ker) > 2
    fprintf('Kernel found to be 3D array. Converting to 2D.\n')
    ker = rgb2gray(ker);
end

%Above are error checks in case the function is used without the rest of
%the wrapper.

pd = 600; %%%%% Padding distance because fft2 hates edges %%%%%
img = img ./ max(img(:));
img = padarray(img,[pd,pd],1,'both'); %Pads image

[r,c] = size(img); %Sizes image 

cm = ((1 + tanh((abs(0:pd) - pd/2)/(pd/3)))/2)'; %I hear what you're saying "Josh, you can't just tanh the shit out of all your data!" Well, watch me.
rm = cm';
cm = (repmat(cm, [1, c]));
rm = (repmat(rm, [r, 1]));

img(1:(1+pd),:) = img(1:(1+pd),:) .* cm;
img(:,1:(1+pd)) = img(:,1:(1+pd)) .* rm;
cm = flipud(cm);
rm = fliplr(rm);
img((end-(pd)):(end),:) = img((end-(pd)):(end),:) .* cm;
img(:,(end-(pd)):(end)) = img(:,(end-(pd)):(end)) .* rm;

[rk,ck] = size(ker); %Both find the matrix size for later use
kerm = zeros([r,c]); %Create a kernel matrix of zeros
kerm(1:rk,1:ck) = ker; %Assign the kernel matrix

imf = fft2(img); 
kerf = fft2(kerm); %Performs the Fourier Transform of the image and kernel matrix

convf = imf.*kerf; %Multiplies the transforms to convelute in real space

conv = ifft2(convf); %Returns the conveluted matrix to real space.
conv = conv(pd:(end-pd-1),pd:(end-pd-1));

end