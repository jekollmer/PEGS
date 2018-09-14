function [ syncir ] = synthcirccreate( dat, mapv, np1, np2 )
% Creates a snythetic image of our circles using the map of each kernel and
% then creating circles at center locations.
%
%	Input:  dat - Image being analyzed.
%
%           mapv - Map variable that is created from the optimization
%
%           np1 - Number of particles of size 1
%
%           np2 - Number of particles of size 2
%
%   Output: syncir - Synthetic circles.

[r,c] = size(dat);
pad = 240; %%%%% Assigns a pad size. This should usually be twice your largest diameter. %%%%%
rect=zeros(r+pad,c+pad); %Creates a blank image

for l = 1:np1 %Runs for every size 1 particle. 
    l1 = 0:(ceil(mapv(2*np1+l))); %Creates a vector the size of the diameter

    d = sqrt((l1 - mapv(2*np1+l)/2).^2 + (l1' - mapv(2*np1+l)/2).^2); %Distance from center point of a circle
    in = d < mapv(2*np1+l)/2; %Logical index for distances inside radius
    rects = zeros(size(in));  %Section of our blank image
    rects(in) = (1 - tanh((abs(d(in)) - mapv(2*np1+l)/2)/mapv(3*np1+l)))/2; %Assigns the section of our blank image to be a circle with fancy lighting (tanh function)
    rect(round(l1 + mapv(l) - mapv(2*np1+l)/2) + pad/2,round(l1 + mapv(np1+l) - mapv(2*np1+l)/2) + pad/2) ... 
        = rect(round(l1 + mapv(l) - mapv(2*np1+l)/2) + pad/2,round(l1 + mapv(np1+l) - mapv(2*np1+l)/2) + pad/2) + rects;
    %Reassigns the circle in our original image by adding it.
end

for l = 1:np2
    l2 = 0:ceil(mapv(4*np1+2*np2+l)); 
 
    d = sqrt((l2 - mapv(4*np1 + 2*np2 + l)/2).^2 + (l2' - mapv(4*np1 + 2*np2 + l)/2).^2); 
    in = d < mapv(4*np1+2*np2+l)/2; 
    rects = zeros(size(in)); 
    rects(in) = (1 - tanh((abs(d(in)) - mapv(4*np1+2*np2+l)/2)/mapv(4*np1+3*np2+l)))/2; 
    rect(round(l2 + mapv(4*np1 + l) - mapv(4*np1+2*np2+l)/2) + pad/2,round(l2 + mapv(4*np1 + np2 + l) - mapv(4*np1+2*np2+l)/2) + pad/2) ... 
        = rect(round(l2 + mapv(4*np1 + l) - mapv(4*np1+2*np2+l)/2) + pad/2,round(l2 + mapv(4*np1 + np2 + l) - mapv(4*np1+2*np2+l)/2) + pad/2) + rects;
end

syncir = rect(pad/2:(end-1-pad/2), pad/2:(end-1-pad/2)); %Cuts off the extra padding of our image and assigns it to our output

end

