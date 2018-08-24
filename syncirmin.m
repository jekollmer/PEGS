function [ hol ] = syncirmin( dat, mapv )
% syncirmin exists to be used in the syncirdif function. This function
% takes in a small section of data and a portion of the map variable (used
% in optimization) and calculates the difference between the data section
% and the circle made from the map variable. This then gets repeated using
% lsqnonlin until it is optimal.
%
%   Input:  dat - Refered to as datt1 or datt2 in syncirdif. A small
%           section of the total data (image) cut around a center location
%           so the synthetic circle is accurate to that pixel.
%
%           mapv - Map variable. This is essentially the x0 guess and is
%           ordered to be [row,column,radius,weight]. This value is optimized.
%
%   Output: hol - The synthetic circle subtracted from the data. This value
%           is minimized and has no real need to be saved as output.

l1 = (1):(length(dat)); % Kernel size

rect = zeros(size(dat)); %Creates a fake image

d = sqrt((-l1 + mapv(2)).^2 + (-l1' + mapv(1)).^2); %Calculates distance array
in = d < mapv(3)/2; %Creates index based on values within the radius

rect(in) = (1 - tanh((abs(d(in)) - mapv(3)/2)/mapv(4)))/2; %Assigns a smooth circle

hol = abs((dat/(max(dat(:)))) - rect); %Subtracts the fake circle from the data. This value is minimized.

end

