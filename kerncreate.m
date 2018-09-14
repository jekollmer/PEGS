function [ kernel ] = kerncreate( opt, kdim, in1, varargin )
%kerncreate - Creates a kernel that is either a filled circle, a circle
%with an empty inner ring, or a circle with a cross or 'X' in the middle.
%Output is a matrix with the points of interest assigned to a value of 1
%and background values set to zero.
%
%   Input:  opt - Option for circle type 
%               0 - Filled circle
%               1 - Circle with inner empty ring
%               2 - Circle (ring) with an 'X' in middle.
%
%           kdim - Dimension of kernel matrix (i.e. side padding)
%
%           If opt == 0
%               in1 - Radius of circle
%               in2 - Not used, don't put anything in.
%
%           If opt == 1
%               in1 - Smaller (inner) Radius
%               varargin - Larger (outer) Radius
%
%           If opt == 2
%               kdim - kdim is no longer padding but the diameter of our
%               circle
%               in1 - The offput angle from the horizontal (counter
%                   clockwise)
%               varargin - Width of the x within the circle. Additionally,
%               this is the width of our ring.
%
%   Output: kernel - Requested kernel

if nargin == 4 && opt ~= 0
    fprintf('Error: Not enough inputs provided. \n')
    kernel = NaN;
    return;
elseif nargin == 4 && opt == 0 || nargin > 4
    fprintf('Error: Excessive input provided. \n')
    kernel = NaN;
    return;
end %Error checks (^)


kernel = zeros(kdim); %Assigns the kernel to initially be a zero square of zeros
c = (kdim+1)/2; %Assigns the center of our circle


if opt == 0
    l1 = 1:kdim; %Assigns our total size of our image (should be diameter + 1)
    d = sqrt((l1 - c).^2 + (l1' - c).^2); %Creates distance array
    in = d < in1; %Creates a logical index for distances within the radius
    kernel(in) = 1; %Makes the index locations 1.
    
elseif opt == 1
    in2 = varargin{1}(1); %Assigns our second radius to be the varargin input
    
    if in1 < in2 %Fixes radii if misordered
        tr = in1; % Assigns a temporary radius
        in1 = in2; % Makes radius 1 into radius 2
        in2 = tr; % Reassigns radius 2
    end
    
    l1 = 1:kdim; %Total kernel size
    d = sqrt((l1 - c).^2 + (l1' - c).^2); %Distance array
    in1 = d < in1; %Logical index for first radius
    in2 = d < in2; %Logical index for smaller radius
    kernel(in1) = 1; %Sets all values within the large radius to be 1
    kernel(in2) = 0; %Sets inner circle to be 0
    
% opt 2 and 3 aren't really used so I'm not going to bother making them better at the moment.    

elseif opt == 2
    
    w = varargin{1}(1); %Pulls the width from our variable input
    ang(1) = in1*180/pi; %Assigns our first angle in radians
    ang(2) = ang + pi/2; %Shifts our first angle pi/2 radians
    r2 = kdim/2; %Assigns our larger radius to be the half kdim
    r1 = r2 - w; %Subtracts the width of our 'X' to get the second radius
    
    for l1 = 1:2
        for l2 = 1:kdim
            for l3 = 1:kdim
                r = sqrt((l2-c)^2 + (l3-c)^2);
                if l2 == l3
                    thet = pi/2;
                else
                    thet = -atan2(l3-c,l2-c);
                end
                s = abs(r*sin(thet - ang(l1)));
                if s <= w/2.0
                    kernel(l3,l2) = 1;
                end;
                if (r > r2)
                    kernel(l3,l2) = 0;
                end;
                if ((r < r2) && (r > r1))
                    kernel(l3,l2) = 1;
                end;
            end
        end
    end
    
elseif opt==3
    
    for l1 = 1:kdim
        for l2 = 1:kdim
            if (l1-c)^2 + (l2-c)^2 >= .85*in1^2 && (l1-c)^2 + (l2-c)^2 >= 1.15*in1^2
                kernel(l1,l2) = .5; %Checks if point is at the edge
            end
            if (l1-c)^2 + (l2-c)^2 < in1^2 
                kernel(l1,l2) = 1; %Checks if the point is within the radius
            end
        end
    end
    
else
    fprintf('Error: Invalid option selection. \n')
    kernel=NaN;
end

kernel = kernel - mean(kernel(:)); %Subtracts the mean from the kernel so the background value is negative and the internal value is positive.

end

