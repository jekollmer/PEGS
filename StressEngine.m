%Computes the optical response of any given point in a loaded photoelastic disk
%Originally from James G. Puckett and Karen E. Daniels. "Equilibrating temperaturelike variables in jammed granular subsystems." Physical Review Letters 110 058001 (2013) 
%Ported to Matlab by Jonathan E. Kollmer, 2016-05-06

%result = photoelastic response (intensity)

%xxi = x-position at which the stress is computed (in meters)
%xxy = y-position at which the stress is computed (in meters)
%z = coordination number
%f list of forces (magintude) acting on the particle
%alpha list of angles of the forces, see diagram in J.Puckett, PhD-Thesis, North Carolina State University
%beta list of agles of the forces, see diagram in J.Puckett, PhD-Thesis, North Carolina State University
%fsigma photoelastic stress coefficient
%rm radius of the disk (in meters)



function result = StressEngine(xxi, xxj, z, f, alpha, beta, fsigma, rm)

    beta = -beta+pi/2;
    pioverfsigma = pi/(fsigma);
    oneoverpirm  = 1./(pi*rm);
    twooverpi    = 2./pi;
     
	sigmaxx=0;
    sigmayy=0;
    sigmaxy=0;
    %double aa,a,b,b2,x1,y1,x2,y2,ch0,ch1,chn,costh1,r10,r11,r1n,th1,s1, s2,sr,th; %temporary
    for k=1:z
		b   = beta(k)+pi/2; 		%rotate pi/2, to match input image     
        a   = alpha(k);
        if (a<0) 
           b2 = b+(pi+2*a);
        else
           b2 = b-(pi-2*a);
        end
        x1  = rm*sin(b);
        y1  = rm*cos(b);
        x2  = rm*sin(b2);
        y2  = rm*cos(b2);
        ch0 = x2-x1;                % chord x
        ch1 = y2-y1;                % chord y
        chn = sqrt(ch0^2+ch1^2);
        ch0 = ch0/chn;              % normalize chord x
        ch1 = ch1/chn;              % normalize chord y
        r10 = xxi - x1;             % r vector x coord
        r11 =-xxj - y1; 			% r vector y coord
        r1n = sqrt(r10^2+r11^2);
        costh1 = (r10*ch0+r11*ch1)/r1n;
        if (r11*ch0>r10*ch1) % important!
            signth =  1;
        else
            signth = -1;
        end 
        th1    = signth*acos(costh1);     %faster than cos(asin(stuff));
        s1 	= -(f(k)*twooverpi)/(r1n)*costh1;
        th 	= th1-beta(k)-alpha(k);              %rotate coordinates
        sigmaxx = sigmaxx + s1*((sin(th))^2);
        sigmayy = sigmayy + s1*((cos(th))^2);
        sigmaxy = sigmaxy + 0.5*s1*(sin(2*th));
	
   end %end of k loop over beta, alpha
   aa     = real(sqrt((sigmaxx-sigmayy)^2+4*(sigmaxy)^2));
   result = (sin(pioverfsigma*aa))^2;  %wrap    = sin(2*pi*t/fsigma*a).^2;
   if (isnan(result)) %for some reason result sometimes gets to be NAN, not sure why this happens
       result = 0; %temproary fix is to set it zero then
   end 
   if (result<0) %for some reason, not sure why this happens
       result = 0; %temproary fix is to set it zero then
   end 
end
