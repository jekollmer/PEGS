%%forceBalance.m
%%This imposes force balance on a given set of forces and their alpha and
%%beta angles as defined by fig. 4.11 in James Puckett's thesis.
%%The code for z==2 is adapted from peDiskSolve
%%and the code for z>2 is my reimplementation of eqations 4.22 in James' theis.
%%Last Edit on 2016/09/22 by Jonathan Kollmer


function [alpha,f] = forceBalance(force,alpha,beta)
    beta = -beta+pi/2;
    z = length(beta);
     if (z<2)
          f=force;
     elseif (z==2) 
          dbeta = (beta(1)-beta(2))/2;
          f(1)     = force(1);
          alpha(1) = acos(sin(dbeta));
          if (alpha(1)>pi/2) 
               alpha(1)=acos(sin(-dbeta));
          end
          f(2)     = f(1);
          alpha(2) = - alpha(1);
     else 

         for k = 1:z 
             sum1 = 0;
             sum2 = 0;

             for i = 1:z
                 if(i~=k)
                 sum1 = sum1 + force(i)*sin(alpha(i)+beta(i)-beta(k));
                 sum2 = sum2 + force(i)*cos(alpha(i)+beta(i)-beta(k));
                 end
             end
             f(k) = sqrt(sum1^2+sum2^2); 
         end
     
         for k = 1:z 
             sum3 = 0;
             for i = 1:z
                 if(i~=k)
                 sum3 = sum3 + f(i)*sin(alpha(i));
                 end
             end;
             a(k) = asin(-sum3/f(k));    
         end
     alpha=a;
     end

end