function particle = disksolve(particle, maskradius, scaling, fileName, verbose, fitoptions)
    pres=particle;
    N = length(particle); %number of particles in this frame
    
    disp(['processing file ',fileName, 'containing ' ,num2str(N), 'particles']); %status indicator

    parfor (n=1:N)
        display(['fitting force(s) to particle ',num2str(n)]); %status indicator
        if (particle(n).z > 0 )
                    %bookkeeping
                    fsigma = particle(n).fsigma;
                    %fsigma = 80;
                    rm = particle(n).rm;

%                     %This is the Camera Image
%   %template = im2double(particle(n).forceImage);
                     template = particle(n).forceImage;
                     template = imresize(template,scaling);
%                     template = template-0.1;
%                     %template = (template -0.2); %fine tunes the image, this should happen in preprocessing!
%                     template = template.*(template > 0); %fine tunes the image, this should happen in preprocessing!
%                     template = template*3; %fine tunes the image, this should happen in preprocessing!

                    %size of the force image
                    px = size(template,1); 
                    
                    %plot the experimental image that is going to be fitted
                    %onto
                    if verbose
                        subplot(1,2,1)
                        imshow(template);
                    end
                    
                    %Create initial guesses for each contact force, based
                    %on the gradient squared in the contact region (also
                    %better guess lower than too high)
                    z = particle(n).z;
                    forces = zeros(z,1);
                    cg2s = sum(particle(n).contactG2s);

                    for i=1:z
                            forces(i)=2*particle(n).f*particle(n).contactG2s(i)/cg2s;%initial guess for force
                    end

                    %Initial guesses for the angles of attack of the forces
                    %for each contact
                    alphas = zeros(z,1);
                    beta = particle(n).betas;

                    %apply force balance to the initial guesses
                    [alphas,forces] = forceBalance(forces,alphas,beta);

                    %create a circular mask
                    cx=px/2;cy=px/2;ix=px;iy=px;r=maskradius*px;
                    [x,y]=meshgrid(-(cx-1):(ix-cx),-(cy-1):(iy-cy));
                    c_mask=((x.^2+y.^2)<=r^2);   

                    %This is the function we want to fit, i.e. a synthetic
                    %version of the force image with free parameters f (here called par(1:z)) and
                    %alpha (here called par(z+1:z*z) since we have to stuff
                    %them into one vector.
                    func = @(par) joForceImg(z, par(1:z),par(z+1:z+z), beta(1:z), fsigma, rm, px, verbose); %+par(2*z+1); %this is the function I want to fit (i.e. synthetic stres image), the fitting paramters are in vector par
                    %This is the error function we are actually fitting,
                    %that is, the distance between our fit function and the
                    %real particle image in terms of the sum of squares of the pixelwise differnce.
                    %Also a mask is applied to crop
                    %out only the circular part of the particle. 
                    err = @(par) real(sum(sum( ( c_mask.*(template-func(par)).^2) ))); %BUG: for some reason I sometimes get imaginary results, this should not happen

                    %Set up initial guesses
                    p0 = zeros(2*z, 1);
                    p0(1:z) = forces;
                    p0(z+1:2*z) = alphas;


                    %Do the fit, will also work with other solvers
                    %TODO: make a user defined option to select between
                    %different solvers
                    p=lsqnonlin(err,p0,[],[],fitoptions);

                    %get back the result from fitting
                    forces = p(1:z);
                    alphas = p(z+1:z+z);

                    %resudual
                    fitError = err(p);
                    
                    %generate an image with the fitted parameters
                    imgFit = joForceImg(z, forces, alphas, beta, fsigma, rm, px*(1/scaling), verbose);
                    
                    %since the image gnerator also forces force balance we
                    %have to explicitly do it once more to the values we
                    %are gonna save 
                    [alphas,forces] = forceBalance(forces,alphas,beta);
                    
                    %store the new information into a new particle vector 
                    pres(n).fitError = fitError;
                    pres(n).forces = forces;
                    pres(n).alphas = alphas;
                    pres(n).synthImg = imgFit;
        end
    end
    %save the result
    save([fileName(1:end-17),'solved.mat'],'pres');
end