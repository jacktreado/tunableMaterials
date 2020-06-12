%% Simulation to do 3 strain rate sims

clear;
close all;
clc;

%% Initialize system

N               = 48;       % number of particles
phi0            = 0.75;     % initial packing fraction
dphi0            = 0.001;    % initial delta phi
dt0             = 0.01;     % time step magnitude
sr              = 1.4;      % size ratio
seed            = 1;        % random number seed
T               = 1;        % initial temperature

Ftol            = 1e-12;    % force tolerance
Ktol            = 1e-25;    % kinetic tolerance
Ptol            = 1e-8;     % pressure tolerance

strainAmp       = 1e-1;     % strain amplitude
strainFreq      = 5e-2;     % strain angular frequency
NCYCLES         = 40;       % number of strain cycles
bdamp           = 1.0;      % damping parameter
makeAVideo      = 0;        % whether or not a video should be made

dphi            = 0.01;     % particle size increase

% strain cycle period
tCycle          = (2.0*pi)/strainFreq; 

% number of iterations to skip plotting
plotskip        = 1e4;

% seed random number generator
rng(seed);

% UNCHANGING FIRE VARIABLES
alpha0          = 0.25;
finc            = 1.1;
fdec            = 0.5;
falpha          = 0.99;

NMIN            = 50;
NNEGMAX         = 2000;
NDELAY          = 1000;


%% Initialize system 

% number of small and large particles
Nsmall = round(0.5*N);
Nlarge = N - Nsmall;

% particle radii
r = ones(N,1);
r(Nsmall+1:N) = sr*ones(Nlarge,1);

% box length
L = sqrt(pi*sum(r.^2)/phi0)*ones(2,1);
Lx = L(1);
Ly = L(2);

% positions
x = Lx.*rand(N,1);
y = Ly.*rand(N,1);

% initialize velocities
vx = randn(N,1);
vy = randn(N,1);

% subtract off center of mass
vx = vx - mean(vx);
vy = vy - mean(vy);

% current T
T0 = 0.5*sum(vx.^2 + vy.^2);

% scale to be T
vx = sqrt(T/T0)*vx;
vy = sqrt(T/T0)*vy;

% forces
Fx      = zeros(N,1);
Fy      = zeros(N,1);

% virial stress tensor
vstress = zeros(4,1);

% contacts
cij = zeros(N);

% energies
K = 0;
U = 0;

% time step
dt = dt0;

%% Loop until jamming

% is jammed binary variable
isjammed = 0;

% iterate check
it = 0;
itmax = 1e6;

% FIRE iterate check
fireit = 0;

% initialize packing fraction
phi = phi0;

% binary root search variables
phiH = -1;
phiL = -1;

% rescale velocity parameter
Kr = 1e-6;

% save number of iterations required
NITSMAX = itmax;
nits = zeros(NITSMAX,3);

% loop over state
while (isjammed == 0 && it < itmax)
    % update iterator
    it = it + 1;
    
    % minimize energy using FIRE 2.0 with backstepping
    
    % FIRE VARIABLES
    dtmax       = 5*dt0;
    dtmin       = 0.01*dt0;
    
    npPos       = 0;
    npNeg       = 0;
    npPMin      = 0;
    
    alpha       = alpha0;
    
    fireit      = 0;
    fcheck      = 10*Ftol;
    kcheck      = 10*Ktol;
    
    
    while ((fcheck > Ftol || kcheck > Ktol || npPMin < NMIN) && fireit < itmax)
        % update iterate
        fireit = fireit + 1;
        
        % Velocity-verlet implementation in FIRE
        vx = vx + 0.5*dt*Fx;
        vy = vy + 0.5*dt*Fy;
        
        % Step 1. calculate P, fnorm, vnorm
        vnorm = sqrt(sum(vx.*vx) + sum(vy.*vy));
        fnorm = sqrt(sum(Fx.*Fx) + sum(Fy.*Fy));
        P = sum(Fx.*vx) + sum(Fy.*vy);

        % plot FIRE information 
        if mod(fireit,plotskip) == 0 || fireit == 1
            fprintf('\nOn FIRE step %d\n',fireit);
            fprintf('\t ** fcheck = %0.5g\n',fcheck);
            fprintf('\t ** kcheck = %0.5g\n',kcheck);
            fprintf('\t ** dt = %0.5g\n',dt);
            fprintf('\t ** P = %0.5g\n',P);
            fprintf('\t ** alpha = %0.5g\n',alpha);
            fprintf('\t ** vnorm = %0.5g\n',vnorm);
            fprintf('\t ** fnorm = %0.5g\n\n',fnorm);
        end

        % Step 2. adjust simulation based on net motion of system
        if (P > 0)
            % increase positive counter
            npPos = npPos + 1;

            % reset negative counter
            npNeg = 0;

            % update alphat for the next iteration
            alphat = alpha;

            % alter simulation if enough positive steps have been taken
            if (npPos > NMIN)
                % change time step
                if (dt*finc < dtmax)
                    dt = dt*finc;
                end
            end

            % decrease alpha
            alpha = alpha*falpha;
        else
            % reset positive counter
            npPos = 0;

            % increase negative counter
            npNeg = npNeg + 1;

            % check for stuck simulation
            if (npNeg > NNEGMAX)
                fprintf('Simulation negative for too long, ending program here.\n')
                error('FIRE did not converge');
            end

            % decrease time step if past initial delay
            if (it > NMIN)
                % decrease time step
                if (dt*fdec > dtmin)
                    dt = dt*fdec;
                end

                % change alpha
                alpha = alpha0;
            end

            % take a half step backwards
            x = x - 0.5*dt*vx;
            y = y - 0.5*dt*vy;

            % reset velocities to 0
            vx = zeros(N,1);
            vy = zeros(N,1);
        end

        % update velocities if forces are acting
        if fnorm > 0
            vx = (1 - alpha).*vx + alpha.*(Fx./fnorm)*vnorm;
            vy = (1 - alpha).*vy + alpha.*(Fy./fnorm)*vnorm;
        end
        
        % do first verlet update for vertices (assume unit mass)
        x = x + dt*vx;
        y = y + dt*vy;
        
        % do periodic boundary conditions
        x = mod(x,Lx);
        y = mod(y,Ly);

        % -- Compute new forces
        Fx = zeros(N,1);
        Fy = zeros(N,1);
        U = 0.0;
        vstress = zeros(4,1);
        cij = zeros(N);
        for ii = 1:N
            % particle i information
            xi = x(ii);
            yi = y(ii);
            ri = r(ii);
            
            % loop over particles j = [i+1,N]
            for jj = ii+1:N
                % sigma ij
                sij = ri + r(jj);
                
                % x distance
                dx = x(jj) - xi;
                dx = dx - Lx*round(dx/Lx);
                if dx < sij
                    % y distance
                    dy = y(jj) - yi;
                    dy = dy - Ly*round(dy/Ly);
                    
                    % calculate full distance
                    dist = sqrt(dx*dx + dy*dy);
                    
                    % check overlap
                    if dist < sij
                        % unit vector
                        ux = dx/dist;
                        uy = dy/dist;
                        
                        % force magnitude
                        ftmp = (1 - dist/sij)/sij;
                        
                        % force in given direction
                        fx = ftmp*ux;
                        fy = ftmp*uy;
                        
                        % update forces
                        Fx(ii) = Fx(ii) - fx;
                        Fy(ii) = Fy(ii) - fy;
                        
                        Fx(jj) = Fx(jj) + fx;
                        Fy(jj) = Fy(jj) + fy;
                        
                        % update virial stress
                        vstress(1) = vstress(1) + fx*dx/(Lx*Ly);    % sigmaXX
                        vstress(2) = vstress(2) + fx*dy/(Lx*Ly);    % sigmaXY
                        vstress(3) = vstress(3) + fy*dx/(Lx*Ly);    % sigmaYX
                        vstress(4) = vstress(4) + fy*dy/(Lx*Ly);    % sigmaYY
                        
                        % update contact network
                        cij(ii,jj) = 1;
                        cij(jj,ii) = 1;
                        
                        % update potential energy
                        U = U + 0.5*(1 - (dist/sij))^2;
                    end
                end
            end
        end
        
        % do second verlet update for vertices
        vx = vx + dt*0.5*Fx;
        vy = vy + dt*0.5*Fy;
        
        % update Kcheck and Fcheck
        fcheck = sqrt(sum(Fx.*Fx) + sum(Fy.*Fy))/N;
        kcheck = 0.5*sum(vx.*vx + vy.*vy)/N;
        
        if fcheck < Ftol
            npPMin = npPMin + 1;
            if npPMin >= NMIN && kcheck < Ktol
                fprintf('FIRE has converged in fireit = %d steps...\n',fireit);
                fprintf('fcheck = %0.3g\n',fcheck);
                fprintf('kcheck = %0.3g\n',kcheck);
                
            end
        else
            npPMin = 0;
        end
    end
    if fireit == itmax
        fprintf('FIRE protocol did not converge in less than itmax = %d iterations on it = %d, ending.\n',itmax,it);
        error('FIRE protocol did not converge.');
    end
    
    % save number of fire iteration
    nits(it,1) = phi;
    nits(it,2) = fireit; 
    nits(it,3) = phiH;
    
    % update virial pressure
    pcheck = 2.0*(vstress(1) + vstress(4))/(N*Lx*Ly);
    
    % print jamming information
    fprintf('\n\n** JAMMING DATA:\n');
    fprintf('it = %d\n',it);
    fprintf('pcheck = %0.3g\n',pcheck);
    fprintf('phi = %0.3f\n',phi);
    
    % get contacts
    ztmp = sum(cij);
    Nc = 0.5*sum(ztmp);
    
    % decide to change phi
    undercompressed = (pcheck < 2.0*Ptol && phiH < 0) || (pcheck < Ptol && phiH > 0);
    overcompressed = (pcheck > 2.0*Ptol && Nc > 0);
    jammed = (pcheck > Ptol & pcheck < 2.0*Ptol & phiH > 0);
    
    if phiH < 0
        if undercompressed
            % update packing fraction
            phiNew = phi + dphi0;
        elseif overcompressed
            % current = upper bound packing fraction
            phiH = phi;
            
            % save this packing fraction
            phiH0 = phi;
            
            % old = old upper bound packing fraction
            phiL = phiH - dphi0;

            % save upper bound positions
            xH = x;
            yH = y;

            % reset velocities
            vx = zeros(N,1);
            vy = zeros(N,1);
            
            % compute new packing fraction
            phiNew = phi - 0.5*dphi0;
            
            % print to console
            fprintf('\n\n ** OVERCOMPRESSED for the first time, phiNew = %0.8f\n',phiNew);
        end
    else
        if undercompressed
            % current = new lower bound packing fraction
            phiL = phi;
            
            % reset positions to initial positions
            x = xH;
            y = yH;
            
            % reset velocities
            vx = zeros(N,1);
            vy = zeros(N,1);
            
            % compute new packing fraction
            phiNew = 0.5*(phiH + phiL);
            
            % print to console
            fprintf('\n\n ** now undercompressed, phiNew = %0.8f\n',phiNew);
        elseif overcompressed
            % current = new upper bound packing fraction
            phiH = phi;
            
            % reset positions to initial positions
            x = xH;
            y = yH;
            
            % reset velocities
            vx = zeros(N,1);
            vy = zeros(N,1);
            
            % compute new packing fraction
            phiNew = 0.5*(phiH + phiL);
            
            % print to console
            fprintf('\n\n ** now overcompressed, phiNew = %0.8f\n',phiNew);
        elseif jammed
            fprintf('\t** At it = %d, jamming found!\n',it);
            fprintf('\t** phiJ = %0.3f\n',phi);
            fprintf('\t** pressure per particle = %0.3g\n',pcheck);
            fprintf('\t** Nc = %d\n',Nc);
            break;
        end
    end
    
    % compute new packing fraction
    phiOld = phi;
    rscale = sqrt(phiNew/phiOld);
    
    % update particle size based on packing fraction
    r = rscale.*r;
    
    % update this phi
    phi = pi*sum(r.^2)/(Lx*Ly);
end
if it == itmax
    fprintf('Jamming protocol did not converge in less than itmax = %d iterations, ending.\n',itmax);
    error('Jamming protocol did not converge.');
end

% remove extra iterations
nits(it+1:end,:) = [];
fprintf('\n\n** Number of iterations when phiH > 0: %d\n',sum(nits((nits(:,3)>0),2)));

% save initial jammed state coordinates
x1 = x;
y1 = y;
r1 = r;
vx1 = vx;
vy1 = vy;
U1 = U;
Fx1 = Fx;
Fy1 = Fy;
vstress1 = vstress;
cij1 = cij;

% initial packing fraction
phi = sum(pi*r.^2)/(Lx*Ly);
phi1 = phi;

%% Compress system by dphi, relax to new state

% increase particle sizes
rs  = sqrt((phi + dphi)/phi);
r   = rs.*r;
r2  = r;

% FIRE VARIABLES
dtmax       = 5*dt0;
dtmin       = 0.01*dt0;

npPos       = 0;
npNeg       = 0;
npPMin      = 0;

alpha       = alpha0;

fireit      = 0;
fcheck      = 10*Ftol;
kcheck      = 10*Ktol;

% relax new state using FIRE
while ((fcheck > Ftol || kcheck > Ktol || npPMin < NMIN) && fireit < itmax)
    % update iterate
    fireit = fireit + 1;

    % Velocity-verlet implementation in FIRE
    vx = vx + 0.5*dt*Fx;
    vy = vy + 0.5*dt*Fy;

    % Step 1. calculate P, fnorm, vnorm
    vnorm = sqrt(sum(vx.*vx) + sum(vy.*vy));
    fnorm = sqrt(sum(Fx.*Fx) + sum(Fy.*Fy));
    P = sum(Fx.*vx) + sum(Fy.*vy);

    % plot FIRE information 
    if mod(fireit,plotskip) == 0 || fireit == 1
        fprintf('\nOn FIRE step %d\n',fireit);
        fprintf('\t ** fcheck = %0.5g\n',fcheck);
        fprintf('\t ** kcheck = %0.5g\n',kcheck);
        fprintf('\t ** dt = %0.5g\n',dt);
        fprintf('\t ** P = %0.5g\n',P);
        fprintf('\t ** alpha = %0.5g\n',alpha);
        fprintf('\t ** vnorm = %0.5g\n',vnorm);
        fprintf('\t ** fnorm = %0.5g\n\n',fnorm);
    end

    % Step 2. adjust simulation based on net motion of system
    if (P > 0)
        % increase positive counter
        npPos = npPos + 1;

        % reset negative counter
        npNeg = 0;

        % update alphat for the next iteration
        alphat = alpha;

        % alter simulation if enough positive steps have been taken
        if (npPos > NMIN)
            % change time step
            if (dt*finc < dtmax)
                dt = dt*finc;
            end
        end

        % decrease alpha
        alpha = alpha*falpha;
    else
        % reset positive counter
        npPos = 0;

        % increase negative counter
        npNeg = npNeg + 1;

        % check for stuck simulation
        if (npNeg > NNEGMAX)
            fprintf('Simulation negative for too long, ending program here.\n')
            error('FIRE did not converge');
        end

        % decrease time step if past initial delay
        if (it > NMIN)
            % decrease time step
            if (dt*fdec > dtmin)
                dt = dt*fdec;
            end

            % change alpha
            alpha = alpha0;
        end

        % take a half step backwards
        x = x - 0.5*dt*vx;
        y = y - 0.5*dt*vy;

        % reset velocities to 0
        vx = zeros(N,1);
        vy = zeros(N,1);
    end

    % update velocities if forces are acting
    if fnorm > 0
        vx = (1 - alpha).*vx + alpha.*(Fx./fnorm)*vnorm;
        vy = (1 - alpha).*vy + alpha.*(Fy./fnorm)*vnorm;
    end

    % do first verlet update for vertices (assume unit mass)
    x = x + dt*vx;
    y = y + dt*vy;

    % do periodic boundary conditions
    x = mod(x,Lx);
    y = mod(y,Ly);

    % -- Compute new forces
    Fx = zeros(N,1);
    Fy = zeros(N,1);
    U = 0.0;
    vstress = zeros(4,1);
    cij = zeros(N);
    for ii = 1:N
        % particle i information
        xi = x(ii);
        yi = y(ii);
        ri = r(ii);

        % loop over particles j = [i+1,N]
        for jj = ii+1:N
            % sigma ij
            sij = ri + r(jj);

            % x distance
            dx = x(jj) - xi;
            dx = dx - Lx*round(dx/Lx);
            if dx < sij
                % y distance
                dy = y(jj) - yi;
                dy = dy - Ly*round(dy/Ly);

                % calculate full distance
                dist = sqrt(dx*dx + dy*dy);

                % check overlap
                if dist < sij
                    % unit vector
                    ux = dx/dist;
                    uy = dy/dist;

                    % force magnitude
                    ftmp = (1 - dist/sij)/sij;

                    % force in given direction
                    fx = ftmp*ux;
                    fy = ftmp*uy;

                    % update forces
                    Fx(ii) = Fx(ii) - fx;
                    Fy(ii) = Fy(ii) - fy;

                    Fx(jj) = Fx(jj) + fx;
                    Fy(jj) = Fy(jj) + fy;

                    % update virial stress
                    vstress(1) = vstress(1) + fx*dx/(Lx*Ly);    % sigmaXX
                    vstress(2) = vstress(2) + fx*dy/(Lx*Ly);    % sigmaXY
                    vstress(3) = vstress(3) + fy*dx/(Lx*Ly);    % sigmaYX
                    vstress(4) = vstress(4) + fy*dy/(Lx*Ly);    % sigmaYY

                    % update contact network
                    cij(ii,jj) = 1;
                    cij(jj,ii) = 1;

                    % update potential energy
                    U = U + 0.5*(1 - (dist/sij))^2;
                end
            end
        end
    end

    % do second verlet update for vertices
    vx = vx + dt*0.5*Fx;
    vy = vy + dt*0.5*Fy;

    % update Kcheck and Fcheck
    fcheck = sqrt(sum(Fx.*Fx) + sum(Fy.*Fy))/N;
    kcheck = 0.5*sum(vx.*vx + vy.*vy)/N;

    if fcheck < Ftol
        npPMin = npPMin + 1;
        if npPMin >= NMIN && kcheck < Ktol
            fprintf('FIRE has converged in fireit = %d steps...\n',fireit);
            fprintf('fcheck = %0.3g\n',fcheck);
            fprintf('kcheck = %0.3g\n',kcheck);
        end
    else
        npPMin = 0;
    end
end
if fireit == itmax
    fprintf('FIRE protocol did not converge in less than itmax = %d iterations on it = %d, ending.\n',itmax,it);
    error('FIRE protocol did not converge.');
end


% save new jammed state
x2 = x;
y2 = y;
vx2 = vx;
vy2 = vy;
U2 = U;
Fx2 = Fx;
Fy2 = Fy;
vstress2 = vstress;
cij2 = cij;

% initial packing fraction
phi = sum(pi*r.^2)/(Lx*Ly);
phi2 = phi;

%% Oscillatory strain with discontinuous volume fraction change
% FIRST JAMMED STATE

% print to console
fprintf('\n\n\n\n &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n\n');
fprintf('           Now shearing first jammed state\n\n');
fprintf(' &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n\n\n\n');

% reset coordinates to initial jammed state
r   = r1;
x   = x1;
y   = y1;
vx  = vx1;
vy  = vy1;
U   = U1;
Fx  = Fx1;
Fy  = Fy1;

Fxold = Fx;
Fyold = Fy;

% anonymous function for strains
gamma = @(g0,w,t) 0.5*g0*(1 - cos(w*t));
gammaDot = @(g0,w,t) 0.5*g0*w*sin(w*t);
lastStrain = 0;
currStrain = 0;
strainRate = 0;

% time variables
t       = 0;
ss      = 1;
totalT  = NCYCLES*tCycle;
dt      = dt0;
NSKIP   = 500;

if makeAVideo == 1
    vobj    = VideoWriter('strainRate1.mp4','MPEG-4');
    open(vobj);
end

% store data during cyclic shear
NSTEPS          = ceil(totalT/dt);
virialStress1    = zeros(NSTEPS,4);
potEnergy1       = zeros(NSTEPS,1);
kinEnergy1       = zeros(NSTEPS,1);
instStrain1      = zeros(NSTEPS,1);

% loop over oscillatory shear
while t < totalT
    % do first verlet update for vertices
    x = x + dt*vx + 0.5*dt*dt*Fxold;
    y = y + dt*vy + 0.5*dt*dt*Fyold;

    % ADD LEbc (on both positions and velocities)
    im  = floor(y/Ly);
    x   = mod(x,Lx);
    y   = mod(y,Ly);
    x   = x - im*currStrain*Ly;
    vx  = vx - im*strainRate*Ly;

    % -- Compute new forces
    Fx = zeros(N,1);
    Fy = zeros(N,1);
    U = 0.0;
    vstress = zeros(4,1);
    cij = zeros(N);
    for ii = 1:N
        % particle i information
        xi = x(ii);
        yi = y(ii);
        ri = r(ii);

        % loop over particles j = [i+1,N]
        for jj = ii+1:N
            % sigma ij
            sij = ri + r(jj);

            % y distance
            dy = y(jj) - yi;
            im = round(dy/Ly);
            dy = dy - Ly*im;
            if dy < sij
                % x distance
                dx = x(jj) - xi;

                % LEbc
                dx = x(jj) - xi;
                dx = dx - Ly*im*currStrain;
                dx = dx - Lx*round(dx/Lx);

                % calculate full distance
                dist = sqrt(dx*dx + dy*dy);

                % check overlap
                if dist < sij
                    % unit vector
                    ux = dx/dist;
                    uy = dy/dist;

                    % force magnitude
                    ftmp = (1 - dist/sij)/sij;

                    % force in given direction
                    fx = ftmp*ux;
                    fy = ftmp*uy;

                    % update forces
                    Fx(ii) = Fx(ii) - fx;
                    Fy(ii) = Fy(ii) - fy;

                    Fx(jj) = Fx(jj) + fx;
                    Fy(jj) = Fy(jj) + fy;

                    % update virial stress
                    vstress(1) = vstress(1) + fx*dx/(Lx*Ly);    % sigmaXX
                    vstress(2) = vstress(2) + fx*dy/(Lx*Ly);    % sigmaXY
                    vstress(3) = vstress(3) + fy*dx/(Lx*Ly);    % sigmaYX
                    vstress(4) = vstress(4) + fy*dy/(Lx*Ly);    % sigmaYY

                    % update contact network
                    cij(ii,jj) = 1;
                    cij(jj,ii) = 1;

                    % update potential energy
                    U = U + 0.5*(1 - (dist/sij))^2;
                end
            end
        end
    end
    
    % implement damping with velocity verlet
    dampingNumX = bdamp*(vx - 0.5*Fxold*dt);
    dampingNumY = bdamp*(vy - 0.5*Fyold*dt);
    
    dampingDenom = 1.0 - 0.5*bdamp*dt;
    
    Fx = (Fx - dampingNumX)./dampingDenom;
    Fy = (Fy - dampingNumY)./dampingDenom;

    % do second verlet update for vertices
    vx = vx + dt*0.5*(Fx + Fxold);
    vy = vy + dt*0.5*(Fy + Fyold);
    
    Fxold = Fx;
    Fyold = Fy;
    
    % store state functions
    virialStress1(ss,:) = vstress';
    potEnergy1(ss) = U;
    kinEnergy1(ss) = 0.5*sum(vx.^2 + vy.^2);
    
    % output
    if mod(ss,NSKIP) == 0
        
        if makeAVideo == 1
            % Draw system with deformed boundary
            figure(1), clf, hold on, box on;

            % draw particles
            for xx = -1:1
                for yy = -1:1
                    for nn = 1:N
                        xplot = x(nn) - r(nn);
                        yplot = y(nn) - r(nn);
                        rectangle('Position',[xplot + xx*Lx + currStrain*yy*Ly yplot+yy*Ly 2*r(nn) 2*r(nn)],'Curvature',[1 1],'FaceColor',[0 0.3 1],'EdgeColor','k');
                    end
                end
            end

            for xx = -1:1
                for yy = -1:1
                    % draw contacts
                    for nn = 1:N
                        for mm = nn+1:N
                            if cij(nn,mm) == 1
                                dy = y(mm) - y(nn);
                                im = round(dy/Ly);
                                dy = dy - Ly*im;

                                dx = x(mm) - x(nn);
                                dx = dx - Ly*im*currStrain;
                                dx = dx - Lx*round(dx/Lx);

                                line([x(nn) + xx*Lx + currStrain*yy*Ly x(nn) + dx + xx*Lx + currStrain*yy*Ly],[y(nn) + yy*Ly y(nn) + dy + yy*Ly],'linestyle','-','color','k','linewidth',1.5);
                            end
                        end
                    end
                end
            end


            % draw both boundaries
            plot([0 1 1+currStrain currStrain 0]*Lx,[0 0 Ly Ly 0],'r');
            plot([0 1 1 0 0]*Lx,[0 0 Ly Ly 0],'k');

            % axes
            axis('equal'); box on; set(gca,'XTick',[],'YTick',[]);
            axis([-Lx/4 5/4*Lx -Ly/4 5/4*Ly]);

            currFrame = getframe(gcf);
            writeVideo(vobj,currFrame);
        end
        
        % print to console
        fprintf('From jammed state 1: t = %0.3g/%0.3g, ss = %d, gam = %0.3g\n',t,totalT,ss,currStrain);
    end
    
    % increment time
    t = t + dt;
    ss = ss + 1;
    
    % update current strain
    currStrain  = gamma(strainAmp,strainFreq,t);
    dStrain     = currStrain - lastStrain;
    lastStrain  = currStrain;
    
    % save strain
    instStrain1(ss) = currStrain;
    
    % current strain rate
    strainRate  = gammaDot(strainAmp,strainFreq,t);
    
    % affine displacement due to shear strain
    x = x + dStrain.*y;
end

if makeAVideo == 1
    close(vobj);
end


%% Oscillatory strain with discontinuous volume fraction change
% SECOND JAMMED STATE

% print to console
fprintf('\n\n\n\n &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n\n');
fprintf('           Now shearing second jammed state\n\n');
fprintf(' &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n\n\n\n');

% reset coordinates to initial jammed state
r   = r2;
x   = x2;
y   = y2;
vx  = vx2;
vy  = vy2;
U   = U2;
Fx  = Fx2;
Fy  = Fy2;

Fxold = Fx;
Fyold = Fy;

% anonymous function for strains
gamma = @(g0,w,t) 0.5*g0*(1 - cos(w*t));
gammaDot = @(g0,w,t) 0.5*g0*w*sin(w*t);
lastStrain = 0;
currStrain = 0;
strainRate = 0;

% time variables
t       = 0;
ss      = 1;
totalT  = NCYCLES*tCycle;
dt      = dt0;
NSKIP   = 500;

if makeAVideo == 1
    vobj    = VideoWriter('strainRate2.mp4','MPEG-4');
    open(vobj);
end

% store data during cyclic shear
NSTEPS          = ceil(totalT/dt);
virialStress2    = zeros(NSTEPS,4);
potEnergy2       = zeros(NSTEPS,1);
kinEnergy2       = zeros(NSTEPS,1);
instStrain2      = zeros(NSTEPS,1);

% loop over oscillatory shear
while t < totalT
    % do first verlet update for vertices
    x = x + dt*vx + 0.5*dt*dt*Fxold;
    y = y + dt*vy + 0.5*dt*dt*Fyold;

    % ADD LEbc (on both positions and velocities)
    im  = floor(y/Ly);
    x   = mod(x,Lx);
    y   = mod(y,Ly);
    x   = x - im*currStrain*Ly;
    vx  = vx - im*strainRate*Ly;

    % -- Compute new forces
    Fx = zeros(N,1);
    Fy = zeros(N,1);
    U = 0.0;
    vstress = zeros(4,1);
    cij = zeros(N);
    for ii = 1:N
        % particle i information
        xi = x(ii);
        yi = y(ii);
        ri = r(ii);

        % loop over particles j = [i+1,N]
        for jj = ii+1:N
            % sigma ij
            sij = ri + r(jj);

            % y distance
            dy = y(jj) - yi;
            im = round(dy/Ly);
            dy = dy - Ly*im;
            if dy < sij
                % x distance
                dx = x(jj) - xi;

                % LEbc
                dx = x(jj) - xi;
                dx = dx - Ly*im*currStrain;
                dx = dx - Lx*round(dx/Lx);

                % calculate full distance
                dist = sqrt(dx*dx + dy*dy);

                % check overlap
                if dist < sij
                    % unit vector
                    ux = dx/dist;
                    uy = dy/dist;

                    % force magnitude
                    ftmp = (1 - dist/sij)/sij;

                    % force in given direction
                    fx = ftmp*ux;
                    fy = ftmp*uy;

                    % update forces
                    Fx(ii) = Fx(ii) - fx;
                    Fy(ii) = Fy(ii) - fy;

                    Fx(jj) = Fx(jj) + fx;
                    Fy(jj) = Fy(jj) + fy;

                    % update virial stress
                    vstress(1) = vstress(1) + fx*dx/(Lx*Ly);    % sigmaXX
                    vstress(2) = vstress(2) + fx*dy/(Lx*Ly);    % sigmaXY
                    vstress(3) = vstress(3) + fy*dx/(Lx*Ly);    % sigmaYX
                    vstress(4) = vstress(4) + fy*dy/(Lx*Ly);    % sigmaYY

                    % update contact network
                    cij(ii,jj) = 1;
                    cij(jj,ii) = 1;

                    % update potential energy
                    U = U + 0.5*(1 - (dist/sij))^2;
                end
            end
        end
    end
    
    % implement damping with velocity verlet
    dampingNumX = bdamp*(vx - 0.5*Fxold*dt);
    dampingNumY = bdamp*(vy - 0.5*Fyold*dt);
    
    dampingDenom = 1.0 - 0.5*bdamp*dt;
    
    Fx = (Fx - dampingNumX)./dampingDenom;
    Fy = (Fy - dampingNumY)./dampingDenom;

    % do second verlet update for vertices
    vx = vx + dt*0.5*(Fx + Fxold);
    vy = vy + dt*0.5*(Fy + Fyold);
    
    Fxold = Fx;
    Fyold = Fy;
    
    % store state functions
    virialStress2(ss,:) = vstress';
    potEnergy2(ss) = U;
    kinEnergy2(ss) = 0.5*sum(vx.^2 + vy.^2);
    
    % output
    if mod(ss,NSKIP) == 0
        
        if makeAVideo == 1
            % Draw system with deformed boundary
            figure(1), clf, hold on, box on;

            % draw particles
            for xx = -1:1
                for yy = -1:1
                    for nn = 1:N
                        xplot = x(nn) - r(nn);
                        yplot = y(nn) - r(nn);
                        rectangle('Position',[xplot + xx*Lx + currStrain*yy*Ly yplot+yy*Ly 2*r(nn) 2*r(nn)],'Curvature',[1 1],'FaceColor',[0 0.3 1],'EdgeColor','k');
                    end
                end
            end

            for xx = -1:1
                for yy = -1:1
                    % draw contacts
                    for nn = 1:N
                        for mm = nn+1:N
                            if cij(nn,mm) == 1
                                dy = y(mm) - y(nn);
                                im = round(dy/Ly);
                                dy = dy - Ly*im;

                                dx = x(mm) - x(nn);
                                dx = dx - Ly*im*currStrain;
                                dx = dx - Lx*round(dx/Lx);


                                if (im ~= 0 && xx == 0)
                                    test = 1;
                                end

                                line([x(nn) + xx*Lx + currStrain*yy*Ly x(nn) + dx + xx*Lx + currStrain*yy*Ly],[y(nn) + yy*Ly y(nn) + dy + yy*Ly],'linestyle','-','color','k','linewidth',1.5);
                            end
                        end
                    end
                end
            end


            % draw both boundaries
            plot([0 1 1+currStrain currStrain 0]*Lx,[0 0 Ly Ly 0],'r');
            plot([0 1 1 0 0]*Lx,[0 0 Ly Ly 0],'k');

            % axes
            axis('equal'); box on; set(gca,'XTick',[],'YTick',[]);
            axis([-Lx/4 5/4*Lx -Ly/4 5/4*Ly]);

            currFrame = getframe(gcf);
            writeVideo(vobj,currFrame);
        end
        
        % print to console
        fprintf('From jammed state 2: t = %0.3g/%0.3g, ss = %d, gam = %0.3g\n',t,totalT,ss,currStrain);
    end
    
    % increment time
    t = t + dt;
    ss = ss + 1;
    
    % update current strain
    currStrain  = gamma(strainAmp,strainFreq,t);
    dStrain     = currStrain - lastStrain;
    lastStrain  = currStrain;
    
    % save strain
    instStrain2(ss) = currStrain;
    
    % current strain rate
    strainRate  = gammaDot(strainAmp,strainFreq,t);
    
    % affine displacement due to shear strain
    x = x + dStrain.*y;
end

if makeAVideo == 1
    close(vobj);
end

%% Oscillatory strain with discontinuous volume fraction change
% WITH VOLUME FRACTION CHANGE

% print to console
fprintf('\n\n\n\n &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n\n');
fprintf('           Now shearing with phi change impulse to first jammed state\n\n');
fprintf(' &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n\n\n\n');

% reset coordinates to initial jammed state
r   = r1;
x   = x1;
y   = y1;
vx  = vx1;
vy  = vy1;
U   = U1;
Fx  = Fx1;
Fy  = Fy1;
phi = phi1;

Fxold = Fx;
Fyold = Fy;

% anonymous function for strains
gamma = @(g0,w,t) 0.5*g0*(1 - cos(w*t));
gammaDot = @(g0,w,t) 0.5*g0*w*sin(w*t);
lastStrain = 0;
currStrain = 0;
strainRate = 0;

% time variables
t       = 0;
ss      = 1;
totalT  = NCYCLES*tCycle;
dt      = dt0;
NSKIP   = 500;

% can update makeAVideo here if desired
makeAVideo = 0;

if makeAVideo == 1
    vobj    = VideoWriter('strainRateImpulse.mp4','MPEG-4');
    open(vobj);
end

% store data during cyclic shear
NSTEPS              = ceil(totalT/dt);
virialStress3       = zeros(NSTEPS,4);
potEnergy3          = zeros(NSTEPS,1);
kinEnergy3          = zeros(NSTEPS,1);
instStrain3         = zeros(NSTEPS,1);

% phi impulse variables
tOn                 = 0.25*totalT;
tOff                = 0.75*totalT;
instphi             = zeros(NSTEPS,1);

% loop over oscillatory shear
while t < totalT
    % do first verlet update for vertices
    x = x + dt*vx + 0.5*dt*dt*Fxold;
    y = y + dt*vy + 0.5*dt*dt*Fyold;

    % ADD LEbc (on both positions and velocities)
    im  = floor(y/Ly);
    x   = mod(x,Lx);
    y   = mod(y,Ly);
    x   = x - im*currStrain*Ly;
    vx  = vx - im*strainRate*Ly;

    % -- Compute new forces
    Fx = zeros(N,1);
    Fy = zeros(N,1);
    U = 0.0;
    vstress = zeros(4,1);
    cij = zeros(N);
    for ii = 1:N
        % particle i information
        xi = x(ii);
        yi = y(ii);
        ri = r(ii);

        % loop over particles j = [i+1,N]
        for jj = ii+1:N
            % sigma ij
            sij = ri + r(jj);

            % y distance
            dy = y(jj) - yi;
            im = round(dy/Ly);
            dy = dy - Ly*im;
            if dy < sij
                % x distance
                dx = x(jj) - xi;

                % LEbc
                dx = x(jj) - xi;
                dx = dx - Ly*im*currStrain;
                dx = dx - Lx*round(dx/Lx);

                % calculate full distance
                dist = sqrt(dx*dx + dy*dy);

                % check overlap
                if dist < sij
                    % unit vector
                    ux = dx/dist;
                    uy = dy/dist;

                    % force magnitude
                    ftmp = (1 - dist/sij)/sij;

                    % force in given direction
                    fx = ftmp*ux;
                    fy = ftmp*uy;

                    % update forces
                    Fx(ii) = Fx(ii) - fx;
                    Fy(ii) = Fy(ii) - fy;

                    Fx(jj) = Fx(jj) + fx;
                    Fy(jj) = Fy(jj) + fy;

                    % update virial stress
                    vstress(1) = vstress(1) + fx*dx/(Lx*Ly);    % sigmaXX
                    vstress(2) = vstress(2) + fx*dy/(Lx*Ly);    % sigmaXY
                    vstress(3) = vstress(3) + fy*dx/(Lx*Ly);    % sigmaYX
                    vstress(4) = vstress(4) + fy*dy/(Lx*Ly);    % sigmaYY

                    % update contact network
                    cij(ii,jj) = 1;
                    cij(jj,ii) = 1;

                    % update potential energy
                    U = U + 0.5*(1 - (dist/sij))^2;
                end
            end
        end
    end
    
    % implement damping with velocity verlet
    dampingNumX = bdamp*(vx - 0.5*Fxold*dt);
    dampingNumY = bdamp*(vy - 0.5*Fyold*dt);
    
    dampingDenom = 1.0 - 0.5*bdamp*dt;
    
    Fx = (Fx - dampingNumX)./dampingDenom;
    Fy = (Fy - dampingNumY)./dampingDenom;

    % do second verlet update for vertices
    vx = vx + dt*0.5*(Fx + Fxold);
    vy = vy + dt*0.5*(Fy + Fyold);
    
    Fxold = Fx;
    Fyold = Fy;
    
    % store state functions
    virialStress3(ss,:) = vstress';
    potEnergy3(ss) = U;
    kinEnergy3(ss) = 0.5*sum(vx.^2 + vy.^2);
    
    % output
    if mod(ss,NSKIP) == 0
        
        if makeAVideo == 1
            % Draw system with deformed boundary
            figure(1), clf, hold on, box on;

            % draw particles
            for xx = -1:1
                for yy = -1:1
                    for nn = 1:N
                        xplot = x(nn) - r(nn);
                        yplot = y(nn) - r(nn);
                        rectangle('Position',[xplot + xx*Lx + currStrain*yy*Ly yplot+yy*Ly 2*r(nn) 2*r(nn)],'Curvature',[1 1],'FaceColor',[0 0.3 1],'EdgeColor','k');
                    end
                end
            end

            for xx = -1:1
                for yy = -1:1
                    % draw contacts
                    for nn = 1:N
                        for mm = nn+1:N
                            if cij(nn,mm) == 1
                                dy = y(mm) - y(nn);
                                im = round(dy/Ly);
                                dy = dy - Ly*im;

                                dx = x(mm) - x(nn);
                                dx = dx - Ly*im*currStrain;
                                dx = dx - Lx*round(dx/Lx);


                                if (im ~= 0 && xx == 0)
                                    test = 1;
                                end

                                line([x(nn) + xx*Lx + currStrain*yy*Ly x(nn) + dx + xx*Lx + currStrain*yy*Ly],[y(nn) + yy*Ly y(nn) + dy + yy*Ly],'linestyle','-','color','k','linewidth',1.5);
                            end
                        end
                    end
                end
            end


            % draw both boundaries
            plot([0 1 1+currStrain currStrain 0]*Lx,[0 0 Ly Ly 0],'r');
            plot([0 1 1 0 0]*Lx,[0 0 Ly Ly 0],'k');

            % axes
            axis('equal'); box on; set(gca,'XTick',[],'YTick',[]);
            axis([-Lx/4 5/4*Lx -Ly/4 5/4*Ly]);

            currFrame = getframe(gcf);
            writeVideo(vobj,currFrame);
        end
        
        % print to console
        fprintf('Impulse protocol: phi %0.3g, t = %0.3g/%0.3g, ss = %d, gam = %0.3g\n',phi,t,totalT,ss,currStrain);
    end
    
    % set radii based on t
    if t > tOn && t < tOff
        r = r2;
        instphi(ss) = phi2;
        phi = phi2;
    else
        r = r1;
        instphi(ss) = phi1;
        phi = phi1;
    end
    
    % increment time
    t = t + dt;
    ss = ss + 1;
    
    % update current strain
    currStrain  = gamma(strainAmp,strainFreq,t);
    dStrain     = currStrain - lastStrain;
    lastStrain  = currStrain;
    
    % save strain
    instStrain3(ss) = currStrain;
    
    
    % current strain rate
    strainRate  = gammaDot(strainAmp,strainFreq,t);
    
    % affine displacement due to shear strain
    x = x + dStrain.*y;
end

if makeAVideo == 1
    close(vobj);
end

%% plot things

% virial stress during strain
vp1 = 0.5*(virialStress1(:,1) + virialStress1(:,4))./(Lx*Ly);
vp2 = 0.5*(virialStress2(:,1) + virialStress2(:,4))./(Lx*Ly);
vp3 = 0.5*(virialStress3(:,1) + virialStress3(:,4))./(Lx*Ly);
vss1 = virialStress1(:,2)./(Lx*Ly);
vss2 = virialStress2(:,2)./(Lx*Ly);
vss3 = virialStress3(:,2)./(Lx*Ly);

% time
timeVals = (0:dt:totalT)./tCycle;
cycleInds = floor(timeVals);

% stress over time
figure(2), clf, hold on, box on;
plot(timeVals,vss1,'k:','linewidth',1.2);
plot(timeVals,vss2,'k--','linewidth',1.2);
plot(timeVals,vss3,'k-','linewidth',1.5);
xlabel('$2\pi t/\omega$','Interpreter','latex');
ylabel('$\sigma_{xy}$','Interpreter','latex');
ax = gca;
ax.FontSize = 22;

figure(3), clf, hold on, box on;
plot(timeVals,vp1,'k:','linewidth',1.2);
plot(timeVals,vp2,'k--','linewidth',1.2);
plot(timeVals,vp3,'k-','linewidth',1.5);
xlabel('$2\pi t/\omega$','Interpreter','latex');
ylabel('$p$','Interpreter','latex');
ax = gca;
ax.FontSize = 22;

% potential energy over time
figure(4), clf, hold on, box on;
plot(timeVals,potEnergy1,'b:','linewidth',1.2);
plot(timeVals,potEnergy2,'b--','linewidth',1.2);
plot(timeVals,potEnergy3,'b-','linewidth',1.5);
xlabel('$2\pi t/\omega$','Interpreter','latex');
ylabel('$U$','Interpreter','latex');
ax = gca;
ax.FontSize = 22;


% plot difference in potential energies over time
du1 = abs(potEnergy3 - potEnergy1);
du2 = abs(potEnergy3 - potEnergy2);
figure(5), clf, hold on, box on;
plot(timeVals,du1,'b-','linewidth',1.2);
plot(timeVals,du2,'b-.','linewidth',1.2);
xlabel('$2\pi t/\omega$','Interpreter','latex');
ylabel('$\Delta U$','Interpreter','latex');
ax = gca;
ax.FontSize = 22;
ax.YScale = 'log';
legend({'$|U_{\rm test} - U_1|$','$|U_{\rm test} - U_2|$'},'Interpreter','latex',...
    'fontsize',16,'location','best');



% plot shear stress vs strain over time
NS = length(vss1);
figure(10), clf, hold on, box on;
plot(instStrain3(1:NS),vss3,'k--','linewidth',1.2);
plot(instStrain1(1:NS),vss1,'r-','linewidth',2);
xlabel('$\gamma$','Interpreter','latex');
ylabel('$\sigma_{xy}$','Interpreter','latex');
ax = gca;
ax.FontSize = 22;
legend({'test','packing $1$'},'Interpreter','latex',...
    'fontsize',16,'location','best');

figure(11), clf, hold on, box on;
plot(instStrain3(1:NS),vss3,'k--','linewidth',1.2);
plot(instStrain2(1:NS),vss2,'g-','linewidth',2);
xlabel('$\gamma$','Interpreter','latex');
ylabel('$\sigma_{xy}$','Interpreter','latex');
ax = gca;
ax.FontSize = 22;
legend({'test','packing $2$'},'Interpreter','latex',...
    'fontsize',16,'location','best');

% plot phi impulse over time
figure(20), clf, hold on, box on;
plot(timeVals,instphi,'g-','linewidth',2.5);
xlabel('$2\pi t/\omega$','Interpreter','latex');
ylabel('$\phi$','Interpreter','latex');
ax = gca;
ax.FontSize = 22;