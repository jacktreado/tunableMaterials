%% Simulation to do 3 strain rate sims

clear;
close all;
clc;

%% Initialize system

N               = 48;       % number of particles
phi0            = 0.50;     % initial packing fraction
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
NCYCLES         = 20;       % number of strain cycles
bdamp           = 1.0;      % damping parameter

dphi            = 0.1;     % particle size increase

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


% Save variables for animation at the end
NSTEPS     = ceil(totalT/dt);
Xanimation = zeros(N,NSTEPS,3);
Yanimation = zeros(N,NSTEPS,3);
Ranimation = zeros(N,NSTEPS,3);


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


% store data during cyclic shear
virialStress1    = zeros(NSTEPS,4);
potEnergy1       = zeros(NSTEPS,1);
kinEnergy1       = zeros(NSTEPS,1);
instStrain      = zeros(NSTEPS,1);

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
        % print to console
        fprintf('From jammed state 1: t = %0.3g/%0.3g, ss = %d, gam = %0.3g\n',t,totalT,ss,currStrain);
    end
    
    %Save variables for animation
    Xanimation(:,ss,1) = x;
    Yanimation(:,ss,1) = y;
    Ranimation(:,ss,1) = r;
    
    % increment time
    t = t + dt;
    ss = ss + 1;
    
    % update current strain
    currStrain  = gamma(strainAmp,strainFreq,t);
    dStrain     = currStrain - lastStrain;
    lastStrain  = currStrain;
    
    % save strain
    instStrain(ss) = currStrain;
    
    % current strain rate
    strainRate  = gammaDot(strainAmp,strainFreq,t);
    
    % affine displacement due to shear strain
    x = x + dStrain.*y;
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


% time variables
t       = 0;
ss      = 1;
dt      = dt0;

% store data during cyclic shear
virialStress2    = zeros(NSTEPS,4);
potEnergy2       = zeros(NSTEPS,1);
kinEnergy2       = zeros(NSTEPS,1);

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
        
        % print to console
        fprintf('From jammed state 2: t = %0.3g/%0.3g, ss = %d, gam = %0.3g\n',t,totalT,ss,currStrain);
    end
    
    %Save variables for animation
    Xanimation(:,ss,2) = x;
    Yanimation(:,ss,2) = y;
    Ranimation(:,ss,2) = r;
    
    % increment time
    t = t + dt;
    ss = ss + 1;
    
    % update current strain
    currStrain  = gamma(strainAmp,strainFreq,t);
    dStrain     = currStrain - lastStrain;
    lastStrain  = currStrain;
    
    % current strain rate
    strainRate  = gammaDot(strainAmp,strainFreq,t);
    
    % affine displacement due to shear strain
    x = x + dStrain.*y;
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
lastStrain = 0;
currStrain = 0;
strainRate = 0;

% time variables
t       = 0;
ss      = 1;
dt      = dt0;

% store data during cyclic shear
virialStress3       = zeros(NSTEPS,4);
potEnergy3          = zeros(NSTEPS,1);
kinEnergy3          = zeros(NSTEPS,1);

% phi impulse variables
tActivation         = pi/3* (tCycle/2/pi);
tDeactivation       = 4*pi/3 * (tCycle/2/pi);
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
        
        % print to console
        fprintf('Impulse protocol: phi %0.3g, t = %0.3g/%0.3g, ss = %d, gam = %0.3g\n',phi,t,totalT,ss,currStrain);
    end
    
    %Save variables for animation
    Xanimation(:,ss,3) = x;
    Yanimation(:,ss,3) = y;
    Ranimation(:,ss,3) = r;
    
    % set radii based on t
    if mod(t-tActivation,tCycle) >= 0  && mod(t-tActivation,tCycle) < mod(tDeactivation-tActivation,tCycle)
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
    
    
    % current strain rate
    strainRate  = gammaDot(strainAmp,strainFreq,t);
    
    % affine displacement due to shear strain
    x = x + dStrain.*y;
end

%% Stress calculation
% virial stress during strain
vp1 = 0.5*(virialStress1(:,1) + virialStress1(:,4))./(Lx*Ly);
vp2 = 0.5*(virialStress2(:,1) + virialStress2(:,4))./(Lx*Ly);
vp3 = 0.5*(virialStress3(:,1) + virialStress3(:,4))./(Lx*Ly);
vss1 = virialStress1(:,2)./(Lx*Ly);
vss2 = virialStress2(:,2)./(Lx*Ly);
vss3 = virialStress3(:,2)./(Lx*Ly);






%% Animations
% time
timeVals = (0:dt:totalT)./tCycle;
cycleInds = floor(timeVals);

stop = length(vss1);
start = stop - 10*ceil(tCycle/dt);

%Animate loops
handle = figure(1);clf
set(handle,'Position',[0,0,1080,360],'units','points')



%Particle animation
subplot(1,3,1)
title(sprintf('phi = %2.3f' , phi1 ));

axis equal
box on; set(gca,'XTick',[],'YTick',[]);
axis([-Lx/4 5/4*Lx -Ly/4 5/4*Ly]);

AXES = gca;
scale = handle.Position(3)*AXES.Position(3)/(AXES.XLim(2) - AXES.XLim(1));
RMax1 = max(Ranimation(:,1,1));
RMin1 = min(Ranimation(:,1,1));

ParticleRMinEdgeSplice1 = animatedline;
ParticleRMaxEdgeSplice1 = animatedline;
ParticleRMinCenterSplice1 = animatedline;
ParticleRMaxCenterSplice1 = animatedline;
set(ParticleRMaxEdgeSplice1  ,'Marker','o','MarkerEdgeColor','k','linestyle','none')
set(ParticleRMinEdgeSplice1  ,'Marker','o','MarkerEdgeColor','k','linestyle','none')
set(ParticleRMaxCenterSplice1,'Marker','+','MarkerEdgeColor','k','Markersize',8,'linestyle','none')
set(ParticleRMinCenterSplice1,'Marker','+','MarkerEdgeColor','k','Markersize',8,'linestyle','none')

ParticleRMinEdge1 = animatedline;
ParticleRMaxEdge1 = animatedline;
ParticleRMinCenter1 = animatedline;
ParticleRMaxCenter1 = animatedline;
set(ParticleRMaxEdge1  ,'Marker','o','MarkerEdgeColor','r','Markersize',2*RMax1*scale,'linestyle','none')
set(ParticleRMinEdge1  ,'Marker','o','MarkerEdgeColor','r','Markersize',2*RMin1*scale,'linestyle','none')
set(ParticleRMaxCenter1,'Marker','+','MarkerEdgeColor','r','Markersize',4,'linestyle','none')
set(ParticleRMinCenter1,'Marker','+','MarkerEdgeColor','r','Markersize',4,'linestyle','none')


Boundary1 = animatedline;
set(Boundary1,'MaximumNumPoints',5,'color','m','linewidth',2)



subplot(1,3,2)
title(sprintf('phi = %2.3f' , phi2 ));

axis equal
box on; set(gca,'XTick',[],'YTick',[]);
axis([-Lx/4 5/4*Lx -Ly/4 5/4*Ly]);

RMax2 = max(Ranimation(:,1,2));
RMin2 = min(Ranimation(:,1,2));

ParticleRMinEdgeSplice2 = animatedline;
ParticleRMaxEdgeSplice2 = animatedline;
ParticleRMinCenterSplice2 = animatedline;
ParticleRMaxCenterSplice2 = animatedline;
set(ParticleRMaxEdgeSplice2  ,'Marker','o','MarkerEdgeColor','k','linestyle','none')
set(ParticleRMinEdgeSplice2  ,'Marker','o','MarkerEdgeColor','k','linestyle','none')
set(ParticleRMaxCenterSplice2,'Marker','+','MarkerEdgeColor','k','Markersize',8,'linestyle','none')
set(ParticleRMinCenterSplice2,'Marker','+','MarkerEdgeColor','k','Markersize',8,'linestyle','none')

ParticleRMinEdge2 = animatedline;
ParticleRMaxEdge2 = animatedline;
ParticleRMinCenter2 = animatedline;
ParticleRMaxCenter2 = animatedline;
set(ParticleRMaxEdge2  ,'Marker','o','MarkerEdgeColor','g','Markersize',2*RMax2*scale,'linestyle','none')
set(ParticleRMinEdge2  ,'Marker','o','MarkerEdgeColor','g','Markersize',2*RMin2*scale,'linestyle','none')
set(ParticleRMaxCenter2,'Marker','+','MarkerEdgeColor','g','Markersize',4,'linestyle','none')
set(ParticleRMinCenter2,'Marker','+','MarkerEdgeColor','g','Markersize',4,'linestyle','none')

Boundary2 = animatedline;
set(Boundary2,'MaximumNumPoints',5,'color','m','linewidth',2)

%Loop Animation
subplot(1,3,3)
miny = min( [min(vss1),min(vss2),min(vss3)]);
maxy = max( [max(vss1),max(vss2),max(vss3)]);

xlim([min(instStrain)-0.01,max(instStrain)+0.01])
ylim([miny,maxy])
MaxNumPoints = 5*round(stop/NCYCLES);
h1 = animatedline; set(h1,'color','r','linewidth',2, 'MaximumNumPoints',MaxNumPoints);
h2 = animatedline; set(h2,'color','g','linewidth',2, 'MaximumNumPoints',MaxNumPoints);
h3 = animatedline; set(h3,'color','k','linewidth',2, 'MaximumNumPoints',MaxNumPoints);

g1 = animatedline; set(g1,'color','r','marker','o','MarkerFaceColor','r','Markersize',10,'MaximumNumPoints',1);
g2 = animatedline; set(g2,'color','g','marker','o','MarkerFaceColor','g','Markersize',10,'MaximumNumPoints',1);
g3 = animatedline; set(g3,'color','k','marker','o','MarkerFaceColor','k','Markersize',10,'MaximumNumPoints',1);


% AnimateVideo = VideoWriter('SplicingAnimation.avi','Motion JPEG AVI');
% open(AnimateVideo)


i = 1;
drawjump = round(stop/NCYCLES/30);
while i <= stop-drawjump
    currStrain  = gamma(strainAmp,strainFreq,i*dt);
    
    %Update animation
    clearpoints(ParticleRMaxEdge1)
    clearpoints(ParticleRMinEdge1)
    clearpoints(ParticleRMaxCenter1)
    clearpoints(ParticleRMinCenter1)
    
    clearpoints(ParticleRMaxEdgeSplice1)
    clearpoints(ParticleRMinEdgeSplice1)
    clearpoints(ParticleRMaxCenterSplice1)
    clearpoints(ParticleRMinCenterSplice1)
    
    clearpoints(ParticleRMaxEdge2)
    clearpoints(ParticleRMinEdge2)
    clearpoints(ParticleRMaxCenter2)
    clearpoints(ParticleRMinCenter2)
    
    clearpoints(ParticleRMaxEdgeSplice2)
    clearpoints(ParticleRMinEdgeSplice2)
    clearpoints(ParticleRMaxCenterSplice2)
    clearpoints(ParticleRMinCenterSplice2)
    
    RMax3 = max(Ranimation(:,i+drawjump-1,3));
    RMin3 = min(Ranimation(:,i+drawjump-1,3));
    set(ParticleRMaxEdgeSplice1  ,'Markersize',2*RMax3*scale)
    set(ParticleRMinEdgeSplice1  ,'Markersize',2*RMin3*scale)
    set(ParticleRMaxEdgeSplice2  ,'Markersize',2*RMax3*scale)
    set(ParticleRMinEdgeSplice2  ,'Markersize',2*RMin3*scale)
    
    bool1 = Ranimation(:,i+drawjump-1,1) == RMax1;
    bool2 = Ranimation(:,i+drawjump-1,2) == RMax2;
    bool3 = Ranimation(:,i+drawjump-1,3) == RMax3;
    
    for xx = -1:1
        for yy = -1:1
            %Animate state 1
            addpoints(ParticleRMaxEdge1,Xanimation(bool1,i+drawjump-1,1)+xx*Lx + currStrain*yy*Ly,...
                Yanimation(bool1,i+drawjump-1,1)+yy*Ly);
            addpoints(ParticleRMinEdge1,Xanimation(~bool1,i+drawjump-1,1)+xx*Lx + currStrain*yy*Ly,...
                Yanimation(~bool1,i+drawjump-1,1)+yy*Ly);
            addpoints(ParticleRMaxCenter1,Xanimation(bool1,i+drawjump-1,1)+xx*Lx + currStrain*yy*Ly,...
                Yanimation(bool1,i+drawjump-1,1)+yy*Ly);
            addpoints(ParticleRMinCenter1,Xanimation(~bool1,i+drawjump-1,1)+xx*Lx + currStrain*yy*Ly,...
                Yanimation(~bool1,i+drawjump-1,1)+yy*Ly);
            
            %Animate state 2
            addpoints(ParticleRMaxEdge2,Xanimation(bool2,i+drawjump-1,1)+xx*Lx + currStrain*yy*Ly,...
                Yanimation(bool2,i+drawjump-1,1)+yy*Ly);
            addpoints(ParticleRMinEdge2,Xanimation(~bool2,i+drawjump-1,1)+xx*Lx + currStrain*yy*Ly,...
                Yanimation(~bool2,i+drawjump-1,1)+yy*Ly);
            addpoints(ParticleRMaxCenter2,Xanimation(bool2,i+drawjump-1,1)+xx*Lx + currStrain*yy*Ly,...
                Yanimation(bool2,i+drawjump-1,1)+yy*Ly);
            addpoints(ParticleRMinCenter2,Xanimation(~bool2,i+drawjump-1,1)+xx*Lx + currStrain*yy*Ly,...
                Yanimation(~bool2,i+drawjump-1,1)+yy*Ly);
            
            %Animate splicing in each subplot
            addpoints(ParticleRMaxEdgeSplice1,Xanimation(bool3,i+drawjump-1,1)+xx*Lx + currStrain*yy*Ly,...
                Yanimation(bool3,i+drawjump-1,1)+yy*Ly);
            addpoints(ParticleRMinEdgeSplice1,Xanimation(~bool3,i+drawjump-1,1)+xx*Lx + currStrain*yy*Ly,...
                Yanimation(~bool3,i+drawjump-1,1)+yy*Ly);
            addpoints(ParticleRMaxCenterSplice1,Xanimation(bool3,i+drawjump-1,1)+xx*Lx + currStrain*yy*Ly,...
                Yanimation(bool3,i+drawjump-1,1)+yy*Ly);
            addpoints(ParticleRMinCenterSplice1,Xanimation(~bool3,i+drawjump-1,1)+xx*Lx + currStrain*yy*Ly,...
                Yanimation(~bool3,i+drawjump-1,1)+yy*Ly);
            
            addpoints(ParticleRMaxEdgeSplice2,Xanimation(bool3,i+drawjump-1,1)+xx*Lx + currStrain*yy*Ly,...
                Yanimation(bool3,i+drawjump-1,1)+yy*Ly);
            addpoints(ParticleRMinEdgeSplice2,Xanimation(~bool3,i+drawjump-1,1)+xx*Lx + currStrain*yy*Ly,...
                Yanimation(~bool3,i+drawjump-1,1)+yy*Ly);
            addpoints(ParticleRMaxCenterSplice2,Xanimation(bool3,i+drawjump-1,1)+xx*Lx + currStrain*yy*Ly,...
                Yanimation(bool3,i+drawjump-1,1)+yy*Ly);
            addpoints(ParticleRMinCenterSplice2,Xanimation(~bool3,i+drawjump-1,1)+xx*Lx + currStrain*yy*Ly,...
                Yanimation(~bool3,i+drawjump-1,1)+yy*Ly);
        end
    end
    addpoints(Boundary1,[0 1 1+currStrain currStrain 0 1]*Lx,[0 0 1 1 0 0]*Ly)
    addpoints(Boundary2,[0 1 1+currStrain currStrain 0 1]*Lx,[0 0 1 1 0 0]*Ly)
    
    
    
    addpoints(h1,instStrain(i:i+drawjump-1),vss1(i:i+drawjump-1))
    addpoints(h2,instStrain(i:i+drawjump-1),vss2(i:i+drawjump-1))
    addpoints(h3,instStrain(i:i+drawjump-1),vss3(i:i+drawjump-1))
    addpoints(g1,instStrain(i+drawjump-1),vss1(i+drawjump-1));
    addpoints(g2,instStrain(i+drawjump-1),vss2(i+drawjump-1));
    addpoints(g3,instStrain(i+drawjump-1),vss3(i+drawjump-1));
    
    title(sprintf('Loop # %d' , round(i*NCYCLES/stop) ));
    
    drawnow
    i = i + drawjump;
    
    if exist('AnimateVideo','var'), writeVideo(AnimateVideo,getframe(gcf)); end
end

if exist('AnimateVideo','var')
    close(AnimateVideo); clear AnimateVideo
end
%
% %
% % draw particles using curved rectangles
%
%
% % Draw system with deformed boundary
% figure(2), clf, hold on, box on;
%
% % draw particles
% for xx = -1:1
%     for yy = -1:1
%         for nn = 1:N
%             xplot = x(nn) - r(nn);
%             yplot = y(nn) - r(nn);
%             rectangle('Position',[xplot + xx*Lx + currStrain*yy*Ly yplot+yy*Ly 2*r(nn) 2*r(nn)],'Curvature',[1 1],'FaceColor',[0 0.3 1],'EdgeColor','k');
%         end
%     end
% end
%
% for xx = -1:1
%     for yy = -1:1
%         % draw contacts
%         for nn = 1:N
%             for mm = nn+1:N
%                 if cij(nn,mm) == 1
%                     dy = y(mm) - y(nn);
%                     im = round(dy/Ly);
%                     dy = dy - Ly*im;
%
%                     dx = x(mm) - x(nn);
%                     dx = dx - Ly*im*currStrain;
%                     dx = dx - Lx*round(dx/Lx);
%
%
%                     if (im ~= 0 && xx == 0)
%                         test = 1;
%                     end
%
%                     line([x(nn) + xx*Lx + currStrain*yy*Ly x(nn) + dx + xx*Lx + currStrain*yy*Ly],[y(nn) + yy*Ly y(nn) + dy + yy*Ly],'linestyle','-','color','k','linewidth',1.5);
%                 end
%             end
%         end
%     end
% end
%
%
% % draw both boundaries
% plot([0 1 1+currStrain currStrain 0]*Lx,[0 0 Ly Ly 0],'r');
% plot([0 1 1 0 0]*Lx,[0 0 Ly Ly 0],'k');
%
% % axes
% axis('equal'); box on; set(gca,'XTick',[],'YTick',[]);
% axis([-Lx/4 5/4*Lx -Ly/4 5/4*Ly]);
