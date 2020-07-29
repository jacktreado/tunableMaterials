%% Simulation to do 3 strain rate sims

clear;
close all;
clc;

%% Initialize system

N               = 48;       % number of particles
phi0            = 0.75;     % initial packing fraction
dphi0           = 0.001;    % initial delta phi
dt0             = 0.01;     % time step magnitude
sr              = 1.4;      % size ratio
seed            = 1;        % random number seed
T               = 1;        % initial temperature

Ftol            = 1e-12;    % force tolerance
Ktol            = 1e-25;    % kinetic tolerance
Ptol            = 1e-8;     % pressure tolerance

strainAmp       = 1e-1;     % strain amplitude
strainFreq      = 5e-2;     % strain angular frequency
NCYCLES         = 100;       % number of strain cycles
bdamp           = 1.0;      % damping parameter

dphi1            = -0.01;    % particle size change for the first state
dphi2            = 0.05;    % particle size change for the second state

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

% initial packing fraction
phi = sum(pi*r.^2)/(Lx*Ly);

%% Create two states by changing phi by dphi1 and dphi2. And then relaxing
% print to console
fprintf('\n\n\n\n &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n\n');
fprintf('           Now compressing new states\n\n');
fprintf(' &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n\n\n\n');

% FIRE VARIABLES
dtmax       = 5*dt0;
dtmin       = 0.01*dt0;
dtstart     = dt;

%copy state variables to a two column vector
X   = repmat(x,[1,2]);
Y   = repmat(y,[1,2]);

VX  = repmat(vx,[1,2]);
VY  = repmat(vy,[1,2]);

FX  = repmat(Fx,[1,2]);
FY  = repmat(Fy,[1,2]);

%adjust particle sizes
R   = repmat(r ,[1,2]);
R(:,1) = sqrt( (phi+dphi1)/phi) * R(:,1);
R(:,2) = sqrt( (phi+dphi2)/phi) * R(:,2);

for k = 1:2
    x = X(:,k);
    y = Y(:,k);
    
    vx = VX(:,k);
    vy = VY(:,k);
    
    Fx = FX(:,k);
    Fy = FY(:,k);
    
    r = R(:,k);
    % relax new state using FIRE
    npPos       = 0;
    npNeg       = 0;
    npPMin      = 0;
    
    alpha       = alpha0;
    
    fireit      = 0;
    fcheck      = 10*Ftol;
    kcheck      = 10*Ktol;
    dt = dtstart;
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
            fprintf('\nOn State %.0f\n',k);
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
    
    X(:,k) = x;
    Y(:,k) = y;
    
    VX(:,k) = vx;
    VY(:,k) = vy;
    
    FX(:,k) = Fx;
    FY(:,k) = Fy;
end

%Volume fractions of states
phi1 = sum(pi*R(:,1).^2)/Lx/Ly;
phi2 = sum(pi*R(:,2).^2)/Lx/Ly;




%% Splicing - Oscillatory strain with discontinuous volume fraction change
% Copy the first state to a third column for the splicing stimulation
X   = [X,X(:,1)];
Y   = [Y,Y(:,1)];

VX  = [VX,VX(:,1)];
VY  = [VY,VY(:,1)];

R   = [R,R(:,1)];

FXold = [FX,FX(:,1)];
FYold = [FY,FY(:,1)];

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
NSTEPS   = ceil(totalT/dt);
NSKIP   = round(NSTEPS/NCYCLES/30);

% print to console
fprintf('\n\n\n\n &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n\n');
fprintf('           Now shearing jammed states\n\n');
fprintf(' &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n\n\n\n');

%Activation and deactivation timing
tActivation         = pi/2* (tCycle/2/pi);
tDeactivation       = 3*pi/2 * (tCycle/2/pi);

% store data during cyclic shear
instStrain      = zeros(NSTEPS,1);
potEnergy       = zeros(NSTEPS,3);
kinEnergy       = zeros(NSTEPS,3);
ShearStress     = zeros(NSTEPS,3);





%Animation set up
AnimateVideo = VideoWriter('SplicingAnimation.avi','Motion JPEG AVI');
open(AnimateVideo)

handle = figure(1);clf
set(handle,'Position',[0,0,1080,270],'units','points')

ax1 = subplot(1,4,1);
quiver1 = quiver(1,1,1,1,0,'LineWidth',3,'color','k','ShowArrowHead',0);
title(sprintf('phi = %2.3f' , phi1 ));
axis equal
box on; set(gca,'XTick',[],'YTick',[]);
axis([-Lx/4 5/4*Lx -Ly/4 5/4*Ly]);
ax1.Position = (ax1.Position - [0.10,0.13,0,0]).*[1,1,1.3,1.3];
scale = handle.Position(3)*ax1.Position(3)/(ax1.XLim(2) - ax1.XLim(1));

phi1BigR = animatedline;
set(phi1BigR    ,'Marker','o','MarkerEdgeColor','r','linestyle','none','MaximumNumPoints',9*Nlarge,'Markersize',2*R(end,1)*scale)
phi1SmallR = animatedline;
set(phi1SmallR  ,'Marker','o','MarkerEdgeColor','r','linestyle','none','MaximumNumPoints',9*Nsmall,'Markersize',2*R(1,1)*scale)
Boundary1 = animatedline;
set(Boundary1,'MaximumNumPoints',5,'color','m','linewidth',2)

ax2 = subplot(1,4,2);
quiver2 = quiver(1,1,1,1,0,'LineWidth',3,'color','k','ShowArrowHead',0);
title(sprintf('phi = %2.3f' , phi2 ));
axis equal
box on; set(gca,'XTick',[],'YTick',[]);
axis([-Lx/4 5/4*Lx -Ly/4 5/4*Ly]);
ax2.Position = (ax2.Position - [0.06,0.13,0,0]).*[1,1,1.3,1.3];
scale = handle.Position(3)*ax2.Position(3)/(ax2.XLim(2) - ax2.XLim(1));

phi2BigR = animatedline;
set(phi2BigR    ,'Marker','o','MarkerEdgeColor',[0 0.5 0],'linestyle','none','MaximumNumPoints',9*Nlarge,'Markersize',2*R(end,2)*scale)
phi2SmallR = animatedline;
set(phi2SmallR  ,'Marker','o','MarkerEdgeColor',[0 0.5 0],'linestyle','none','MaximumNumPoints',9*Nsmall,'Markersize',2*R(1,2)*scale)
Boundary2 = animatedline;
set(Boundary2,'MaximumNumPoints',5,'color','m','linewidth',2)

ax3 = subplot(1,4,3);
quiver3 = quiver(1,1,1,1,0,'LineWidth',3,'color','k','ShowArrowHead',0);
ax3title = title(sprintf('phi = %2.3f' , phi2 ));
axis equal
box on; set(gca,'XTick',[],'YTick',[]);
axis([-Lx/4 5/4*Lx -Ly/4 5/4*Ly]);
ax3.Position = (ax3.Position - [0.02,0.13,0,0]).*[1,1,1.3,1.3];
scale = handle.Position(3)*ax3.Position(3)/(ax3.XLim(2) - ax3.XLim(1));

SplicedBigR = animatedline;
set(SplicedBigR    ,'Marker','o','MarkerEdgeColor','b','linestyle','none','MaximumNumPoints',9*Nlarge)
SplicedSmallR = animatedline;
set(SplicedSmallR  ,'Marker','o','MarkerEdgeColor','b','linestyle','none','MaximumNumPoints',9*Nsmall)
Boundary3 = animatedline;
set(Boundary3,'MaximumNumPoints',5,'color','m','linewidth',2)

ax4 = subplot(1,4,4);
ax4title = title(sprintf('Loop # %d' , 0 ));
ylim(ax4,[-5e-6,5e-6])
xlim(ax4,[0,0.1])

ax4.Position = (ax4.Position - [-0.02,0,0,0]).*[1,1,1.4,1];

MaxNumPoints = 5*round(NSTEPS/NCYCLES);
Path1 = animatedline; set(Path1,'color','r','linewidth',2, 'MaximumNumPoints',MaxNumPoints);
Marker1 = animatedline; set(Marker1,'color','r','marker','o','MarkerFaceColor','r','Markersize',10,'MaximumNumPoints',1);
Path2 = animatedline; set(Path2,'color',[0,0.5,0],'linewidth',2, 'MaximumNumPoints',MaxNumPoints);
Marker2 = animatedline; set(Marker2,'color',[0,0.5,0],'marker','o','MarkerFaceColor',[0,0.5,0],'Markersize',10,'MaximumNumPoints',1);
Path3 = animatedline; set(Path3,'color','b','linewidth',2, 'MaximumNumPoints',MaxNumPoints);
Marker3 = animatedline; set(Marker3,'color','b','marker','o','MarkerFaceColor','b','Markersize',10,'MaximumNumPoints',1);


QUIVER = [quiver1,quiver2,quiver3];






% loop over oscillatory shear
while t < totalT
    %Particle radius of third column depends on time
    if mod(t-tActivation,tCycle) >= 0  && mod(t-tActivation,tCycle) < mod(tDeactivation-tActivation,tCycle)
        R(:,3) = R(:,1);
        phi3 = phi1;
    else
        R(:,3) = R(:,2);
        phi3 = phi2;
    end
    
    cijk = zeros(N,N,k);
    for k = 1:3
        x = X(:,k);
        y = Y(:,k);
        
        vx = VX(:,k);
        vy = VY(:,k);
        
        r = R(:,k);
        Fxold = FXold(:,k);
        Fyold = FYold(:,k);
        
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
                        cijk(ii,jj,k) = 1;
                        cijk(jj,ii,k) = 1;
                        
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
        potEnergy(ss,k) = U;
        kinEnergy(ss,k) = 0.5*sum(vx.^2 + vy.^2);
        
        % store shear stress
        ShearStress(ss,k) = vstress(2)/Lx/Ly;
        
        
        X(:,k) = x;
        Y(:,k) = y;
        
        VX(:,k) = vx;
        VY(:,k) = vy;
        
        FXold(:,k) = Fxold;
        FYold(:,k) = Fyold;
    end
    
    % save strain
    instStrain(ss) = currStrain;
    
    % output
    if mod(ss,NSKIP) == 0
        % print to console
        fprintf('t = %0.3g/%0.3g, ss = %d, gam = %0.3g\n',t,totalT,ss,currStrain);
    end
    
    %Animate
    if ~mod(ss,NSKIP)
        for k = 1:3
            q = 1;
            Xquiver = zeros(N^2,1);
            Yquiver = zeros(N^2,1);
            Uquiver = zeros(N^2,1);
            Vquiver = zeros(N^2,1);
            for xx = -1:1
                for yy = -1:1
                    %Animate state 1
                    addpoints(phi1BigR  ,X(Nsmall+1:end,1)+xx*Lx + currStrain*yy*Ly, Y(Nsmall+1:end,1)+yy*Ly);
                    addpoints(phi1SmallR,X(1:Nsmall,1)    +xx*Lx + currStrain*yy*Ly, Y(1:Nsmall,1)    +yy*Ly);
                    
                    %Animate state 2
                    addpoints(phi2BigR  ,X(Nsmall+1:end,2)+xx*Lx + currStrain*yy*Ly, Y(Nsmall+1:end,2)+yy*Ly);
                    addpoints(phi2SmallR,X(1:Nsmall,2)    +xx*Lx + currStrain*yy*Ly, Y(1:Nsmall,2)    +yy*Ly);
                    
                    %Animate Spliced
                    set(ax3title,'String',sprintf('phi = %2.3f' , phi3 ))
                    set(SplicedBigR  ,'MarkerSize',2*R(end,3)*scale)
                    set(SplicedSmallR,'MarkerSize',2*R(1  ,3)*scale)
                    addpoints(SplicedBigR  ,X(Nsmall+1:end,3)+xx*Lx + currStrain*yy*Ly, Y(Nsmall+1:end,3)+yy*Ly);
                    addpoints(SplicedSmallR,X(1:Nsmall,3)    +xx*Lx + currStrain*yy*Ly, Y(1:Nsmall,3)    +yy*Ly);
                    
                    %draw contacts
                    for nn = 1:N
                        for mm = nn+1:N
                            if cijk(nn,mm,k) == 1
                                dy = Y(mm,k) - Y(nn,k);
                                im = round(dy/Ly);
                                dy = dy - Ly*im;
                                
                                dx = X(mm,k) - X(nn,k);
                                dx = dx - Ly*im*currStrain;
                                dx = dx - Lx*round(dx/Lx);
                       
                                Xquiver(q) = X(nn,k) + xx*Lx + currStrain*yy*Ly;
                                Yquiver(q) = Y(nn,k) + yy*Ly;
                                Uquiver(q) = dx;
                                Vquiver(q) = dy;
                                q = q+1;
                            end
                        end
                    end
                end
            end
            QUIVER(k).XData = Xquiver;
            QUIVER(k).YData = Yquiver;
            QUIVER(k).UData = Uquiver;
            QUIVER(k).VData = Vquiver;
        end
        addpoints(Boundary1,[0 1 1+currStrain currStrain 0 1]*Lx,[0 0 1 1 0 0]*Ly)
        addpoints(Boundary2,[0 1 1+currStrain currStrain 0 1]*Lx,[0 0 1 1 0 0]*Ly)
        addpoints(Boundary3,[0 1 1+currStrain currStrain 0 1]*Lx,[0 0 1 1 0 0]*Ly)
        
        %Work loop
        set(ax4title,'String',sprintf('Loop # %d' , floor(ss*NCYCLES/NSTEPS) +1));
        addpoints(Path1,instStrain(ss-NSKIP+1:ss),ShearStress(ss-NSKIP+1:ss,1))
        addpoints(Marker1,currStrain,ShearStress(ss,1));
        addpoints(Path2,instStrain(ss-NSKIP+1:ss),ShearStress(ss-NSKIP+1:ss,2))
        addpoints(Marker2,currStrain,ShearStress(ss,2));
        addpoints(Path3,instStrain(ss-NSKIP+1:ss),ShearStress(ss-NSKIP+1:ss,3))
        addpoints(Marker3,currStrain,ShearStress(ss,3));
        
        %Update figure
        drawnow limitrate
        writeVideo(AnimateVideo,getframe(gcf))
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
    X = X + dStrain.*Y;
end
close(AnimateVideo); clear AnimateVideo