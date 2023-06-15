function sol = plug_flow(feedProf, L, Di, epsb, Q, Cfeed, Dax, tpulse, tfinal, opt)
% Plug flow model - no adsorption, no mass transfer
% Change feed profile between pulse (e.g.: chromatografic peak) and step (e.g.: breakthrough experiment)
% Single component
% Uses finite diferences to solve the system of partial differential equations


%% Default arguments (Example)
if nargin == 0
    feedProf =  'pulse';                % feed profile, can be 'pulse' (e.g.: chromatografic peak) or 'step' (e.g.: breakthrough experiment)
    L =         10;                     % cm, column length
    Di =        1;                      % cm, column internal diameter
    epsb =      0.708;                  % column bulk porosity
    Q =         4;                      % mL/min, flow rate
    Cfeed = 	0.6*300/50;             % g/L, feed concentration
    Dax =       5.57e-3;                % cm2/min, axial dispersion coefficient
    tpulse =    50*0.001/4;             % min, feed pulse duration. For a step injection set tpulse = tfinal
    tfinal =    3;                      % min, final time for calculation
    opt.npz =   2000;                   % number of discretization points in z
    opt.npt =   2000;                   % number of discretization points in t
    opt.fig =   true;                   % true - show figures; false - do not show figures
end


%% Calculations
addpath('../')
data = struct('feedProf',feedProf,'Cfeed',Cfeed,'npz',opt.npz,'npt',opt.npt,'Q',Q,'epsb',epsb,'tpulse',tpulse,'Dax',Dax);
nc = length(Cfeed); data.nc = nc;  % number of components
A = pi()*Di^2/4;  % cm2
data.F = (1-epsb)/epsb;
data.ui = Q/epsb/A;  % cm/min

tspan = 0:tfinal/(data.npt-1):tfinal; % time span (min)
y0 = zeros(nc*data.npz,1);
data.h = L/(data.npz-1);

sol = ode45(@sedo, tspan, y0, [], data);

t = sol.x;
for i = 1:nc
    C(:,:,i) = sol.y(data.npz*i-data.npz+1:data.npz*i,:);
end


% Plot figures
if isfield(opt,'fig') && (opt.fig == 1)
    % Concentration history at column exit (Chromatogram)
    figure;
    for i = 1:nc
        plot(t, C(end,:,i), 'LineWidth',1.5);
        hold on;
    end
    hold off;
    axis([0 tfinal  0 inf]) % fix the axes
    xlabel('t')
    ylabel('C')
end


%% Solving PDE using pdepe matlab function
function DyDt = sedo(t, y, data)

N=data.npz;
nc=data.nc;

Cinj = setFeedProfile(data.feedProf, t, data.tpulse, data.Cfeed);

DyDt = zeros(nc*N,1);

for j=1:nc
    
    y(N*j-N+1) = (2*data.h*data.ui*Cinj(j) - data.Dax(j)*(-y(N*j-N+3)+4*y(N*j-N+2)))/(2*data.h*data.ui-3*data.Dax(j)); % <<< From the boundary condition: z = 0 , C = Cinj + Dax/u * dC/dz
    y(N*j) = 4/3*y(N*j-1)-1/3*y(N*j-2); % <<< From the boundary condition: z = L , dC/dz=0
      
    % Forward finite differences
    DyDt(N*j-N+1) = data.Dax(j) * 1/(data.h^2)*(-y(N*j-N+4)+4*y(N*j-N+3)-5*y(N*j-N+2)+2*y(N*j-N+1)) - data.ui * 1/(2*data.h)*(-y(N*j-N+3)+4*y(N*j-N+2)-3*y(N*j-N+1));
        
    % Central finite differences
    for i=2:N-1
        DyDt(N*j-N+i) = data.Dax(j) * 1/(data.h^2)*(y(N*j-N+i+1)-2*y(N*j-N+i)+y(N*j-N+i-1)) - data.ui * 1/(2*data.h)*(y(N*j-N+i+1)-y(N*j-N+i-1));
    end
    
    % Backward finite differences
    DyDt(N*j) = 4/3*DyDt(N*j-1) - 1/3*DyDt(N*j-2); % <<< From the boundary condition: z = L , dC/dz=0;
end


%% Feed profile 
function Cinj=setFeedProfile(feedProf, t, tpulse, Cfeed)
% Defines the feed profile. Can be 'pulse' (e.g.: chromatografic peak) or 'step' (e.g.: breakthrough experiment)

if strcmp(feedProf,'pulse')
    if t<=tpulse
        Cinj = Cfeed;
    else
        Cinj = zeros(1,length(Cfeed));
    end
    
elseif strcmp(feedProf,'step')
    Cinj = Cfeed;
    
else
    error('Invalid feed profile. feedProf must be "step" or "pulse"')
end