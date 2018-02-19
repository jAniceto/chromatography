function sol=EDMlinear_pdepe_1c(feedProf,H,L,Di,epsb,epsp,Q,Cfeed,Dapp,tpulse,tfinal,opt)
% Equilibrium-Dispersive Model (EDM) 
% Linear isotherm model
% Change feed profile between pulse (e.g.: chromatografic peak) and step (e.g.: breakthrough experiment)
% Single component
% Uses pdepe functin to solve the system of partial differential equations


global data res
clc

%% Default arguments (Example)
if nargin == 0
    feedProf =  'pulse';                % feed profile, can be 'pulse' (e.g.: chromatografic peak) or 'step' (e.g.: breakthrough experiment)
    H =         5.5;                    % linear isotherm constant
    L =         10;                     % cm, column length
    Di =        1;                      % cm, column internal diameter
    epsb =      0.4;                    % column bulk porosity
    epsp =      0.5;                    % particle porosity
    Q =         4;                      % mL/min, flow rate
    Cfeed = 	0.6*300/50;             % g/L, feed concentration
    KLDF =      13.3*60/10;             % min-1, linear driving force (LDF) mass transfer coefficient
    Dapp =       5.57e-3;                % cm2/min, axial dispersion coefficient
    tpulse =    50*0.001/4;             % min, feed pulse duration. For a step injection set tpulse = tfinal
    tfinal =    7;                      % min, final time for calculation
    opt.npz =   150;                    % number of discretization points in z
    opt.npt =   100;                    % number of discretization points in t
    opt.fig =   true;                   % true - show figures; false - do not show figures
end


%% Calculations
data = struct('feedProf',feedProf,'H',H,'Cfeed',Cfeed,'npz',opt.npz,'npt',opt.npt,'Q',Q,'epsb',epsb,'epsp',epsp,'tpulse',tpulse,'Dapp',Dapp);
nc = length(Cfeed);
data.nc=nc;
A = pi()*Di^2/4; % cm2
data.F = (1-epsb)/epsb;
data.ui = Q/epsb/A; % cm/min

% Run pdpe
m = 0;
x = linspace(0,L,opt.npz);
t = linspace(0,tfinal,opt.npt);
options = odeset('MaxOrder',2);
sol = pdepe_mod(m,@pde,@pde_ic,@pde_bc,x,t,options);
u1 = sol(:,:,1);

C1 = u1';
res.t = t'; res.z = x'; res.C = C1;


% Plot figures
if isfield(opt,'fig') && (opt.fig == 1)
    % Concentration history at column exit (Chromatogram)
    h1 = figure;
    plot(t,C1(end,:));
    axis([0 tfinal  0 inf]) % fix the axes
    xlabel('t')
    ylabel('C')
    
    figure
    surf(x,t,u1)
    title('Concentration C(x,t)')
    xlabel('Distance z')
    ylabel('Time t')
    
    figure
    surf(x,t,u2)
    title('Solid loading q(x,t)')
    xlabel('Distance z')
    ylabel('Time t')

    % Concentration over time with slider figure    
    slidermin = 1; % t = 0
    slidermax = opt.npt; % t = tfinal
    res.hslider = plot(x',C1(:,end));
	res.ymax = str2double(num2str(max(C1(:,2)),1));
    axis([0 tfinal  0 res.ymax])
    xlabel('Position, z')
    ylabel('Concentration, C')
    anno_text = sprintf('t = %.1f', 0);
    res.time_box = annotation('textbox',[.15 .48 .4 .5],'String',anno_text,'EdgeColor', 'none','FitBoxToText','on');
    sld = uicontrol('Style', 'slider',...
        'Min',slidermin,'Max',slidermax,'Value',slidermin,...
        'SliderStep',[1 1]./(slidermax-slidermin),...
        'Position', [150 390 120 20],...
        'Callback', @sliderplot); 

end


%% Auxiliary functions
function sliderplot(source,event)
% Handles the UI slider in the Concentration vs Position plot
global res

val = round(get(source,'Value'));
set(res.hslider,'YData',res.C(:,val));
anno_text = sprintf('t = %.2f',res.t(val));
set(res.time_box,'String',anno_text);


%% pdepe functions
function [c,f,s] = pde(x,t,u,DuDx)
% Main pdepe function describing the system of partial diferential equations
global data

c = 1 + (1-data.epsb)*(data.epsp+data.H)/data.epsb ; 

f =  data.Dapp*DuDx; 

s = -data.ui*DuDx ; 


function u0 = pde_ic(x)
% pdepe function describing initial conditions

u0 = 0 ; 


function [pl,ql,pr,qr] = pde_bc(xl,ul,xr,ur,t)
% pdepe function describing boundary conditions
global data

Cin = setFeedProfile(data.feedProf, t, data.tpulse, data.Cfeed);

% Left boundary conditions (z = 0)
pl =  Cin-ul(1) ; 
ql =  0 ; 

% Right boundary conditions (z = L)
pr =  0 ;  
qr =  1/data.Dapp ; 


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
