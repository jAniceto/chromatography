function [sol, t, x] = transport_dispersive_ldf_2c(isoType,feedProf,parameter,L,Di,epsb,Q,Cfeed,KLDF,Dax,tpulse,tfinal,opt)
% Transport-Dispersive Model (TDM) considering mass transfer resistence in the solid to be dominant and 
% using the Linear Driving Force Model (LDF) approach (Glauckauf and Coates, 1947)
% Change the isotherm model by modifying the isotherm function
% Change feed profile between pulse (e.g.: chromatografic peak) and step (e.g.: breakthrough experiment)
% Two components
% Uses pdepe functin to solve the system of partial differential equations


global data res

%% Default arguments (Example)
if nargin == 0
    isoType =   'linear-langmuir';      % isotherm type, can be 'linear' or 'linear-langmuir'
    feedProf =  'pulse';                % feed profile, can be 'pulse' (e.g.: chromatografic peak) or 'step' (e.g.: breakthrough experiment)
    parameter = [5.5*0.13 0.13 1.99 ; 6.5*0.13 0.13 2.50]';   % isotherm parameters (depends on the isotherm model chosen)
    L =         10;                     % cm, column length
    Di =        1;                      % cm, column internal diameter
    epsb =      0.708;                  % column bulk porosity
    Q =         4;                      % mL/min, flow rate
    Cfeed = 	[0.6*300/50 ; 0.6*300/50];             % g/L, feed concentration
    KLDF =      [13.3*60/10 ; 8*60/10];             % min-1, linear driving force (LDF) mass transfer coefficient
    Dax =       [5.57e-3 ; 6.57e-3];              % cm2/min, axial dispersion coefficient
    tpulse =    50*0.001/4;             % min, feed pulse duration. For a step injection set tpulse = tfinal
    tfinal =    7;                      % min, final time for calculation
    opt.npz =   150;                    % number of discretization points in z
    opt.npt =   100;                    % number of discretization points in t
    opt.fig =   true;                   % true - show figures; false - do not show figures
end


%% Calculations
addpath('../')
data = struct('isoType',isoType,'feedProf',feedProf,'parameter',parameter,'Cfeed',Cfeed,'npz',opt.npz,'npt',opt.npt,'Q',Q,'epsb',epsb,'tpulse',tpulse,'KLDF',KLDF,'Dax',Dax);
nc = length(Cfeed);
data.nc=nc;
A = pi()*Di^2/4; % cm2
data.F = (1-epsb)/epsb;
data.ui = Q/epsb/A; % cm/min

% Run pdpe
m = 0;
x = linspace(0,L,opt.npz);
t = linspace(0,tfinal,opt.npt);
sol = pdepe(m,@pde,@pde_ic,@pde_bc,x,t);
u1 = sol(:,:,1);
u2 = sol(:,:,2);

C1 = u1'; C2 = u2';
res.t = t'; res.z = x'; res.C = C1;


% Plot figures
if isfield(opt,'fig') && (opt.fig == 1)
    % Concentration history at column exit (Chromatogram)
    h1 = figure;
    plot(t,C1(end,:), t,C2(end,:));
    axis([0 tfinal  0 inf]) % fix the axes
    xlabel('t')
    ylabel('C')
    
%     figure
%     surf(x,t,u1)
%     title('Concentration C(x,t)')
%     xlabel('Distance z')
%     ylabel('Time t')
%     
%     figure
%     surf(x,t,u2)
%     title('Solid loading q(x,t)')
%     xlabel('Distance z')
%     ylabel('Time t')
% 
%     % Concentration over time with slider figure    
%     slidermin = 1; % t = 0
%     slidermax = opt.npt; % t = tfinal
%     res.hslider = plot(x',C1(:,end));
% 	res.ymax = str2double(num2str(max(C1(:,2)),1));
%     axis([0 data.L  0 res.ymax])
%     xlabel('Position, z')
%     ylabel('Concentration, C')
%     anno_text = sprintf('t = %.1f', 0);
%     res.time_box = annotation('textbox',[.15 .48 .4 .5],'String',anno_text,'EdgeColor', 'none','FitBoxToText','on');
%     sld = uicontrol('Style', 'slider',...
%         'Min',slidermin,'Max',slidermax,'Value',slidermin,...
%         'SliderStep',[1 1]./(slidermax-slidermin),...
%         'Position', [150 390 120 20],...
%         'Callback', @sliderplot); 

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

c = [1 ; 1 ; 1 ; 1]; 

f = [data.Dax(1) ; data.Dax(2) ; 0 ; 0].* DuDx; 

qast = isotherm(data.isoType, u(1), u(2), data.parameter);

s = [-data.ui*DuDx(1)-(1-data.epsb)/data.epsb*data.KLDF(1)*(qast(1)-u(3)) ;
     -data.ui*DuDx(2)-(1-data.epsb)/data.epsb*data.KLDF(2)*(qast(2)-u(4)) ;
     data.KLDF(1)*(qast(1)-u(3)) ;
     data.KLDF(2)*(qast(2)-u(4)) ]; 


function u0 = pde_ic(x)
% pdepe function describing initial conditions

u0 = [0 ; 0 ; 0 ; 0]; 


function [pl,ql,pr,qr] = pde_bc(xl,ul,xr,ur,t)
% pdepe function describing boundary conditions
global data

Cin = setFeedProfile(data.feedProf, t, data.tpulse, data.Cfeed);

% Left boundary conditions (z = 0)
pl = [ Cin(1)-ul(1) ; Cin(2)-ul(2) ; 0 ; 0 ]; 
ql = [ 1/data.ui ; 1/data.ui ; 3.14 ; 3.14 ]; 

% Right boundary conditions (z = L)
pr = [0 ; 0 ; 0 ; 0]; 
qr = [ 1/data.Dax(1) ; 1/data.Dax(2) ; 3.14 ; 3.14 ];


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


%% Isotherm
function q=isotherm(isoType, c1, c2, parameter)
% Defines the isotherm
%       Linear isotherm: q1 = H1*C1  , parameter = H
%       Linear-Langmuir isotherm:   q1(C1,C2) = m1*C1 + a1*C1/(1+b1*C1) , parameter = [a1 b1 m1]


if strcmp(isoType,'linear') 
    q = parameter*c;
    
elseif strcmp(isoType,'langmuir')
    q = parameter(1)*c/(1+parameter(2)*c);
    
elseif strcmp(isoType,'linear-langmuir')
%         q = parameter(3)*c + parameter(1)*c/(1+parameter(2)*c);
        q1 = parameter(3,1)*c1 + parameter(1,1)*c1/(1+parameter(2,1)*c1+parameter(2,2)*c2); 
        q2 = parameter(3,2)*c2 + parameter(1,2)*c2/(1+parameter(2,1)*c1+parameter(2,2)*c2); 
    
else
    error('Invalid isotherm type. isoType must be "linear" or "langmuir" or "linear-langmuir"')
end

q = [q1 ; q2];
