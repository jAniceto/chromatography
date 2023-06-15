function [sol, t, x, C] = transport_dispersive_ldf(isoType, feedProf, parameter, L, Di, epsb, Q, Cfeed, KLDF, Dax, tpulse, tfinal, opt)
% Transport-Dispersive Model (TDM) considering mass transfer resistence in the solid to be dominant and 
% using the Linear Driving Force Model (LDF) approach (Glauckauf and Coates, 1947)
% Change the isotherm model by modifying the isotherm function
% Change feed profile between pulse (e.g.: chromatografic peak) and step (e.g.: breakthrough experiment)
% Multicomponent
% Uses pdepe functin to solve the system of partial differential equations


global data res

%% Default arguments (Example)
if nargin == 0
    isoType =   'linear-langmuir';                                                  % isotherm type, can be 'linear' or 'linear-langmuir'
    feedProf =  'pulse';                                                            % feed profile, can be 'pulse' (e.g.: chromatografic peak) or 'step' (e.g.: breakthrough experiment)
    parameter = [5.5*0.13 0.13 1.99 ; 6.5*0.13 0.13 2.50  ; 6.8*0.13 0.13 2.80];    % isotherm parameters (depends on the isotherm model chosen)
    L =         10;                                                                 % cm, column length
    Di =        1;                                                                  % cm, column internal diameter
    epsb =      0.708;                                                              % column bulk porosity
    Q =         4;                                                                  % mL/min, flow rate
    Cfeed = 	[0.6*300/50 ; 0.6*300/50 ; 0.6*300/50];                             % g/L, feed concentration
    KLDF =      [13.3*60/10 ; 8*60/10  ; 6*60/10];                                  % min-1, linear driving force (LDF) mass transfer coefficient
    Dax =       [5.57e-3 ; 6.57e-3  ; 7.57e-3];                                     % cm2/min, axial dispersion coefficient
    tpulse =    50*0.001/4;                                                         % min, feed pulse duration. For a step injection set tpulse = tfinal
    tfinal =    7;                                                                  % min, final time for calculation
    opt.npz =   150;                                                                % number of discretization points in z
    opt.npt =   100;                                                                % number of discretization points in t
    opt.fig =   true;                                                               % true - show figures; false - do not show figures
end


%% Calculations
addpath('../')
data = struct('isoType',isoType,'feedProf',feedProf,'parameter',parameter','Cfeed',Cfeed,'npz',opt.npz,'npt',opt.npt,'Q',Q,'epsb',epsb,'tpulse',tpulse,'KLDF',KLDF,'Dax',Dax);
nc = length(Cfeed);
data.nc = nc;
A = pi()*Di^2/4; % cm2
data.F = (1-epsb)/epsb;
data.ui = Q/epsb/A; % cm/min

% Run pdepe
m = 0;
x = linspace(0,L,opt.npz);
t = linspace(0,tfinal,opt.npt);
sol = pdepe(m,@pde,@pde_ic,@pde_bc,x,t);

for i = 1:data.nc
    C{i} = sol(:,:,i);    
end
res.t = t'; res.z = x'; res.C = C;


% Plot figures
if isfield(opt,'fig') && (opt.fig == 1)
    % Concentration history at column exit (Chromatogram)
    h1 = figure;
    hold all
    for i = 1:data.nc
        plot(t,C{i}(:,end));
    end
    hold off;
    axis([0 tfinal  0 inf]) % fix the axes
    xlabel('t')
    ylabel('C')
    
    for i = 1:data.nc
        figure
        surf(x,t,C{i})
        title(sprintf('Concentration C%i(x,t)',i))
        xlabel('Distance z')
        ylabel('Time t')
    end

    % Concentration over time with slider figure 
    figure
    hold all
    slidermin = 1; % t = 0
    slidermax = opt.npt; % t = tfinal
    ymax_aux = [];
    for i = 1:data.nc
        res.hslider(i) = plot(x',C{i}(end,:));
        ymax_aux = [ymax_aux max(C{i}(2,:))];
    end
    hold off
	res.ymax = str2double(num2str(max(ymax_aux),1));
    axis([0 L  0 res.ymax])
    xlabel('Position, z')
    ylabel('Concentration, C')
    anno_text = sprintf('t = %.1f', 0);
    res.time_box = annotation('textbox',[.55 .48 .4 .5],'String',anno_text,'EdgeColor', 'none','FitBoxToText','on');
    sld = uicontrol('Style', 'slider',...
        'Min',slidermin,'Max',slidermax,'Value',slidermin,...
        'SliderStep',[1 1]./(slidermax-slidermin),...
        'Position', [390 390 120 20],...
        'Callback', @sliderplot); 
end


%% Auxiliary functions
function sliderplot(source,event)
% Handles the UI slider in the Concentration vs Position plot
global data res

val = round(get(source,'Value'));
for i = 1:data.nc
    set(res.hslider(i),'YData',res.C{i}(val,:));
end
anno_text = sprintf('t = %.2f',res.t(val));
set(res.time_box,'String',anno_text);


%% pdepe functions
function [c,f,s] = pde(x,t,u,DuDx)
% Main pdepe function describing the system of partial diferential equations
global data

qast = isotherm(data.isoType, u(1:data.nc), data.parameter);

c_mbf = []; c_mbs = [];
f_mbf = []; f_mbs = [];
s_mbf = []; s_mbs = [];

for i = 1:data.nc
    c_mbf = [ c_mbf ; 1 ];
    c_mbs = [ c_mbs ; 1 ];
    f_mbf = [ f_mbf ; data.Dax(i) ];
    f_mbs = [ f_mbs ; 0 ];
    s_mbf = [ s_mbf ; -data.ui*DuDx(i)-(1-data.epsb)/data.epsb*data.KLDF(i)*(qast(i)-u(data.nc+i)) ];
    s_mbs = [ s_mbs ; data.KLDF(i)*(qast(i)-u(data.nc+i)) ];
end

c = [ c_mbf ; c_mbs];
f = [ f_mbf ; f_mbs].*DuDx;
s = [ s_mbf ; s_mbs];


function u0 = pde_ic(x)
% pdepe function describing initial conditions
global data

u0_mbf = []; u0_mbs = [];
for i = 1:data.nc 
    u0_mbf = [ u0_mbf ; 0 ];
    u0_mbs = [ u0_mbs ; 0 ];
end
u0 = [ u0_mbf ; u0_mbs];


function [pl,ql,pr,qr] = pde_bc(xl,ul,xr,ur,t)
% pdepe function describing boundary conditions
global data

Cin = setFeedProfile(data.feedProf, t, data.tpulse, data.Cfeed);

% Left boundary conditions (z = 0)
pl_mbf = []; pl_mbs = [];
ql_mbf = []; ql_mbs = [];
for i = 1:data.nc
    pl_mbf = [ pl_mbf ; Cin(i)-ul(i) ];
    pl_mbs = [ pl_mbs ; 0 ];
    ql_mbf = [ ql_mbf ; 1/data.ui ];
    ql_mbs = [ ql_mbs ; 3.14 ];
end
pl = [ pl_mbf ; pl_mbs ]; 
ql = [ ql_mbf ; ql_mbs ];

% Right boundary conditions (z = L)
pr_mbf = []; pr_mbs = [];
qr_mbf = []; qr_mbs = [];
for i = 1:data.nc
    pr_mbf = [ pr_mbf ; 0 ];
    pr_mbs = [ pr_mbs ; 0 ];
    qr_mbf = [ qr_mbf ; 1/data.Dax(i) ];
    qr_mbs = [ qr_mbs ; 3.14 ];
end
pr = [ pr_mbf ; pr_mbs ]; 
qr = [ qr_mbf ; qr_mbs ];


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
