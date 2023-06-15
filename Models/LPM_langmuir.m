function [sol, t, x, C] = LPM_langmuir(feedProf, parameter, L, Di, epsb, epsp, Rp, Q, Cfeed, keff, Dax, tpulse, tfinal, opt)
% Lumped Pore Model (LPM) considering mass transfer resistence and axial dispersion.
% Langmuir isotherm
% Change feed profile between pulse (e.g.: chromatografic peak) and step (e.g.: breakthrough experiment)
% Multicomponent
% Uses pdepe function to solve the system of partial differential equations


global data res

%% Default arguments (Example)
if nargin == 0
    feedProf =  'step';  % feed profile, can be 'pulse' (e.g.: chromatografic peak) or 'step' (e.g.: breakthrough experiment)
    parameter = [2.260 5.35e-3; 2.579 4.21e-3];  % isotherm parameters (depends on the isotherm model chosen)
    L = 15;  % cm, column length
    Di = 1;  % cm, column internal diameter
    epsb = 0.335;  % column bulk porosity
    epsp = 0.475;  % particle porosity
    Rp = 10e-4;  % cm, particle radius
    Q = 1;  % mL/min, flow rate
    Cfeed = [4.35 ; 5.25];  % g/L, feed concentration
    keff = [0.339 ; 0.348];  % cm/min, effective mass transfer coefficient
    Dax = [1.873e-3 ; 1.847e-3];  % cm2/min, axial dispersion coefficient
    tpulse = 1;  % min, feed pulse duration. For a step injection set tpulse = tfinal
    tfinal = 30;  % min, final time for calculation
    opt.npz = 400;  % number of discretization points in z
    opt.npt = 400;  % number of discretization points in t
    opt.fig = true;  % true - show figures; false - do not show figures
end


%% Calculations
addpath('../')
data = struct('feedProf',feedProf,'parameter',parameter','Cfeed',Cfeed,'npz',opt.npz,'npt',opt.npt,'Q',Q,'Rp',Rp,'epsb',epsb,'epsp',epsp,'tpulse',tpulse,'keff',keff,'Dax',Dax);
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

c_mbf = []; c_mbp = [];
f_mbf = []; f_mbp = [];
s_mbf = []; s_mbp = [];

for i = 1:data.nc
    c_mbf = [ c_mbf ; 1 ];
    c_mbp = [ c_mbp ; data.epsb + (1-data.epsp) * data.parameter(i,1) / (1+data.parameter(i,2)*u(data.nc+i))^2 ];
    f_mbf = [ f_mbf ; data.Dax(i) ];
    f_mbp = [ f_mbp ; 0 ];
    s_mbf = [ s_mbf ; -data.ui * DuDx(i) - (1-data.epsb) / data.epsb * data.keff(i) * 3/data.Rp * (u(i) - u(data.nc+i)) ];
    s_mbp = [ s_mbp ; data.keff(i) * 3/data.Rp * (u(i) - u(data.nc+i)) ];
end

c = [ c_mbf ; c_mbp];
f = [ f_mbf ; f_mbp].*DuDx;
s = [ s_mbf ; s_mbp];


function u0 = pde_ic(x)
% pdepe function describing initial conditions
global data

u0_mbf = []; u0_mbp = [];
for i = 1:data.nc 
    u0_mbf = [ u0_mbf ; 0 ];
    u0_mbp = [ u0_mbp ; 0 ];
end
u0 = [ u0_mbf ; u0_mbp];


function [pl,ql,pr,qr] = pde_bc(xl,ul,xr,ur,t)
% pdepe function describing boundary conditions
global data

Cin = setFeedProfile(data.feedProf, t, data.tpulse, data.Cfeed);

% Left boundary conditions (z = 0)
pl_mbf = []; pl_mbp = [];
ql_mbf = []; ql_mbp = [];
for i = 1:data.nc
    pl_mbf = [ pl_mbf ; -data.ui * (ul(i) - Cin(i)) ];
    pl_mbp = [ pl_mbp ; 0 ];
    ql_mbf = [ ql_mbf ; 1 ];
    ql_mbp = [ ql_mbp ; 3.14 ];
end
pl = [ pl_mbf ; pl_mbp ]; 
ql = [ ql_mbf ; ql_mbp ];

% Right boundary conditions (z = L)
pr_mbf = []; pr_mbp = [];
qr_mbf = []; qr_mbp = [];
for i = 1:data.nc
    pr_mbf = [ pr_mbf ; 0 ];
    pr_mbp = [ pr_mbp ; 0 ];
    qr_mbf = [ qr_mbf ; 1/data.Dax(i) ];
    qr_mbp = [ qr_mbp ; 3.14 ];
end
pr = [ pr_mbf ; pr_mbp ]; 
qr = [ qr_mbf ; qr_mbp ];


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
