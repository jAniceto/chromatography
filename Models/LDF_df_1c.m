function sol=LDF_df_1c(isotherm,feedProf,parameter,L,Di,epsb,Q,Cfeed,KLDF,Dax,tpulse,tfinal,opt)
% Transport-Dispersive Model (TDM) considering mass transfer resistence in the solid to be dominant and 
% using the Linear Driving Force Model (LDF) approach (Glauckauf and Coates, 1947)
% Choose the isotherm model: linear, linearlangmuir, langmuir
% Change feed profile between pulse (e.g.: chromatografic peak) and step (e.g.: breakthrough experiment)
% Single component
% Uses finite diferences to solve the system of partial differential equations


%% Default arguments (Example)
if nargin == 0
    isotherm =  'linearlangmuir';      % isotherm model
    feedProf =  'pulse';                % feed profile, can be 'pulse' (e.g.: chromatografic peak) or 'step' (e.g.: breakthrough experiment)
    parameter = [5.5*0.13 0.13 1.99];   % isotherm parameters (depends on the isotherm model chosen)
    L =         10;                     % cm, column length
    Di =        1;                      % cm, column internal diameter
    epsb =      0.708;                  % column bulk porosity
    Q =         4;                      % mL/min, flow rate
    Cfeed = 	0.6*300/50;             % g/L, feed concentration
    KLDF =      13.3*60/10;             % min-1, linear driving force (LDF) mass transfer coefficient
    Dax =       5.57e-3;                % cm2/min, axial dispersion coefficient
    tpulse =    50*0.001/4;             % min, feed pulse duration. For a step injection set tpulse = tfinal
    tfinal =    7;                      % min, final time for calculation
    opt.npz =   150;                    % number of discretization points in z
    opt.npt =   100;                    % number of discretization points in t
    opt.fig =   true;                   % true - show figures; false - do not show figures
end


%% Calculations
addpath('../')
data = struct('feedProf',feedProf,'parameter',parameter,'Cfeed',Cfeed,'npz',opt.npz,'npt',opt.npt,'Q',Q,'epsb',epsb,'tpulse',tpulse,'KLDF',KLDF,'Dax',Dax);
data.isoT = str2func(isotherm);
nc = length(Cfeed);
data.nc=nc;
A = pi()*Di^2/4; % cm2
data.F = (1-epsb)/epsb;
data.ui = Q/epsb/A; % cm/min

tspan = 0:tfinal/(data.npt-1):tfinal; % time span (min)
y0 = zeros(2*nc*data.npz,1);
data.h = L/(data.npz-1);

sol = ode45(@sedo,tspan,y0,[],data);
t = sol.x;
C1 = sol.y(1:data.npz,:);
q1 = sol.y((data.npz+1):end,:);


% Plot figures
if isfield(opt,'fig') && (opt.fig == 1)
    % Concentration history at column exit (Chromatogram)
    h1 = figure;
    plot(t,C1(end,:));
    axis([0 tfinal  0 inf]) % fix the axes
    xlabel('t')
    ylabel('C')
    
%     % Concentration profile inside the column through time
%     h2 = figure;
%     z = 0:data.h:L; z = z';
%     xlabel('z')
%     ylabel('C')
%     axis([0 z(end)  0 ceil(Cfeed)]) % fix the axes
%     set(gca,'NextPlot','replacechildren')
%     
%     for i=1:length(t) 
%         plot(z,C1(:,i));
%         anno_text = sprintf('t = %.1f', t(i));
%         if exist('time_box', 'var') == 1 
%             delete(time_box)
%         end
%         time_box = annotation('textbox',[.15 .48 .4 .5],'String',anno_text,'EdgeColor', 'none','FitBoxToText','on');
% %         drawnow;
%         movieframes(i)=getframe;
%     end
%     
%     h3 = figure;
%     movie(movieframes)
end




%% Solving PDE using pdepe matlab function
function DyDt = sedo(t,y,data)

N=data.npz;
nc=data.nc;

% CC = zeros(N,nc);
% for j=1:nc
%     CC(:,j)=y(N*j-N+1:N*j-N+N);
% end
% OR %
CC = reshape(y(1:nc*N),N,nc);
qast = data.isoT(CC,data.parameter);
qast = reshape(qast,N*nc,1);

Cinj = setFeedProfile(data.feedProf,t,data.tpulse,data.Cfeed);


DyDt = zeros(2*nc*N,1);

for j=1:nc
    
    y(N*j-N+1) = (2*data.h*data.ui*Cinj(j) - data.Dax(j)*(-y(N*j-N+3)+4*y(N*j-N+2)))/(2*data.h*data.ui-3*data.Dax(j)); % <<< From the boundary condition: z = 0 , C = Cinj + Dax/u * dC/dz
    y(N*j) = 4/3*y(N*j-1)-1/3*y(N*j-2); % <<< From the boundary condition: z = L , dC/dz=0
    
%     qast = data.isoT(y( (N*j-N+1):(N*j) ),data.parameter);
    
    % Forward finite differences
    DyDt(N*j-N+1) = data.Dax(j) * 1/(data.h^2)*(-y(N*j-N+4)+4*y(N*j-N+3)-5*y(N*j-N+2)+2*y(N*j-N+1)) - data.ui * 1/(2*data.h)*(-y(N*j-N+3)+4*y(N*j-N+2)-3*y(N*j-N+1)) - data.F * data.KLDF(j) * (qast(N*j-N+1)-y(N*nc+N*j-N+1)); %qast(N*j-N+1)
    DyDt(N*nc+N*j-N+1) = data.KLDF(j)*(qast(N*j-N+1)-y(N*nc+N*j-N+1));%qast(N*j-N+1)
    
    % Central finite differences
    for i=2:N-1
        DyDt(N*j-N+i) = data.Dax(j) * 1/(data.h^2)*(y(N*j-N+i+1)-2*y(N*j-N+i)+y(N*j-N+i-1)) - data.ui * 1/(2*data.h)*(y(N*j-N+i+1)-y(N*j-N+i-1)) - data.F * data.KLDF(j) * (qast(N*j-N+i)-y(N*nc+N*j-N+i)); %qast(N*j-N+i)
        DyDt(N*nc+N*j-N+i) = data.KLDF(j)*(qast(N*j-N+i)-y(N*nc+N*j-N+i));%qast(N*j-N+i)
    end
    
    % Backward finite differences
    DyDt(N*j) = 4/3*DyDt(N*j-1) - 1/3*DyDt(N*j-2); % <<< From the boundary condition: z = L , dC/dz=0;
    DyDt(N*nc+N*j) = data.KLDF(j)*(qast(N*j)-y(N*nc+N*j));%qast(N*j)

end



%% Isotherm: Langmuir
function q=langmuir(c,parameter)
% Multicomponent competitive Langmuir isotherm
%    Isotherm:   q1(C1,C2) = a1*C1/(1+b1*C1+b2*C2)
%    C is the concentration of the components [C1 C2 ... Cn]
%    parameter is an array containing the Langmuir parameters, like so: [a1 b1 ; a2 b2]

nc = length(c(1,:));
ndata = length(c(:,1));
q=zeros(ndata,nc);

langmuirDenominator = 1; % = 1 + b1*C1 + b2*C2 + ... + bn*Cn
for i = 1:nc
    langmuirDenominator = langmuirDenominator +  parameter(i,2)*c(:,i);
end

for i = 1:nc
    q(:,i) = parameter(i,1)*c(:,i) ./ langmuirDenominator;
end



%% Isotherm: Linear-Langmuir
function q=linearlangmuir(c,parameter)
% Multicomponent competitive Linear-Langmuir isotherm
%    Isotherm:   q1(C1,C2) = m1*C1 + a1*C1/(1+b1*C1+b2*C2)
%    C is the concentration of the components [C1 C2 ... Cn]
%    parameter is an array containing the Langmuir parameters, like so: [a1 b1 m1; a2 b2 m2]

nc = length(c(1,:));
ndata = length(c(:,1));
q=zeros(ndata,nc);

langmuirDenominator = 1; % = 1 + b1*C1 + b2*C2 + ... + bn*Cn
for i = 1:nc
    langmuirDenominator = langmuirDenominator +  parameter(i,2)*c(:,i);
end

for i = 1:nc
    q(:,i) = parameter(i,3)*c(:,i) + parameter(i,1)*c(:,i) ./ langmuirDenominator;
end


function q=linear(c,parameter)
% Multicomponent competitive Linear isotherm
%    Isotherm:   q1(C1) = H1*C1
%    C is the concentration of the components [C1 C2 ... Cn]
%    parameter is an array containing the Langmuir parameters, like so: [H1; H2]

nc = length(c(1,:));
ndata = length(c(:,1));
q=zeros(ndata,nc);



for i = 1:nc
    q(:,i) = parameter(i)*c(:,i);
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