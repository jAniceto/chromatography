function [t, C] = klinkenberg(H, L, Di, epsb, Q, Cfeed, KLDF, t, opt)
% Klinkenberg model
% Linear isotherm
% Change feed profile between pulse (e.g.: chromatografic peak) and step (e.g.: breakthrough experiment)
% Multicomponent

global data res

%% Default arguments (Example)
if nargin == 0
    H =         [2.028 2.297];    % isotherm parameters (depends on the isotherm model chosen)
    L =         25;                                                                 % cm, column length
    Di =        0.46;                                                                  % cm, column internal diameter
    epsb =      0.335;                                                              % column bulk porosity
    Q =         1;                                                                  % mL/min, flow rate
    Cfeed = 	[0.0569 0.0679];                             % g/L, feed concentration
    KLDF =      [28.84 11.78];                                  % min-1, linear driving force (LDF) mass transfer coefficient
    t =         {linspace(0,25,100)' linspace(0,25,100)'};
    opt.fig =   true;                                                               % true - show figures; false - do not show figures
end


%% Calculations
addpath('../')
data = struct('H',H','Cfeed',Cfeed,'Q',Q,'epsb',epsb,'KLDF',KLDF);
nc = length(Cfeed);
data.nc = nc;
A = pi()*Di^2/4; % cm2
data.F = (1-epsb)/epsb;
data.ui = Q/epsb/A; % cm/min

% Klinkenberg
for i = 1:data.nc
    xi = H(i)*KLDF(i)*L/data.ui*((1-epsb)/epsb);
    tau = KLDF(i)*(t{i}-L/data.ui);
    
    time_points_to_cut = 1;
    for j = 1:length(tau)
        if tau(j) < 0
            time_points_to_cut = time_points_to_cut + 1;
        end
    end
    t{i} = t{i}(time_points_to_cut:end);
    tau = tau(time_points_to_cut:end);
    
    C{i} = 0.5*(1+erf(sqrt(tau)-sqrt(xi)+1/8./sqrt(tau)+1/8./sqrt(xi)))*Cfeed(i);
end
res.t = t'; res.C = C;


% Plot figures
if isfield(opt,'fig') && (opt.fig == 1)
    % Concentration history at column exit (Chromatogram)
    h1 = figure;
    hold all
    for i = 1:data.nc
        plot(t{i},C{i},'-','LineWidth',1.5);
    end
    hold off;
    xlabel('{\itt} (min)')
    ylabel('{\itC} (mg/mL)')
    box on;
end
