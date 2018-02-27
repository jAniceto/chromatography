function fitModel_isotherm_KLDF(exp_tc, exp_Cfeed, isoType, feedProf, parameter, L, Di, epsb, Q, Dax, tpulse, tfinal, opt)
% Fit the chromatographic model to the data provided
% Two parameters are fitted: H and KLDF
% Any number of diferent experiments can be used by expanding the exp_Cfeed
% and exp_tc to include additional data. Parameters are fitted to all data
% simultaneously.

global spde_count


%% Default arguments (Example)
if nargin == 0
    parameter = [2 10];         % [ H  KLDF(min-1) ], initial estimate for optimization
    isoType =   'linear';       % isotherm model
    feedProf =  'pulse';        % feed profile, can be 'pulse' (e.g.: chromatografic peak) or 'step' (e.g.: breakthrough experiment)
    L =         25;             % cm, column length
    Di =        0.46;           % cm, column internal diameter
    epsb =      0.335;          % column bulk porosity
    Q =         1;              % mL/min, flow rate
    Dax =       4.779E-3;       % cm2/min, axial dispersion coefficient
    tpulse =    28;             % min, feed pulse duration. For a step injection set tpulse = tfinal
    tfinal =    28;             % min, final time for calculation
    opt.npz =   50;             % number of discretization points in z
    opt.npt =   30;             % number of discretization points in t
    
    %% Experimental data:
    exp_Cfeed(1) = 0.100;                % g/L, feed concentration
    exp_tc{1} =    [0.66	0.0000
                    2.36	0.0000
                    4.26	0.0000
                    5.06	0.0000
                    6.16	0.0000
                    6.56	0.0024
                    6.86	0.0194
                    7.16	0.0451
                    7.42	0.0599
                    8.06	0.0796
                    8.66	0.0911
                    9.06	0.0956
                    9.66	0.0967
                    10.67	0.1021
                    12.16	0.1010
                    16.16	0.0989
                    18.66	0.0998
                    21.66	0.1016
                    23.66	0.0978];
    exp_Cfeed(2) = 0.505;                % g/L, feed concentration
    exp_tc{2} =    [0.66	0.0000
                    1.16	0.0000
                    2.16	0.0000
                    3.66	0.0000
                    5.86	0.0021
                    6.16	0.0209
                    6.46	0.1268
                    7.36	0.3810
                    7.91	0.4510
                    9.06	0.4589
                    10.36	0.5003
                    11.46	0.5190
                    12.46	0.5048
                    18.96	0.4691];
    exp_Cfeed(3) = 1.004;                % g/L, feed concentration
    exp_tc{3} =    [1.16	0.0000
                    2.16	0.0000
                    3.66	0.0000
                    5.66	0.0022
                    6.66	0.1834
                    6.96	0.3019
                    7.26	0.4066
                    7.56	0.5131
                    7.86	0.5765
                    8.86	0.7138
                    10.06	0.8132
                    11.06	0.9165
                    12.66	0.9747
                    14.06	0.9569
                    17.26	0.9243];
end


%% Calculations
addpath('Models/');  % add models subfolder to path
tic
data = struct('isoType',isoType, 'feedProf',feedProf, 'L',L, 'Di',Di, 'epsb',epsb, 'Q',Q, 'tpulse',tpulse, 'tfinal',tfinal, 'Dax',Dax);
exp = struct('Cfeed',exp_Cfeed);
data.ndata = length(exp_tc); % Number of data sets ([t C] pairs)   

% Check for invalid data (NaN)
for j = 1:data.ndata
    exp_tc{j} = exp_tc{j}(~any(isnan(exp_tc{j}),2),:);
end

% Separate t and C data
for j = 1:data.ndata
    exp.t{j} = exp_tc{j}(:,1);
    exp.c{j} = exp_tc{j}(:,2);
end

% Parameter optimization using fminsearch
spde_count = 0;
options = optimset('PlotFcns',@optimplotfval,'TolFun',1e-6,'TolX',1e-6);
[parameter, fval, exitflag] = fminsearch(@fobj, parameter, options, exp, data, opt)

% Calculate concentration history with optimized parameters
[sol, t, z, C] = LDF_pdepe(data.isoType, data.feedProf, ones(1,data.ndata)*parameter(1), data.L, data.Di, data.epsb, data.Q, exp.Cfeed, ones(1,data.ndata)*parameter(2), ones(1,data.ndata)*data.Dax, data.tpulse, data.tfinal, opt);

Conc = zeros(length(t),data.ndata);
for i=1:data.ndata
    Conc(:,i) = sol(:,end,i); 
end

% Plot concentration history at column exit
figure
hold on;
for j = 1:data.ndata
    plot(exp_tc{j}(:,1),exp_tc{j}(:,2),'o');
end
for j = 1:data.ndata
    plot(t,Conc(:,j),'-');
end
hold off;
xlabel('t (min)')
ylabel('C (mg/mL)')

% Remove C = 0 values from experimental points in order to calculate AARD
for j = 1:data.ndata
    exp_tc_aard{j} = exp_tc{j}(all(exp_tc{j},2),:);
end

% Evaluate final solution (with optimized parameters) at the experimental times
for j = 1:data.ndata
    [Cout, dCoutdt] = pdeval(0, t, sol(:,end,j), exp_tc_aard{j}(:,1));
    Conc_aard{j} = Cout;
end

for j = 1:data.ndata
    AARD{j} = abs(Conc_aard{j}-exp_tc_aard{j}(:,2))./exp_tc_aard{j}(:,2);
end
for j = 1:data.ndata
    fprintf('AARD%i =\n', j);
    disp(AARD{j});
end

AARDtotal = 0;
for j = 1:data.ndata
    AARDtotal = AARDtotal + sum(AARD{j})/length(AARD{j});
end
fprintf('AARDtotal = %.4f \n\n', AARDtotal);

toc


%% Objective funtion
function f = fobj(parameter, exp, data, opt)
global spde_count

[sol, t, z, C] = LDF_pdepe(data.isoType, data.feedProf, ones(1,data.ndata)*parameter(1), data.L, data.Di, data.epsb, data.Q, exp.Cfeed, ones(1,data.ndata)*parameter(2), ones(1,data.ndata)*data.Dax, data.tpulse, data.tfinal, opt);


% Evaluate ODE solution at the experimental t
for j = 1:data.ndata
    [Cout, dCoutdt] = pdeval(0, t, sol(:,end,j), exp.t{j});
    Conc{j} = Cout;
end

% Objective function for parameter optimization (Least-Squares)
f = 0;
for j = 1:data.ndata
    f = f + sum( (exp.c{j}-Conc{j}).^2 );
end

% Info to show at each iteration
spde_count = spde_count + 1;
fprintf('spde count = %i  |  fobj = %f  |  [H KLDF] = [%f %f]\n',spde_count, f, parameter(1), parameter(2));

% Show current figure
% if mod(spde_count,10)==0
%     t = sol.x';
%     Conc = zeros(length(t),data.ndata);
%     for i=1:data.ndata
%         Conc(:,i) = sol.y(i*opt.npz,:)';
%     end
    
%     figure(5)
%     clf(5)
%     plot(exp_tc1(:,1),exp_tc1(:,2),'o',exp_tc2(:,1),exp_tc2(:,2),'o',exp_tc3(:,1),exp_tc3(:,2),'o');
%     for i=1:data.ndata
%         hold on
%         plot(t,Conc(:,i),'-');
%     end
%     drawnow;
% end
