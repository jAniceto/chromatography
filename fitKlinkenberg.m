function fitKlinkenberg(exp_tc, exp_Cfeed, parameter, L, Di, epsb, Q, tfinal, opt)
% Fit the Klinkenberg model to the data provided
% Two parameters are fitted: H and KLDF
% Any number of diferent experiments can be used by expanding the exp_Cfeed
% and exp_tc to include additional data. Parameters are fitted to all data
% simultaneously.

global spde_count


%% Default arguments (Example)
if nargin == 0
    parameter = [2 29];         % [ H  KLDF(min-1) ], initial estimate for optimization
    L =         25;             % cm, column length
    Di =        0.46;           % cm, column internal diameter
    epsb =      0.335;          % column bulk porosity
    Q =         1;              % mL/min, flow rate
    tfinal =    25;
    opt.fig =   false;
    
    %% Experimental data:
    exp_Cfeed(1) = 0.100;                % g/L, feed concentration
    exp_tc{1} =    [2.36	0
                    4.26	0
                    5.06	0
                    6.16	0
                    6.56	0.002426767
                    6.86	0.019353929
                    7.16	0.045096765
                    7.42	0.059907137
                    8.06	0.079553852
                    8.66	0.091062002
                    9.06	0.095583065
                    9.66	0.096728626
                    10.67	0.102127158
                    12.16	0.101030992
                    16.16	0.098891963
                    18.66	0.099783668
                    21.66	0.101550143
                    23.66	0.097753721];

    exp_Cfeed(2) = 0.502;                % g/L, feed concentration
    exp_tc{2} =    [2.16	0.003
                    3.66	0.002
                    5.86	0.002
                    6.16	0.021
                    6.46	0.127
                    7.36	0.381
                    7.91	0.451
                    9.06	0.459
                    10.36	0.500
                    11.46	0.519
                    12.46	0.505
                    16.03	0.454
                    18.96	0.469
                    23.18	0.484];
end


%% Calculations
addpath('Models/');  % add models subfolder to path
tic
data = struct('L',L, 'Di',Di, 'epsb',epsb, 'Q',Q);
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
options = optimset('PlotFcns',@optimplotfval,'TolFun',1e-8,'TolX',1e-8);
[parameter, fval, exitflag] = fminsearch(@fobj, parameter, options, exp, data, opt)

% Calculate concentration history with optimized parameters
time = {linspace(0,tfinal,500)' linspace(0,tfinal,500)'};
[t, C] = klinkenberg(ones(1,data.ndata)*parameter(1), data.L, data.Di, data.epsb, data.Q, exp.Cfeed, ones(1,data.ndata)*parameter(2), time, opt);

% Plot concentration history at column exit
figure(10)
hold on;
for j = 1:data.ndata
    plot(exp_tc{j}(:,1),exp_tc{j}(:,2),'o','MarkerSize',10);
end
ax = gca; ax.ColorOrderIndex = 1;
for j = 1:data.ndata
    plot(t{j}, C{j}, '-','LineWidth',1.5);
end
hold off;
set(gca,'fontsize',12)
xtickformat('%.0f'); ytickformat('%.2f');
xlabel('{\itt} (min)')
ylabel('{\itC} (mg/mL)')
    
% Remove C = 0 values from experimental points in order to calculate AARD
for j = 1:data.ndata
    exp_tc_aard{j} = exp_tc{j}(all(exp_tc{j},2),:);
end

for j = 1:data.ndata
    time_aard{j} = exp_tc_aard{j}(:,1);
end
[t, Conc_aard] = klinkenberg(ones(1,data.ndata)*parameter(1), data.L, data.Di, data.epsb, data.Q, exp.Cfeed, ones(1,data.ndata)*parameter(2), time_aard, opt);

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

if (parameter(1) <= 0) || (parameter(2) <= 0)
    f = 1000000;

else
    [t, C] = klinkenberg(ones(1,data.ndata)*parameter(1), data.L, data.Di, data.epsb, data.Q, exp.Cfeed, ones(1,data.ndata)*parameter(2), exp.t, opt);

    % Objective function for parameter optimization (Least-Squares)
    f = 0;
    for j = 1:data.ndata
        f = f + sum( (exp.c{j}-C{j}).^2 );
    end
end

% Info to show at each iteration
spde_count = spde_count + 1;
fprintf('spde count = %i  |  fobj = %f  |  [H KLDF] = [%f %f]\n',spde_count, f, parameter(1), parameter(2));
