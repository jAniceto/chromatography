function fit_klinkenberg(exp_tc, exp_Cfeed, H_est, K_est, L, Di, epsb, Q, tfinal, opt)
% Fit the Klinkenberg model (H and K parameters) to the data provided. 
%
% ARGUMENTS
% ---------
% exp_tc : array of n-by-2 matrices
%   Data to fit the model. The array has m elements, each corresponding to
%   a dataset. Each dataset has a n-by-2 matrix with t and c values
%   (columns) and n datapoints (lines).
% exp_Cfeed : vector with m elements (one per dataset)
%   Feed concentration.
% H_est : float
%   Initial estimate for the Klinkenberg linear equilibrium constant.
% K_est : float
%   Initial estimate for the Klinkenberg mass transfer coefficient.
% L : float
%   Column length.
% Di : float
%   Column diameter.
% epsb : float
%   Bed porosity.
% Q : float
%   Flow-rate.
% tfinal : float
%   Final time for calculation. Should be higher than the highest time in
%   the experimental datasets.
% opt : struct
%   Options:
%       opt.fig (true | false) specify if the model function should make a
%       plot.
% 
% NOTES
% -----
% Any number of diferent experiments ([t c] datasets) can be used by 
% expanding the ´exp_Cfeed´ and ´exp_tc´ arrays to include additional data. 
% Parameters are fitted to all data simultaneously.


%% Default arguments (Example)
if nargin == 0
    H_est = 2;  % initial estimate for the linear isotherm parameter
    K_est = 29;  % initial estimate for the mass transfer coefficient
    L = 25;  % cm, column length
    Di = 0.46;  % cm, column internal diameter
    epsb = 0.335;  % column bulk porosity
    Q = 1;  % mL/min, flow rate
    tfinal = 25; 
    opt.fig = false;
    
    % Experimental data
    exp_Cfeed(1) = 0.100; % g/L, feed concentration
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

    exp_Cfeed(2) = 0.502; % g/L, feed concentration
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
global fobj_evals

addpath('models/');  % add models subfolder to path
tic
data = struct('L',L, 'Di',Di, 'epsb',epsb, 'Q',Q);
exp = struct('Cfeed', exp_Cfeed);
data.n = length(exp_tc); % Number of data sets ([t C] pairs)   

% Check for invalid data (NaN)
for j = 1:data.n
    exp_tc{j} = exp_tc{j}(~any(isnan(exp_tc{j}),2),:);
end

% Separate t and C data
for j = 1:data.n
    exp.t{j} = exp_tc{j}(:,1);
    exp.c{j} = exp_tc{j}(:,2);
end

% Parameter optimization using fminsearch
fobj_evals = 0;
params = [H_est, K_est];
options = optimset('PlotFcns',@optimplotfval, 'TolFun',1e-8, 'TolX',1e-8);
[parameter, fval, exitflag, outputs] = fminsearch(@fobj, params, options, exp, data, opt)

H = parameter(1);
K = parameter(2);


%% Plot concentration history at column exit
time = linspace(0,tfinal,500)';

% Run the model for the optimized parameters for a fine t vector
for j = 1:data.n
    [t, C, out] = klinkenberg(H, K, data.L, data.Di, data.epsb, data.Q, exp.Cfeed(j), time, opt);
    tmodel{j} = t;
    Cmodel{j} = C;
end

% Plot experimental data
figure()
hold on;
for j = 1:data.n
    plot(exp.t{j}, exp.c{j}, 'o', 'MarkerSize',10);
end

% Plot model
ax = gca; ax.ColorOrderIndex = 1;  % reset colors
for j = 1:data.n
    plot(tmodel{j}, Cmodel{j}, '-', 'LineWidth',1.5);
end
hold off;
set(gca, 'fontsize',12)
xtickformat('%.0f'); 
ytickformat('%.2f');
xlabel('{\itt} (min)')
ylabel('{\itC} (mg/mL)')


%% Calculate model metrics
% Run the model for the optimized parameters for the experimental values of t vector
for j = 1:data.n
    [t, C, out] = klinkenberg(H, K, data.L, data.Di, data.epsb, data.Q, exp.Cfeed(j), exp.t{j}, opt);
    Cmodel{j} = C;
    Cexp{j} = exp.c{j}(out.n_removed_points+1:end);
end

% AARD (%)
for j = 1:data.n
    aard{j} = abs(Cexp{j} - Cmodel{j}) ./ Cexp{j} * 100;

    fprintf('AARD%i =\n', j);
    disp(aard{j});
end

% AARD global (%)
aard_all = [];
for j = 1:data.n
    aard_tmp = aard{j}(~any( isnan( aard{j} ) | isinf( aard{j} ), 2 ), :);  % remove Inf and NaN values before suming the errors
    aard_all = [aard_all ; aard_tmp];
end
aard_total = mean(aard_all);
fprintf('AARDtotal = %.2f %% \n\n', aard_total);

toc


%% Objective funtion
function f = fobj(parameter, exp, data, opt)
global fobj_evals

H = parameter(1);
K = parameter(2);

if (parameter(1) <= 0) || (parameter(2) <= 0)
    f = 1000000;

else
    
    for j = 1:data.n
        [t, C, out] = klinkenberg(H, K, data.L, data.Di, data.epsb, data.Q, exp.Cfeed(j), exp.t{j}, opt);
        Ccalc{j} = C;
        Cexp{j} = exp.c{j}(out.n_removed_points+1:end);
    end
    
    % Objective function for parameter optimization (least-squares)
    f = 0;
    for j = 1:data.n
        f = f + sum( (Cexp{j} - Ccalc{j}).^2 );
    end
end

% Info to show at each iteration
fobj_evals = fobj_evals + 1;
fprintf('%i  |  fobj = %f  |  [H K] = [%f %f]\n', fobj_evals, f, parameter(1), parameter(2));
