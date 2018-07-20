function [parameter, C, q, AARDev] = fitIsotherm(C, q, parameter, isotherm, opt)
% ISOTHERMFIT  Fits an adsorption isotherm to experimental data. 
%
%   ISOTHERMFIT(C, q, parameter, isotherm, opt) returns the fitted parameters, the fitted isotherm plot, and the AARD(%) of the fitting.
%   
%   C          Concentration in the fluid phase.
%
%   q          Concentration in the solid phase (adsorbent).
%
%   parameter  Matrix containing initial estimations for the isotherm parameters.
%
%   isotherm   string containing the isotherm to use. Default is Langmuir isotherm. Other options are:
%                 'langmuir' for Langmuir isotherm: q = parameter(1)*C/(1+parameter(2)*C)
%                 'freundlich' for Freundlich isotherm: q = parameter(1)*C^(1/parameter(2))
%                 'linearlangmuir' for Linear-Langmuir isotherm: q = parameter(1)*C/(1+parameter(2)*C) + parameter(3)*C;
%                 'bilangmuir' for Bi-Langmuir isotherm: q = parameter(1)*C/(1+parameter(2)*C) + parameter(3)*C/(1+parameter(4)*C);
%                 'langmuirfreundlich' for Langmuir-Freundlich isotherm: q = parameter(1)*C/(1+parameter(2)*C^(1/parameter(3)))
%                 'toth' for Toth isotherm: q = parameter(1)*C/((1+(parameter(2)*C)^parameter(3))^(1/parameter(3)))
%                 'linear' for linear isotherm: q = parameter(1)*C
%
%  opt         options structure
%              opt.objectFun  Select the objective function to use in optimization. Default is an Least-Squares function. Option are:
%                 opt.objectFun = 'aard' to use an Absolute Average Relative Deviation objective function: f = sum(abs(qexp-qcalc)/qexp)/length(qexp)
%                 opt.objectFun = 'ls' to use a Least-Squares objective function: fval = sum((qexp-qcalc)^2)
%              opt.fig  Show figures if true (default) or don't show figures if false
%
% See also ...

%% Validating inputs and setting default options
if nargin < 3
    fprintf(2,'Missing input arguments. Type "help fitIsotherm" for more information.\n');
    return;
elseif nargin == 3
    isotherm = 'langmuir';
    opt.objectFun = 'ls';
    opt.figs = true;
elseif nargin == 4
    opt.objectFun = 'ls';
    opt.fig = true;
end
if ~isfield(opt, 'fig')
    opt.fig = true;
end
if ~isfield(opt, 'objectFun')
    opt.objectFun = 'ls';
end

fprintf('\nIsotherm: %s', isotherm);

%% Optimization
data.C = C;
data.q = q;

% If using the aard objective function, find zeros in q array and substitute by 1e-8
if strcmp(opt.objectFun,'aard')
    data.q(~data.q) = 1e-8;
end

% Creating a function handle to work with any isotherm
isoT = str2func(isotherm);
objF = str2func(opt.objectFun);

% Parameter optimization using fminsearch
[parameter, fval, exitflag] = fminsearch(objF, parameter, [], data, isoT);

C = 0:(max(data.C)/(250-1)):max(data.C);
q = isoT(C, parameter);

% Ploting the isotherm with the fitted parameter
if opt.fig
    figure
    plot(data.C, data.q, 'o');
    hold on;
    plot(C, q, '-', 'LineWidth', 2);
    title(isotherm);
    xlabel('\itC')
    ylabel('\itq')
    hold off;
end


% Calculation of AARD
qcalc = isoT(data.C,parameter);
data.q(~data.q) = 1e-8;
qcalc(~qcalc) = 1e-8;
AARDev = sum(abs(data.q-qcalc)./data.q)/length(data.C)*100;

% Print results
if exitflag==1
    fprintf('\nOptimization successful.\n');
else
    fprintf('\nOptimization failed (exitflag = %i)\n', exitflag);
end

if strcmp(isotherm,'langmuir')
    fprintf('a = %f\nb = %f\n',parameter(1),parameter(2));
elseif strcmp(isotherm,'freundlich')
    fprintf('K = %f\nn = %f\n',parameter(1),parameter(2));
elseif strcmp(isotherm,'linearlangmuir')
    fprintf('m = %f\na = %f\nb = %f\n',parameter(3),parameter(1),parameter(2));
elseif strcmp(isotherm,'bilangmuir')
    fprintf('a1 = %f\nb1 = %f\na2 = %f\nb2 = %f\n',parameter(1),parameter(2),parameter(3),parameter(4));
elseif strcmp(isotherm,'langmuirfreundlich')
    fprintf('a = %f\nb = %f\nn = %f\n',parameter(1),parameter(2),parameter(3));
elseif strcmp(isotherm,'toth')
    fprintf('a = %f\nb = %f\nv = %f\n',parameter(1),parameter(2),parameter(3));
elseif strcmp(isotherm,'linear')
    fprintf('H = %f\n',parameter(1));
end

fprintf('AARD = %.2f %%\n\n',AARDev)



%% Possible objective funtions to use: AARD; Least squares
function f=aard(parameter,data,isoT)
% Calculated q using the desired isotherm equation
qcalc = isoT(data.C,parameter);

% Objective function for parameter optimization (AARD)
f = sum(abs(data.q-qcalc)./data.q)/length(data.C);


function f=ls(parameter,data,isoT)
% Calculated q using the desired isotherm equation
qcalc = isoT(data.C,parameter);

% Objective function for parameter optimization (Least-Squares)
f = sum((data.q-qcalc).^2);



%% Possible isotherms to adjust
function q=langmuir(C,parameter)
q = parameter(1)*C./(1+parameter(2)*C);


function q=freundlich(C,parameter)
q = parameter(1)*C.^(1/parameter(2));


function q=linearlangmuir(C,parameter)
q = parameter(1)*C./(1+parameter(2)*C) + parameter(3)*C;


function q=bilangmuir(C,parameter)
q = parameter(1)*C./(1+parameter(2)*C) + parameter(3)*C./(1+parameter(4)*C);


function q=langmuirfreundlich(C,parameter)
q = parameter(1)*C./(1+parameter(2)*C.^(1/parameter(3)));


function q=toth(C,parameter)
q = parameter(1)*C./((1+(parameter(2)*C).^parameter(3)).^(1/parameter(3)));


function q=linear(C,parameter)
q = parameter(1)*C;