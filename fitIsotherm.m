function fitIsotherm(C, q, parameter, isotherm, objectFun)
% ISOTHERMFIT  Fits an adsorption isotherm to experimental data. 
%
%   ISOTHERMFIT(C,q,parameter,isotherm,objectFun) returns the fitted parameters, the fitted isotherm plot, and the AARD(%) of the fitting.
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
%  objectFun  string containing the objective function to use in optimization. Default is an Absolute Average Relative Deviation objective function. Other option are:
%                 'aard' to use an Absolute Average Relative Deviation objective function: f = sum(abs(qexp-qcalc)/qexp)/length(qexp)
%                 'ls' to use a Least-Squares objective function: fval = sum((qexp-qcalc)^2)
%
% See also ...

%% Validating inputs and setting default options
if nargin < 3
    fprintf(2,'\nMissing input arguments. Type "help fitIsotherm" for more information.\n');
    return;
elseif nargin == 3
    isotherm = 'langmuir';
    objectFun = 'ls';
elseif nargin == 4
    objectFun = 'ls';
end
fprintf('\n\n\nIsotherm: %s\n',isotherm);

%% Optimization
data.C = C;
data.q = q;

% If using the aard objective function, find zeros in q array and substitute by 1e-8
if strcmp(objectFun,'aard')
    data.q(~data.q) = 1e-8;
end

% Creating a function handle to work with any isotherm
isoT = str2func(isotherm);
objF = str2func(objectFun);

% Parameter optimization using fminsearch
[parameter, fval, exitflag] = fminsearch(objF, parameter, [], data, isoT);

% Ploting the isotherm with the fitted parameter
C = 0:(max(data.C)/(100-1)):max(data.C);
q = isoT(C,parameter);
if ishandle(51)==0
    figure(51)
    plot(data.C,data.q,'o');
    hold all
    plot(C,q,'-');

elseif ishandle(51)==1
    hold all
    plot(C,q,'-')
end
xlabel('C')
ylabel('q')

% Calculation of AARD
qcalc = isoT(data.C,parameter);
data.q(~data.q)=1e-8;
qcalc(~qcalc)=1e-8;
AARDev = sum(abs(data.q-qcalc)./data.q)/length(data.C)*100;

% Print results
if exitflag==1
    fprintf('\nOptimization successful.\n');
else
    fprintf('\nOptimization failed (exitflag = %i)\n',exitflag);
end

% fprintf(strcat('\nparameter = [ ',num2str(parameter),' ]'));

if strcmp(isotherm,'langmuir')
    fprintf('\na = %f\nb = %f\n',parameter(1),parameter(2));
elseif strcmp(isotherm,'freundlich')
    fprintf('\nK = %f\nn = %f\n',parameter(1),parameter(2));
elseif strcmp(isotherm,'linearlangmuir')
    fprintf('\nm = %f\na = %f\nb = %f\n',parameter(3),parameter(1),parameter(2));
elseif strcmp(isotherm,'bilangmuir')
    fprintf('\na1 = %f\nb1 = %f\na2 = %f\nb2 = %f\n',parameter(1),parameter(2),parameter(3),parameter(4));
elseif strcmp(isotherm,'langmuirfreundlich')
    fprintf('\na = %f\nb = %f\nn = %f\n',parameter(1),parameter(2),parameter(3));
elseif strcmp(isotherm,'toth')
    fprintf('\na = %f\nb = %f\nv = %f\n',parameter(1),parameter(2),parameter(3));
elseif strcmp(isotherm,'linear')
    fprintf('\nH = %f\n',parameter(1));
end

fprintf('\nAARD = %.2f\n',AARDev)



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