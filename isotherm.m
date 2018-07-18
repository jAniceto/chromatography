function q=isotherm(isoType, c, parameters)
% ISOTHERM  Defines the isotherm function
% Inputs:
% - isoType: String, selecting the isotherms type ('linear', 'langmuir', or 'linear-langmuir')
% - c: Matrix (or column vector), containing the concentrations in the liquid phase
% - parameters: Matrix (or line vector) containing the isotherm parameters. Add
% lines if more than one component.


%% VALIDATE INPUTS
% Check for valid isotherm type
availableIsotherms = {'linear','langmuir','linear-langmuir'};
if ~any(strcmp(availableIsotherms, isoType))
    error_str = strcat(sprintf('"%s" is not a valid isotherm type. Please select one from:', isoType), sprintf(' "%s"', availableIsotherms{:}));
    error(error_str);
end

% Check for valid number of isotherm parameters for the selected isotherm type
isothermParameters = containers.Map(availableIsotherms, [1 2 3]);  % specifies the number of parameters per isotherm type
if length(parameters(1,:)) ~= isothermParameters(isoType)
    error('Incorrect number of parameters provided in "parameters" for the isotherm type selected in "isoType".');
end


nc = length(parameters(:,1));

%% ISOTHERM LIBRARY
if strcmp(isoType,'linear')
% Linear isotherm 
% Single component: q = H * C
    parameters = parameters';
    q = parameters.*c;        

    
    
elseif strcmp(isoType,'langmuir')
% Single component: q = a * C / (1 + b * C)
% Multicompoent (competitive): q1 = a1 * C1 / (1 + b1 * C1 + b2 * C2)
    langmuirDenominator = ones(size(c));
    q = zeros(size(c));
    for i = 1:nc
        langmuirDenominator(:,i) = langmuirDenominator(:,i) +  parameters(i,2).*c(:,i);
    end
    for i = 1:nc
        q(:,i) = parameters(i,1).*c(:,i)./langmuirDenominator(:,i);
    end
    
    
elseif strcmp(isoType,'linear-langmuir') 
% Single component: q = H * C + a * C / (1 + b * C)
% Multicompoent (competitive): q1 = H1 * C1 + a1 * C1 / (1 + b1 * C1 + b2 * C2)
    langmuirDenominator = ones(size(c));
    q = zeros(size(c));
    for i = 1:nc
        langmuirDenominator(:,i) = langmuirDenominator(:,i) +  parameters(i,3).*c(:,i);
    end
    for i = 1:nc
        q(:,i) = parameters(i,1).*c(:,i) + parameters(i,2).*c(:,i)./langmuirDenominator(:,i);
    end
end
