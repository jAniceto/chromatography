%% Isotherm Library

function q=isotherm(isoType, c, parameter)
% Defines the isotherm

nc = length(c);

if strcmp(isoType,'linear')
% Linear isotherm:   q1(C1) = H*C1, parameter = H
    for i = 1:nc
        q(i) = parameter(i)*c(i);        
    end 
    
    
elseif strcmp(isoType,'langmuir')
% Langmuir isotherm:   q1(C1,C2) = a1*C1/(1+b1*C1) , parameter = [a1 b1]
    langmuirDenominator = 1;
    for i = 1:nc
        langmuirDenominator = langmuirDenominator +  parameter(2,i)*c(i);
    end
    for i = 1:nc
        q(i) = parameter(1,i)*c(i)/langmuirDenominator;        
    end
    
    
elseif strcmp(isoType,'linear-langmuir') 
% Linear-Langmuir isotherm:   q1(C1,C2) = m1*C1 + a1*C1/(1+b1*C1) , parameter = [a1 ; b1 ; m1]
    langmuirDenominator = 1;
    for i = 1:nc
        langmuirDenominator = langmuirDenominator +  parameter(2,i)*c(i);
    end
    for i = 1:nc
        q(i) = parameter(3,i)*c(i) + parameter(1,i)*c(i)/langmuirDenominator;        
    end
   
else
    error('Invalid isotherm type. isoType must be "linear" or "langmuir" or "linear-langmuir"')
end