function plotIsotherms()

x = linspace(0, 10, 100)';


%% Linear isotherm 
% Single component: q = H * C

parameters = [2.5 ; 3.5 ; 4];

c = [x x x];
nc = length(parameters(:,1));
q = isotherm('linear', c, parameters);

figure
plot(c, q, '-');
xlabel('\it{C_{i}}')
ylabel('\it{q_{i}}')
title('Linear isotherms');
legend_str{1} = '';
for i = 1:nc
    legend_str{i} = sprintf('Component %i: H = %.1f', i, parameters(i,1));
end
legend(legend_str, 'location', 'northwest')



%% Langmuir isotherm
% Single component: q = a * C / (1 + b * C)
% Multicompoent (competitive): q1 = a1 * C1 / (1 + b1 * C1 + b2 * C2)

parameters = [2.5 0.5;  % isotherm constants for component 1: [a1 b1]
              3.5 0.4;  % isotherm constants for component 2: [a2 b2]
              4 0.2];   % isotherm constants for component 3: [a3 b3]

c = [x x x];
nc = length(parameters(:,1));
q = isotherm('langmuir', c, parameters);

figure
plot(c, q, '-');
xlabel('{\itC}_{i}')
ylabel('{\itq}_{i}')
title('Langmuir isotherms (competitive)');
legend_str{1} = '';
for i = 1:nc
    legend_str{i} = sprintf('Component %i: a = %.1f, b = %.1f', i, parameters(i,1), parameters(i,2));
end
legend(legend_str, 'location', 'northwest')


%% Linear-Langmuir isotherm
% Single component: q = H * C + a * C / (1 + b * C)
% Multicompoent (competitive): q1 = H1 * C1 + a1 * C1 / (1 + b1 * C1 + b2 * C2)

parameters = [1.5 2.5 0.5; % isotherm constants for component 1: [H1 a1 b1]
              2.5 3.5 0.3; % isotherm constants for component 2: [H2 a2 b2]
              3.0 4 0.2];  % isotherm constants for component 2: [H3 a3 b3]

c = [x x x];
nc = length(parameters(:,1));
q = isotherm('linear-langmuir', c, parameters);

figure
plot(c, q, '-');
xlabel('{\itC}_{i}')
ylabel('{\itq}_{i}')
title('Linear-Langmuir isotherms (competitive)');
legend_str{1} = '';
for i = 1:nc
    legend_str{i} = sprintf('Component %i: H = %.1f, a = %.1f, b = %.1f', i, parameters(i,1), parameters(i,2), parameters(i,3));
end
legend(legend_str, 'location', 'northwest')


