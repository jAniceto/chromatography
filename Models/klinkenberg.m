function [t_, C, output] = klinkenberg(H, K, L, Di, epsb, Q, Cfeed, t, opt)
% Klinkenberg model
%
% Klinkenberg provides an useful approximation to the analytical solution 
% of the Convection-Dispersion model proposed by Anzelius for the case of 
% a single solute, an initially clean bed, frontal loading and negligible 
% axial dispersion. According to the Klinkenberg approximation the solute 
% concentration respect to axial distance and time is given by:
%
% $$\frac{C}{C_\mathrm{F}} \approx \frac{1}{2} [ 1 + \text{erf}(\sqrt{\tau} - \sqrt{\xi} + \frac{1}{8 \sqrt{\tau}} + \frac{1}{8 \sqrt{\xi}} ) ]$$
%
% $$\tau = K (t - \frac{z}{u_i})$$
% 
% $$\xi = \frac{K H z}{u_i} (\frac{1 - \varepsilon_b}{\varepsilon_b})$$
%
% ARGUMENTS
% ---------
% H : float
%   Linear isotherm parameter.
% K : float
%   Mass transfer coefficient.
% L : float
%   Column length.
% Di : float
%   Column diameter.
% epsb : float
%   Bed porosity.
% Q : float
%   Flow-rate.
% Cfeed : float
%   Feed concentration.
% t : float
%   Values of the time coordinate for which the model should be calculated.
% opt : struct
%   Options:
%       opt.fig (true | false) Whether to show plots.
%       plot.
% 
% REFERENCES
% ----------
% A. Klinkenberg, Ind. Eng. Chem., 46, 2285–2289 (1954).
%
% J.D. Seader, E.J. Henley, D.K. Roper. Separation Process Principles. 
% John Wiley & Sons, Inc., 3rd edition, 2011.


%% Default arguments (Example)
% If you run this function without arguments, by pressing ´Run´ or running
% ´klinkenberg()´ the following arguments will apply.
if nargin == 0
    H = 2.028;  % linear isotherm parameter
    K = 28.84;  % min-1, mass transfer coefficient
    L = 25;  % cm, column length
    Di = 0.46;  % cm, column internal diameter
    epsb = 0.335;  % column bulk porosity
    Q = 1;  % mL/min, flow rate
    Cfeed = 0.0569;  % g/L, feed concentration
    t = linspace(0, 25, 100)';  % min, time coordinate
    opt.fig = true;  % true - show figures; false - do not show figures
end


%% Data preparation
global data res

data = struct('H',H, 'K',K, 'L',L, 'Di',Di, 'epsb',epsb, 'Q',Q, 'Cfeed',Cfeed);
data.ndp = length(t); 
data.A = pi()*Di^2/4;  % cm2
data.F = (1-epsb) / epsb;
data.ui = Q / epsb / data.A;  % cm/min


%% Klinkenberg model
% Klinkenberg can not calculate value for t < L/ui. Remove those cases.
t_ = [];
for i = 1:data.ndp
    if t(i) > L / data.ui
        t_ = [t_ ; t(i)];
    end
end
output.n_removed_points = length(t) - length(t_);
output.n_final_points = length(t_);

% Klinkenberg
xi = H * K * L/data.ui * ((1-epsb)/epsb);
tau = K * (t_ - L/data.ui);
C = 0.5 * (1 + erf(sqrt(tau) - sqrt(xi) + 1/8./sqrt(tau) + 1/8./sqrt(xi))) * Cfeed;


%% Plot figures
if isfield(opt,'fig') && (opt.fig == 1)
    % Concentration history at column exit (Chromatogram)
    plot(t_, C, '-', 'LineWidth',1.5);
    xlabel('{\itt} (min)')
    ylabel('{\itC} (mg/mL)')
    box on;
end
