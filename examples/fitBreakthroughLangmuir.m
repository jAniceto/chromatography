function fitBreakthroughLangmuir()


parameters = [5.906	6.243 3.76];        % [ H  KLDF(min-1) ], initial estimate for optimization
isoType =   'langmuir';       % isotherm model
feedProf =  'pulse';        % feed profile, can be 'pulse' (e.g.: chromatografic peak) or 'step' (e.g.: breakthrough experiment)
L =         35.4;             % cm, column length
Di =        1.95;           % cm, column internal diameter
epsb =      0.34;          % column bulk porosity
Q =         21.8;              % mL/min, flow rate
Dax =       1.98;       % cm2/min, axial dispersion coefficient
tpulse =    15;             % min, feed pulse duration. For a step injection set tpulse = tfinal
tfinal =    15;             % min, final time for calculation
opt.npz =   80;             % number of discretization points in z
opt.npt =   50;             % number of discretization points in t
opt.optim_fig = true;

%% Experimental data:
exp_Cfeed(1) = 2.83;                % g/L, feed concentration
exp_tc{1} =    [0.907	0.011
                1.411	0.026
                1.915	0.026
                2.419	0.026
                2.894	0.027
                3.403	0.293
                3.908	2.647
                4.413	3.001
                4.919	2.870
                5.392	3.072
                5.898	3.076
                6.403	3.009
                6.908	2.899
                7.413	3.027
                7.887	3.058
                8.392	3.025
                8.897	3.011
                9.413	2.895
                9.918	3.019
                10.392	3.100
                10.897	2.899
                11.402	3.024
                11.907	2.869
                12.402	2.727
                12.897	2.811
                13.402	2.825
                13.907	2.864
                14.349	2.774];

            
%% RUN EXAMPLE ABOVE
addpath('../', '../Models/');  % add parent folder to path

fitModel_isotherm_KLDF(exp_tc, exp_Cfeed, isoType, feedProf, parameters, L, Di, epsb, Q, Dax, tpulse, tfinal, opt)