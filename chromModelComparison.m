function chromModelComparison()
% This function compares the results obtained from diferent chromatographic
% models. Comment in/out different sections of this programm to test
% various cases.
clc
addpath('Models/');  % add models subfolder to path


%% INPUT DATA
isoType =   'linear';               % isotherm type, can be 'linear' or 'linear-langmuir'
feedProf =  'pulse';                % feed profile, can be 'pulse' or 'step'
parameter = 5.5;                    % isotherm parameters (depends on the isotherm model chosen)
L =         10;                     % cm, column length
Di =        1;                      % cm, column internal diameter
epsb =      0.4;                    % column bulk porosity
epsp =      0.5;                    % particle porosity
Q =         4;                      % mL/min, flow rate
Cfeed =     0.6*300/50;             % g/L, feed concentration
KLDF =      13.3*60/10;             % min-1, linear driving force (LDF) mass transfer coefficient
Dax =       5.57e-3;                % cm2/min, axial dispersion coefficient
tpulse =    50*0.001/4;             % min, feed pulse duration. For a step injection set tpulse = tfinal
tfinal =    10;                     % min, final time for calculation
opt.npz =   150;                    % number of discretization points in z
opt.npt =   100;                    % number of discretization points in t
opt.fig =   true;                   % true - show figures; false - do not show figures


%% Influence of the placeholder values used in the pdepe function for the boundary conditions
% Change the placeholder values in LDF_pdepe_1c function
% LDF_pdepe_1c(isoType,feedProf,parameter,L,Di,epsb,Q,Cfeed,KLDF,Dax,tpulse,tfinal,opt)


%% Comparing pdepe vs finite diferences for the LDF model
% tic
% LDF_df_1c(isoType,feedProf,parameter,L,Di,epsb,Q,Cfeed,KLDF,Dax,tpulse,tfinal,opt)
% fprintf('Calculation time: %.1f s \n', toc)
% 
% tic
% LDF_pdepe_1c(isoType,feedProf,parameter,L,Di,epsb,Q,Cfeed,KLDF,Dax,tpulse,tfinal,opt)
% fprintf('Calculation time: %.1f s \n', toc)


%% Comparing LDF model with TDM model
LDF_pdepe_1c(isoType,feedProf,parameter,L,Di,epsb,Q,Cfeed,KLDF,Dax,tpulse,tfinal,opt)

% LDF and TDM models are equivelent if:
H_TDM = (parameter-epsp)/(1-epsp);
KLDF_TDM = KLDF*(epsp+(1-epsp)*H_TDM);
TDMlinear_pdepe_1c(feedProf,H_TDM,L,Di,epsb,epsp,Q,Cfeed,KLDF_TDM,Dax,tpulse,tfinal,opt)


%% Comparing TDM model with Dax -> 0 with TM model (Dax neglected)
% Dax = 5.57e-1;
% TDMlinear_pdepe_1c(feedProf,parameter,L,Di,epsb,epsp,Q,Cfeed,KLDF,Dax,tpulse,tfinal,opt)
% 
% TMlinear_pdepe_1c(feedProf,parameter,L,Di,epsb,epsp,Q,Cfeed,KLDF,tpulse,tfinal,opt)


%% EDM model - DOES NOT WORK
% Dax =       5.57e-3;
% opt.npz =   350;
% opt.npt =   200; 
% ui = Q/(pi*Di^2/4);
% 
% % LDF and EDM models are equivelent if:
% H_EDM = (parameter-epsp)/(1-epsp);
% kmod = (1-epsb)/epsb*(epsp+(1-epsp)*H_EDM);
% Dapp = Dax + ((kmod/(1+kmod))^2) * epsb/(1-epsb)/KLDF*ui^2;
% EDMlinear_pdepe_1c(feedProf,H_EDM,L,Di,epsb,epsp,Q,Cfeed,Dapp,tpulse,tfinal,opt)


%% Comparing LDF_pdepe with LDF_pdepe_1c
% LDF_pdepe_1c(isoType,feedProf,parameter,L,Di,epsb,Q,Cfeed,KLDF,Dax,tpulse,tfinal,opt)

% LDF_pdepe(isoType,feedProf,parameter,L,Di,epsb,Q,Cfeed,KLDF,Dax,tpulse,tfinal,opt)

