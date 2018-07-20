function fitIsotherms()

%% Experimental data
Cq{1} = [   0.111	0.325
            0.127	0.387
            0.532	0.727
            0.914	0.856
            1.423	0.830
            1.470	0.850
            2.893	0.882
            2.996	0.870];
        
Cq{2} = [   0	0
            0.142	0.459
            0.153	0.374
            0.649	0.599
            0.672	0.572
            1.400	0.737
            2.073	0.820
            2.162	0.796
            2.945	0.809
            3.530	0.775
            3.596	0.827
            3.696	0.779
            3.729	0.799 ];

Cq{3} = [   0	0
            0.290	0.392
            0.306	0.336
            0.567	0.550
            1.141	0.543
            1.274	0.590
            2.554	0.587
            2.593	0.591
            3.868	0.642
            5.624	0.683 ];

        
%% Run isotherm fitting
addpath('../');  % add parent folder to path

% Fit Langmuir isotherm to all data series
for i = 1:length(Cq)
    C = Cq{i}(:,1);
    q = Cq{i}(:,2);
    fitIsotherm(C, q, [1 2], 'langmuir');
end

% Fit all isotherms to Cq{3} data a plot results together
opt.fig = false;
C = Cq{3}(:,1);
q = Cq{3}(:,2);

figure
plot(C, q, 'o');
hold on;
[parameter, Ccalc, qcalc, AARDev] = fitIsotherm(C, q, 0.8, 'linear', opt);
plot(Ccalc, qcalc, '-', 'LineWidth', 2);
[parameter, Ccalc, qcalc, AARDev] = fitIsotherm(C, q, [1 2], 'langmuir', opt);
plot(Ccalc, qcalc, '-', 'LineWidth', 2);
[parameter, Ccalc, qcalc, AARDev] = fitIsotherm(C, q, [1 2 3], 'linearlangmuir', opt);
plot(Ccalc, qcalc, '-', 'LineWidth', 2);
[parameter, Ccalc, qcalc, AARDev] = fitIsotherm(C, q, [1 2], 'freundlich', opt);
plot(Ccalc, qcalc, '-', 'LineWidth', 2);
[parameter, Ccalc, qcalc, AARDev] = fitIsotherm(C, q, [1 2 3 4], 'bilangmuir', opt);
plot(Ccalc, qcalc, '-', 'LineWidth', 2);
[parameter, Ccalc, qcalc, AARDev] = fitIsotherm(C, q, [1 2 3], 'langmuirfreundlich', opt);
plot(Ccalc, qcalc, '-', 'LineWidth', 2);
[parameter, Ccalc, qcalc, AARDev] = fitIsotherm(C, q, [1 2 3], 'toth', opt);
plot(Ccalc, qcalc, '-', 'LineWidth', 2);
xlabel('\itC')
ylabel('\itq')
legend({'experimental', 'linear', 'langmuir', 'linearlangmuir', 'freundlich', 'bilangmuir', 'langmuirfreundlich', 'toth'}, 'Location','NorthWest')
