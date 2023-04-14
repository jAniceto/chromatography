function gaussPeakFit(xy, peakSplit, a, b, c)
% GAUSSPEAKFIT  Fits a gaussian curve to the data provided. Works for any
% number of peaks.
%
% Gaussian function : f(x) = a * exp( -(x-b)^2 / (2*c^2) ) where,
% a is the max height of the peak,
% b is the position of the center of the peak (mean), and
% c controls the width of the peak (standard deviation)
% 
% INPUTS
% ======
% xy : matrix in the [x y] format containing x, y data points
% 
% peakSplit : vector containing the x points around which the data is
% separated into ist respective peaks. For 2 peaks, peakSplit must contain
% 1 element. To fit a single peak peakSplit must be an empty vector.
%
% a : vector cointainig the initial estimations for the max height of the
% peak (a) for each peak.
%
% b : vector cointainig the initial estimations for the mean of the peak 
% (b) for each peak.
%
% c : vector cointainig the initial estimations for the standard deviation 
% of the peak (c) for each peak. Note: a good way to provide a estimate for
% c is to give the value of peak width at the baseline divided by 6.
%
% OUTPUTS
% =======
% Exitflag : optimization status. If exitflag = 1, optimization converged
%
% Peak height (a),  Peak mean (b), and Peak standard deviation (c) : optimizaed
% parameters
%
% R : correlation coefficient calculated using Matlab corrcoef function
%
% AARD : absolute average relative deviation
%
% A : Area under peak calculated using Matlab trapz function
%
% Figure containing the original and fited data


%% Default arguments (Example)
if nargin == 0
    peakSplit = [4.9];    % x point around which data is split for the peaks
    a = [1 1];           % max height of the peak
    b = [4.5 5.5];         % position of the center of the peak (mean)
    c = [0.5 0.5];         % standard deviation (width of the peak)

    % Peak data points (matrix in the [x y] format)
    xy = [3.400862695	0.000621997
    3.556320916	0.000770183
    3.730181839	0.000121869
    3.850869582	0.000603344
    3.941791752	0.001993669
    4.006435059	0.013544969
    4.047465804	0.027979439
    4.085800072	0.046092345
    4.113247889	0.066926177
    4.129714733	0.079206249
    4.154501394	0.107931563
    4.171020316	0.126423484
    4.179338628	0.142689194
    4.1918213	0.167708429
    4.204281423	0.190038086
    4.229034169	0.214717919
    4.245609834	0.239978242
    4.259499693	0.270189692
    4.274720013	0.296441438
    4.291354831	0.328757416
    4.298334058	0.347954171
    4.320442784	0.382569505
    4.332920804	0.407033897
    4.351532463	0.431161879
    4.366060997	0.456225998
    4.382634941	0.481280855
    4.403271092	0.502904633
    4.422863894	0.519646799
    4.458171398	0.531125828
    4.501642088	0.525286441
    4.526052464	0.509128451
    4.550407652	0.486387578
    4.56661501	0.467716124
    4.59324601	0.442825607
    4.607096851	0.416677396
    4.623238413	0.390157775
    4.635333725	0.368972856
    4.655534891	0.338711535
    4.667683792	0.323918569
    4.67980462	0.305777165
    4.695947558	0.279421689
    4.712140638	0.259047234
    4.724248422	0.239349876
    4.740439525	0.218739463
    4.752587073	0.203785201
    4.768810199	0.186994571
    4.785023606	0.169044669
    4.801203069	0.147045866
    4.817467823	0.135220612
    4.83363705	0.112000983
    4.890728916	0.090320894
    4.964297835	0.082138625
    5.005342193	0.098196874
    5.040955612	0.117081108
    5.058807016	0.131854364
    5.083555672	0.156046465
    5.106960837	0.182643978
    5.128976478	0.206156179
    5.162013166	0.243002094
    5.182678235	0.26807519
    5.203316728	0.289978242
    5.211601626	0.302258314
    5.236322462	0.323132042
    5.269277022	0.350181713
    5.285792432	0.368254723
    5.328652998	0.39181974
    5.382592094	0.39602649
    5.441358785	0.381300803
    5.480049838	0.360642515
    5.511971302	0.345793123
    5.54728124	0.328478142
    5.567617032	0.31427507
    5.601495748	0.288901043
    5.62457257	0.276333718
    5.636734108	0.26304826
    5.663008959	0.24415011
    5.685511212	0.225525818
    5.742465745	0.187464777
    5.775006893	0.165202659
    5.803506737	0.148865136
    5.823841191	0.134502479
    5.878181525	0.109933775
    5.908034763	0.092372018
    5.962380386	0.068434256
    6.038636529	0.049006623
    6.15416043	0.028697572
    6.314992975	0.012803532
    6.521171177	0.005739514
    6.583660977	0.003149451
    6.721103882	0.001504274
    6.780012519	0.00134013
    6.898641279	0.000211635
    7.010470242	0.001288834
    7.0827441	0.001288834
    7.246380684	0.001009561
    7.385473101	0.000929768
    7.502352887	0.000314226
    7.683086881	0.000391168]; 
end


%% CALCULATIONS
% Data treatment
Npeaks = length(peakSplit)+1; % number of peaks 

if Npeaks == 1
    data.x{1} = xy(:,1);
    data.y{1} = xy(:,2);
else
    % If there are more than one peak, split data points acording to respective peaks using peakSplit input
    n=0;
    for i = 1:length(peakSplit)
        for j = 1:length(xy(:,1))
            if xy(n+j,1) < peakSplit(i)
                data.x{i}(j,1) = xy(n+j,1);
                data.y{i}(j,1) = xy(n+j,2);
            else
                n=n+j-1;
                break
            end
        end
    end
    data.x{i+1} = xy(n+1:end,1);
    data.y{i+1} = xy(n+1:end,2);
end

% Parameter optimization using fminsearch
for i = 1:Npeaks
    
    parameters = [a(i) b(i) c(i)];
    options = optimset('PlotFcns',@optimplotfval);
    [parameters, fval, exitflag] = fminsearch(@fobj, parameters, options, [data.x{i} data.y{i}]);

    % Calculate Gaussian curve using optimized parameters
    x = floor(min(xy(:,1))):(ceil(max(xy(:,1)))-floor(min(xy(:,1))))/100:ceil(max(xy(:,1))); x = x';
    aOpt = parameters(1);
    bOpt = parameters(2);
    cOpt = parameters(3);
    y = aOpt .* exp( -(x-bOpt).^2 ./ (2*cOpt.^2) );

    % Calculate quality of fit statistics
    yCalc = aOpt .* exp( -(data.x{i}-bOpt).^2 ./ (2*cOpt.^2) );
    AARD = sum(abs(data.y{i}-yCalc)./data.y{i})/length(data.y{i});
    [R,P] = corrcoef(data.y{i},yCalc);

    % Calculate are under curve (integral)
    Area = trapz(yCalc);

    % Plot data
    figure(100)
    plot(data.x{i}, data.y{i}, 'o', x, y, '-')
    hold on

    % Output results
    fprintf('Optimized parameters for peak %i: \nExitflag = %i \nPeak height (a) = %.4f \nPeak mean (b) = %.4f \nPeak standard deviation (c) = %.4f \n\n', i, exitflag, aOpt, bOpt, cOpt)
    fprintf('Fit quality for peak %i: \nCorrelation coefficient (R) = %.4f , for a p-value of %.4f \nAARD = %.4f \n\n', i, R(1,2), P(1,2), AARD)
    fprintf('Area under peak %i (A) = %.4f \n\n', i, Area)
    if i ~= Npeaks
        fprintf('----- \n\n')
    end

end
hold off


%% Objective funtion
function f = fobj(parameters, data)

xExp = data(:,1);
yExp = data(:,2);

a = parameters(1);
b = parameters(2);
c = parameters(3);
yCalc = a .* exp( -(xExp-b).^2 ./ (2*c.^2) );

% Objective function for parameter optimization (Least-Squares)
f = sum( (yExp-yCalc).^2 );
