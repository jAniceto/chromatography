function gaussPeakFit_1p()
% gaussPeakFit_1p fits a gaussian curve to the data provided. Only one peak
% must be present.
% 
% Gaussian function: 
% f(x) = a * exp( -(x-b)^2 / (2*c^2) ) where,
% a is the max height of the peak,
% b is the position of the center of the peak (mean), and
% c controls the width of the peak (standard deviation)

clc


%% INPUT DATA (initial estimations)
a = 1;  % max height of the peak
b = 4.5;  % position of the center of the peak (mean)
c = 0.5;  % width of the peak (standard deviation)


%% INPUT PEAK DATA (matrix in the [x y] format)
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
4.890728916	0.090320894]; 


%% CALCULATIONS
data.x = xy(:,1);
data.y = xy(:,2);

% Parameter optimization using fminsearch
parameters = [a b c];
options = optimset('PlotFcns',@optimplotfval);
[parameters, fval, exitflag] = fminsearch(@fobj, parameters, options, data);

% Calculate Gaussian curve using optimized parameters
x = floor(min(data.x)):(ceil(max(data.x))-floor(min(data.x)))/100:ceil(max(data.x)); x = x';
a = parameters(1);
b = parameters(2);
c = parameters(3);
y = a .* exp( -(x-b).^2 ./ (2*c.^2) );

% Calculate quality of fit statistics
ycalc = a .* exp( -(data.x-b).^2 ./ (2*c.^2) );
AARD = sum(abs(data.y-ycalc)./data.y)/length(data.y);
[R,P] = corrcoef(data.y,ycalc);

% Calculate are under curve (integral)
Area = trapz(ycalc);

% Plot data
figure
plot(data.x, data.y, 'o', x, y, '-')

% Output results
fprintf('Optimized parameters: \nPeak height (a) = %.4f \nPeak mean (b) = %.4f \nPeak width (c) = %.4f \n\n', a, b, c)
fprintf('Fit quality: \nCorrelation coefficient (R) = %.4f , for a p-value of %.4f \nAARD = %.4f \n\n', R(1,2), P(1,2), AARD)
fprintf('Area under peak (A) = %.4f \n\n', Area)


%% Objective funtion
function f = fobj(parameters, data)

a = parameters(1);
b = parameters(2);
c = parameters(3);
ycalc = a .* exp( -(data.x-b).^2 ./ (2*c.^2) );

% Objective function for parameter optimization (Least-Squares)
f = sum( (data.y-ycalc).^2 );
