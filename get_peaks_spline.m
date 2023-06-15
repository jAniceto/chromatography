function get_peaks_spline(data, prominence, verbose)
% Fit Gaussian curves to chromatogram signal.
%
% Gaussian function : f(x) = a * exp( -(x-b)^2 / (2*c^2) ) where,
%   a is the max height of the peak,
%   b is the position of the center of the peak (mean), and
%   c controls the width of the peak (standard deviation)
%
% INPUTS
% ======
% data : (n x 2) array
%   Array with two columns containing time vs signal data.
%
% prominence : float
%   Threshold value of peak height to be considered a peak. If too many 
%   peaks are found, increase the value. If no peaks are found, try to 
%   decrease the value. Default is 200.
%
% verbose : logical
%   Wether or not to show all fitting process data. Default is false.
% 
% EXAMPLE
% =======
% get_peaks(original_data)
% get_peaks(original_data, 300, true)


% Argument handling (requires MATLAB 2022+)
arguments
    data (:,2) double;
    prominence (1,1) double = 200
    verbose (1,1) logical = false
end

fprintf('Processing signal... \n\n')

x = data(:,1);  % time
y = data(:,2);  % absorbance

% Baseline fix (minimum to zero)
ymin = min(y);
y = y - ymin;

% Plot chromatogram
figure(1)
plot(x, y, '-', 'DisplayName','Signal')
xlabel('{\itt} (min)')
ylabel('{\itA} (mAU)')

% Find peaks
[pks, locs, w, p] = findpeaks(y, x, 'MinPeakProminence', prominence);
n_peaks = length(pks);  % number of peaks
text(locs+0.2, pks, num2str((1:numel(pks))'))  % label peaks in the figure

% Terminate is no peaks are found
if n_peaks < 1
    fprintf('NO PEAKS FOUND! Try changing the prominence value.\n\n')
    return;
end

% Create different sets of data (one for each peak)
if n_peaks == 1
    split_data{1} = data;
else
    % Find peak borders (using inverted the signal)
    [border, border_locs] = findpeaks(-y, x, 'MinPeakProminence', prominence);
    
    if verbose
        figure(1)
        hold on;
        plot([border_locs border_locs], [0 max(y)], '-.', 'Color','#808080', 'DisplayName','Split')
    end

    % Split data points acording to respective peaks using border_locs
    n = 0;
    for i = 1:length(border_locs)
        for j = 1:length(x)
            if x(n+j) < border_locs(i)
                split_data{i}(j,1) = x(n+j);
                split_data{i}(j,2) = y(n+j);
            else
                n = n+j-1;
                break
            end
        end
    end
    split_data{i+1}(:,1) = x(n+1:end);
    split_data{i+1}(:,2) = y(n+1:end);
end

% Fit Gaussian peaks to data
for i = 1:n_peaks
    
    x = split_data{i}(:,1);
    y = split_data{i}(:,2);

    % Fit spline
    [curve, goodness, output] = fit(x, y, 'smoothingspline');
    exitflag = output.exitflag;
        
    % Calculate are under curve (integral)
    xspan = linspace(x(1), y(end), 1000)';
    int = integrate(curve, xspan, x(1));
    area = int(end);
     
    % Plot fitted spline
    figure(1)
    hold on;
    peak_name = sprintf('Peak %i', i);
%     plot(x_calc, y_calc, '--', 'DisplayName',peak_name)
    plot(curve, x, y);
    legend;

    % Calculate quality of fit statistics
    R2 = goodness.rsquare;
    R2adj = goodness.adjrsquare;

    % Output results
    fprintf('PEAK %i \n', i)
    fprintf('exitflag = %i \n', exitflag)
    fprintf('Area under fitted peak (A) = %.4f \n', area)
    % fprintf('AARD = %.4f \n', AARD)
    fprintf('Correlation coefficient (R^2) = %.4f \n', R2)
    fprintf('Adjusted correlation coefficient (R^2adj) = %.4f \n\n', R2adj)
    if i ~= n_peaks
        fprintf('----- \n\n')
    end

end



function f = fobj(parameters, data)
% Objective function for parameter optimization (Least-Squares)

x_exp = data(:,1);
y_exp = data(:,2);

a = parameters(1);
b = parameters(2);
c = parameters(3);

y_calc = a .* exp( -(x_exp-b).^2 ./ (2*c.^2) );

f = sum( (y_exp - y_calc).^2 );
