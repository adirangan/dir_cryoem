function outputSignal = max_filter_1d_imdilate_0(inputSignal, filterDiameter)
    % inputSignal: 1-dimensional input signal
    % filterDiameter: diameter of the max-filter window
    
    % Ensure filterDiameter is odd for symmetric filtering
    if mod(filterDiameter, 2) == 0
        filterDiameter = filterDiameter + 1;
    end
    
    % Perform max-filter using imdilate
    filterSize = [1, filterDiameter];
    se = strel('arbitrary', ones(filterSize));
    outputSignal = imdilate(inputSignal, se);
end
%{

% Example usage:
% Define your input signal
inputSignal = [1, 3, 2, 7, 5, 8, 4, 6];

% Define filter diameter (window size)
filterDiameter = 3;

% Call the max_filter_1d_imdilate_0 function
outputSignal = max_filter_1d_imdilate_0(inputSignal, filterDiameter);

% Display the original and filtered signals
disp('Original Signal:');
disp(inputSignal);
disp('Filtered Signal:');
disp(outputSignal);
 %}
