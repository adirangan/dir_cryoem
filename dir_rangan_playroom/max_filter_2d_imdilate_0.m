function outputSignal = max_filter_2d_imdilate_0(inputSignal, filterDiameter)
    % inputSignal: 2-dimensional input signal
    % filterDiameter: diameter of the max-filter window
    % Ensure filterDiameter is odd for symmetric filtering
    if mod(filterDiameter, 2) == 0; filterDiameter = filterDiameter + 1; end;
    outputSignal = imdilate(inputSignal, strel('square', filterDiameter));
%{
 
% Example usage:
% Define your input signal
inputSignal = round(128*rand(32,48));
%inputSignal = zeros(32,48); tmp_p_=randperm(numel(inputSignal),32); inputSignal(tmp_p_) = 1;
% Define filter diameter (window size)
filterDiameter = 3;
% Call the max_filter_1d_imdilate_0 function
outputSignal = max_filter_2d_imdilate_0(inputSignal, filterDiameter);
% Display the original and filtered signals
figure(1);clf;set(gcf,'Position',1+[0,0,1024,512]);fig81s;
subplot(1,2,1);imagesc(inputSignal);title('Original Signal:'); axis image; axisnotick;
subplot(1,2,2);imagesc(outputSignal);title('Filtered Signal:'); axis image; axisnotick;

 %}
