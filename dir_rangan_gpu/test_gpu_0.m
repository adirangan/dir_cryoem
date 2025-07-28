%%%%%%%%;
% https://www.mathworks.com/help/parallel-computing/measuring-gpu-performance.html ;
%%%%%%%%;

%%%%%%%%;
% check gpu setup. ;
%%%%%%%%;
gpu = gpuDevice;
disp(gpu.Name + " GPU detected and available.");

%%%%%%%%;
% Measure Host/GPU Bandwidth. ;
%%%%%%%%;
sizeOfDouble = 8;
maxSize = 0.25*gpu.AvailableMemory;
maxNumTests = 15;
sizes = logspace(4,log10(maxSize),maxNumTests);
sizes(sizes/sizeOfDouble > intmax) = [];

numTests = numel(sizes);
numElements = floor(sizes/sizeOfDouble);
sendTimes = inf(1,numTests);
gatherTimes = inf(1,numTests);

for idx=1:numTests
    disp("Test " + idx + " of " + numTests + ". Timing send and gather for array with " + numElements(idx) + " elements.")

    % Generate random data on GPU and host.
    gpuData = randi([0 9],numElements(idx),1,"gpuArray");
    hostData = gather(gpuData);

    % Time sending data to GPU.
    sendFcn = @() gpuArray(hostData);
    sendTimes(idx) = gputimeit(sendFcn);

    % Time gathering data back from GPU.
    gatherFcn = @() gather(gpuData);
    gatherTimes(idx) = gputimeit(gatherFcn);
end

sendBandwidth = (sizes./sendTimes)/1e9;
gatherBandwidth = (sizes./gatherTimes)/1e9;

[maxSendBandwidth,maxSendIdx] = max(sendBandwidth);
[maxGatherBandwidth,maxGatherIdx] = max(gatherBandwidth);
fprintf("Achieved peak send speed of %.2f GB/s",maxSendBandwidth)

fprintf("Achieved peak gather speed of %.2f GB/s",maxGatherBandwidth)

figure
semilogx(sizes,sendBandwidth,MarkerIndices=maxSendIdx,Marker="o")
hold on
semilogx(sizes,gatherBandwidth,MarkerIndices=maxGatherIdx,Marker="o")
grid on
title("Data Transfer Bandwidth")
xlabel("Array size (bytes)")
ylabel("Transfer speed (GB/s)")
legend(["Send to GPU" "Gather from GPU"],Location="SouthEast")
hold off

%%%%%%%%;
% Measure Read and Write Speed During Memory Intensive Operations. ;
%%%%%%%%;

reset(gpu);
sizes = logspace(4.5,log10(maxSize),maxNumTests);
sizes(sizes/sizeOfDouble > intmax) = [];
numTests = numel(sizes);
numElements = floor(sizes/sizeOfDouble);
memoryTimesGPU = inf(1,numTests);
memoryTimesHost = inf(1,numTests);

for idx=1:numTests
    disp("Test " + idx + " of " + numTests + ". Timing plus operation on GPU and CPU for arrays with " + numElements(idx) + " elements.")

    % Generate random data on GPU and host.
    gpuData = randi([0 9],numElements(idx),1,"gpuArray");
    hostData = gather(gpuData);

    % Time the plus function on GPU.
    plusFcn = @() plus(gpuData,1.0);
    memoryTimesGPU(idx) = gputimeit(plusFcn);

    % Time the plus function on host.
    plusFcn = @() plus(hostData,1.0);
    memoryTimesHost(idx) = timeit(plusFcn);
end

memoryBandwidthGPU = 2*(sizes./memoryTimesGPU)/1e9;
memoryBandwidthHost = 2*(sizes./memoryTimesHost)/1e9;

[maxBWGPU,maxBWIdxGPU] = max(memoryBandwidthGPU);
[maxBWHost,maxBWIdxHost] = max(memoryBandwidthHost);
fprintf("Achieved peak read+write speed on the GPU: %.2f GB/s",maxBWGPU)

fprintf("Achieved peak read+write speed on the host: %.2f GB/s",maxBWHost)

figure
semilogx(sizes,memoryBandwidthGPU,MarkerIndices=maxBWIdxGPU,Marker="o")
hold on
semilogx(sizes,memoryBandwidthHost,MarkerIndices=maxBWIdxHost,Marker="o")

grid on
title("Read+Write Bandwidth")
xlabel("Array size (bytes)")
ylabel("Speed (GB/s)")
legend(["GPU" "Host"],Location="NorthWest")
hold off

%%%%%%%%;
% Measure Processing Power During Computationally Intensive Operations. ;
%%%%%%%%;

reset(gpu)
sizes = logspace(4,log10(maxSize)-1,maxNumTests);
sizes(sizes/sizeOfDouble > intmax) = [];

gpu.SingleDoubleRatio

numTests = numel(sizes);
NDouble = floor(sqrt(sizes/sizeOfDouble));
mmTimesHostDouble = inf(1,numTests);
mmTimesGPUDouble = inf(1,numTests);

for idx=1:numTests
    disp("Test " + idx + " of " + numTests + ". Timing double-precision matrix-matrix multiplication with " + NDouble(idx)^2 + " elements.")
    
    % Generate random data on GPU.
    A = rand(NDouble(idx),"gpuArray");
    B = rand(NDouble(idx),"gpuArray");
    
    % Time the matrix multiplication on GPU.
    mmTimesGPUDouble(idx) = gputimeit(@() A*B);

    % Gather the data and time matrix multiplication on the host.
    A = gather(A);
    B = gather(B);
    mmTimesHostDouble(idx) = timeit(@() A*B);
end

mmFlopsHostDouble = (2*NDouble.^3 - NDouble.^2)./mmTimesHostDouble;
[maxFlopsHostDouble,maxFlopsHostDoubleIdx] = max(mmFlopsHostDouble);
mmFlopsGPUDouble = (2*NDouble.^3 - NDouble.^2)./mmTimesGPUDouble;
[maxFlopsGPUDouble,maxFlopsGPUDoubleIdx] = max(mmFlopsGPUDouble);
fprintf("Achieved peak double-precision processing power on the GPU: %.2f TFLOPS",maxFlopsGPUDouble/1e12)

fprintf("Achieved peak double-precision processing power on the host: %.2f TFLOPS",maxFlopsHostDouble/1e12)

NSingle = floor(sqrt(sizes/(sizeOfDouble/2)));
mmTimesHostSingle = inf(1,numTests);
mmTimesGPUSingle = inf(1,numTests);

for idx=1:numTests
    disp("Test " + idx + " of " + numTests + ". Timing single-precision matrix-matrix multiplication with " + NSingle(idx)^2 + " elements.")

    % Generate random, single-precision data on GPU.
    A = rand(NSingle(idx),"single","gpuArray");
    B = rand(NSingle(idx),"single","gpuArray");

    % Time matrix multiplication on GPU.
    mmTimesGPUSingle(idx) = gputimeit(@() A*B);

    % Gather the data and time matrix multiplication on the host.
    A = gather(A);
    B = gather(B);
    mmTimesHostSingle(idx) = timeit(@() A*B);
end

mmFlopsHostSingle = (2*NSingle.^3 - NSingle.^2)./mmTimesHostSingle;
[maxFlopsHostSingle,maxFlopsHostSingleIdx] = max(mmFlopsHostSingle);
mmFlopsGPUSingle = (2*NSingle.^3 - NSingle.^2)./mmTimesGPUSingle;
[maxFlopsGPUSingle,maxFlopsGPUSingleIdx] = max(mmFlopsGPUSingle);
fprintf("Achieved peak single-precision processing power on the GPU: %.2f TFLOPS",maxFlopsGPUSingle/1e12)

fprintf("Achieved peak single-precision processing power on the host: %.2f TFLOPS",maxFlopsHostSingle/1e12)

figure
loglog(NDouble.^2,mmFlopsGPUDouble,MarkerIndices=maxFlopsGPUDoubleIdx,Marker="o")
hold on
loglog(NSingle.^2,mmFlopsGPUSingle,MarkerIndices=maxFlopsGPUSingleIdx,Marker="o")
loglog(NDouble.^2,mmFlopsHostDouble,MarkerIndices=maxFlopsHostDoubleIdx,Marker="o")
loglog(NSingle.^2,mmFlopsHostSingle,MarkerIndices=maxFlopsHostSingleIdx,Marker="o")

grid on
title("Matrix-Matrix Multiply")
xlabel("Matrix size (numel)")
ylabel("processing power (FLOPS)")
legend(["GPU double" "GPU single" "Host double" "Host single"],Location="SouthEast")
hold off

