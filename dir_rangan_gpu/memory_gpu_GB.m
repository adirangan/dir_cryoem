function [memory_available_gpu_in_GB,memory_total_gpu_in_GB] = memory_gpu_GB(flag_verbose);
if nargin<1; flag_verbose=[]; end;
if isempty(flag_verbose); flag_verbose=1; end;
memory_available_gpu_in_GB = gpuDevice().AvailableMemory/1e9;
memory_total_gpu_in_GB = gpuDevice().TotalMemory/1e9;
memory_used_gpu_in_GB = memory_total_gpu_in_GB - memory_available_gpu_in_GB;
if (flag_verbose>0); disp(sprintf(' %% memory_used_gpu_in_GB: %0.2f/%0.2f memory_available_gpu_in_GB: %0.2f/%0.2f',memory_used_gpu_in_GB,memory_total_gpu_in_GB,memory_available_gpu_in_GB,memory_total_gpu_in_GB)); end;
