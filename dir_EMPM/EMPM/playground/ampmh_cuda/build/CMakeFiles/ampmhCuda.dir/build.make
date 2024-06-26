# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.27

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /mnt/sw/nix/store/b8b9avcp9iqklkzm50a0bplgxnbshq3s-cmake-3.27.9/bin/cmake

# The command to remove a file.
RM = /mnt/sw/nix/store/b8b9avcp9iqklkzm50a0bplgxnbshq3s-cmake-3.27.9/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /mnt/home/wtang/Code/EMPM/ampmh_cuda

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /mnt/home/wtang/Code/EMPM/ampmh_cuda/build

# Include any dependencies generated for this target.
include CMakeFiles/ampmhCuda.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/ampmhCuda.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/ampmhCuda.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/ampmhCuda.dir/flags.make

CMakeFiles/ampmhCuda.dir/src/ampmhCuda.c.o: CMakeFiles/ampmhCuda.dir/flags.make
CMakeFiles/ampmhCuda.dir/src/ampmhCuda.c.o: /mnt/home/wtang/Code/EMPM/ampmh_cuda/src/ampmhCuda.c
CMakeFiles/ampmhCuda.dir/src/ampmhCuda.c.o: CMakeFiles/ampmhCuda.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/mnt/home/wtang/Code/EMPM/ampmh_cuda/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/ampmhCuda.dir/src/ampmhCuda.c.o"
	/mnt/sw/nix/store/l56mmgrmljf1di9ms77jah53db8a90bb-gcc-11.4.0/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/ampmhCuda.dir/src/ampmhCuda.c.o -MF CMakeFiles/ampmhCuda.dir/src/ampmhCuda.c.o.d -o CMakeFiles/ampmhCuda.dir/src/ampmhCuda.c.o -c /mnt/home/wtang/Code/EMPM/ampmh_cuda/src/ampmhCuda.c

CMakeFiles/ampmhCuda.dir/src/ampmhCuda.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing C source to CMakeFiles/ampmhCuda.dir/src/ampmhCuda.c.i"
	/mnt/sw/nix/store/l56mmgrmljf1di9ms77jah53db8a90bb-gcc-11.4.0/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /mnt/home/wtang/Code/EMPM/ampmh_cuda/src/ampmhCuda.c > CMakeFiles/ampmhCuda.dir/src/ampmhCuda.c.i

CMakeFiles/ampmhCuda.dir/src/ampmhCuda.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling C source to assembly CMakeFiles/ampmhCuda.dir/src/ampmhCuda.c.s"
	/mnt/sw/nix/store/l56mmgrmljf1di9ms77jah53db8a90bb-gcc-11.4.0/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /mnt/home/wtang/Code/EMPM/ampmh_cuda/src/ampmhCuda.c -o CMakeFiles/ampmhCuda.dir/src/ampmhCuda.c.s

CMakeFiles/ampmhCuda.dir/src/kernel/ampmhKernel.cu.o: CMakeFiles/ampmhCuda.dir/flags.make
CMakeFiles/ampmhCuda.dir/src/kernel/ampmhKernel.cu.o: CMakeFiles/ampmhCuda.dir/includes_CUDA.rsp
CMakeFiles/ampmhCuda.dir/src/kernel/ampmhKernel.cu.o: /mnt/home/wtang/Code/EMPM/ampmh_cuda/src/kernel/ampmhKernel.cu
CMakeFiles/ampmhCuda.dir/src/kernel/ampmhKernel.cu.o: CMakeFiles/ampmhCuda.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/mnt/home/wtang/Code/EMPM/ampmh_cuda/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CUDA object CMakeFiles/ampmhCuda.dir/src/kernel/ampmhKernel.cu.o"
	/mnt/sw/nix/store/blsyycc4ldkvdxzwyyv1xlk3qn9gblzd-cuda-12.2.1/bin/nvcc -forward-unknown-to-host-compiler $(CUDA_DEFINES) $(CUDA_INCLUDES) $(CUDA_FLAGS) -MD -MT CMakeFiles/ampmhCuda.dir/src/kernel/ampmhKernel.cu.o -MF CMakeFiles/ampmhCuda.dir/src/kernel/ampmhKernel.cu.o.d -x cu -c /mnt/home/wtang/Code/EMPM/ampmh_cuda/src/kernel/ampmhKernel.cu -o CMakeFiles/ampmhCuda.dir/src/kernel/ampmhKernel.cu.o

CMakeFiles/ampmhCuda.dir/src/kernel/ampmhKernel.cu.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CUDA source to CMakeFiles/ampmhCuda.dir/src/kernel/ampmhKernel.cu.i"
	$(CMAKE_COMMAND) -E cmake_unimplemented_variable CMAKE_CUDA_CREATE_PREPROCESSED_SOURCE

CMakeFiles/ampmhCuda.dir/src/kernel/ampmhKernel.cu.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CUDA source to assembly CMakeFiles/ampmhCuda.dir/src/kernel/ampmhKernel.cu.s"
	$(CMAKE_COMMAND) -E cmake_unimplemented_variable CMAKE_CUDA_CREATE_ASSEMBLY_SOURCE

# Object files for target ampmhCuda
ampmhCuda_OBJECTS = \
"CMakeFiles/ampmhCuda.dir/src/ampmhCuda.c.o" \
"CMakeFiles/ampmhCuda.dir/src/kernel/ampmhKernel.cu.o"

# External object files for target ampmhCuda
ampmhCuda_EXTERNAL_OBJECTS =

ampmhCuda: CMakeFiles/ampmhCuda.dir/src/ampmhCuda.c.o
ampmhCuda: CMakeFiles/ampmhCuda.dir/src/kernel/ampmhKernel.cu.o
ampmhCuda: CMakeFiles/ampmhCuda.dir/build.make
ampmhCuda: /mnt/sw/nix/store/blsyycc4ldkvdxzwyyv1xlk3qn9gblzd-cuda-12.2.1/lib64/libcudart.so
ampmhCuda: /usr/lib64/libcuda.so
ampmhCuda: /mnt/sw/nix/store/blsyycc4ldkvdxzwyyv1xlk3qn9gblzd-cuda-12.2.1/lib64/libculibos.a
ampmhCuda: /mnt/sw/nix/store/blsyycc4ldkvdxzwyyv1xlk3qn9gblzd-cuda-12.2.1/lib64/libcublas.so
ampmhCuda: /mnt/sw/nix/store/blsyycc4ldkvdxzwyyv1xlk3qn9gblzd-cuda-12.2.1/lib64/libcufftw.so
ampmhCuda: /mnt/sw/nix/store/blsyycc4ldkvdxzwyyv1xlk3qn9gblzd-cuda-12.2.1/lib64/libculibos.a
ampmhCuda: /mnt/sw/nix/store/blsyycc4ldkvdxzwyyv1xlk3qn9gblzd-cuda-12.2.1/lib64/libcublasLt.so
ampmhCuda: /mnt/sw/nix/store/blsyycc4ldkvdxzwyyv1xlk3qn9gblzd-cuda-12.2.1/lib64/libcufft.so
ampmhCuda: CMakeFiles/ampmhCuda.dir/linkLibs.rsp
ampmhCuda: CMakeFiles/ampmhCuda.dir/objects1.rsp
ampmhCuda: CMakeFiles/ampmhCuda.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/mnt/home/wtang/Code/EMPM/ampmh_cuda/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CUDA executable ampmhCuda"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/ampmhCuda.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/ampmhCuda.dir/build: ampmhCuda
.PHONY : CMakeFiles/ampmhCuda.dir/build

CMakeFiles/ampmhCuda.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/ampmhCuda.dir/cmake_clean.cmake
.PHONY : CMakeFiles/ampmhCuda.dir/clean

CMakeFiles/ampmhCuda.dir/depend:
	cd /mnt/home/wtang/Code/EMPM/ampmh_cuda/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /mnt/home/wtang/Code/EMPM/ampmh_cuda /mnt/home/wtang/Code/EMPM/ampmh_cuda /mnt/home/wtang/Code/EMPM/ampmh_cuda/build /mnt/home/wtang/Code/EMPM/ampmh_cuda/build /mnt/home/wtang/Code/EMPM/ampmh_cuda/build/CMakeFiles/ampmhCuda.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : CMakeFiles/ampmhCuda.dir/depend

