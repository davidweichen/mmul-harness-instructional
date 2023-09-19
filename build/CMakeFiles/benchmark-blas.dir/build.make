# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/gator/mmul-harness-instructional

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/gator/mmul-harness-instructional/build

# Include any dependencies generated for this target.
include CMakeFiles/benchmark-blas.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/benchmark-blas.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/benchmark-blas.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/benchmark-blas.dir/flags.make

CMakeFiles/benchmark-blas.dir/dgemm-blas.cpp.o: CMakeFiles/benchmark-blas.dir/flags.make
CMakeFiles/benchmark-blas.dir/dgemm-blas.cpp.o: ../dgemm-blas.cpp
CMakeFiles/benchmark-blas.dir/dgemm-blas.cpp.o: CMakeFiles/benchmark-blas.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/gator/mmul-harness-instructional/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/benchmark-blas.dir/dgemm-blas.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/benchmark-blas.dir/dgemm-blas.cpp.o -MF CMakeFiles/benchmark-blas.dir/dgemm-blas.cpp.o.d -o CMakeFiles/benchmark-blas.dir/dgemm-blas.cpp.o -c /home/gator/mmul-harness-instructional/dgemm-blas.cpp

CMakeFiles/benchmark-blas.dir/dgemm-blas.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/benchmark-blas.dir/dgemm-blas.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/gator/mmul-harness-instructional/dgemm-blas.cpp > CMakeFiles/benchmark-blas.dir/dgemm-blas.cpp.i

CMakeFiles/benchmark-blas.dir/dgemm-blas.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/benchmark-blas.dir/dgemm-blas.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/gator/mmul-harness-instructional/dgemm-blas.cpp -o CMakeFiles/benchmark-blas.dir/dgemm-blas.cpp.s

# Object files for target benchmark-blas
benchmark__blas_OBJECTS = \
"CMakeFiles/benchmark-blas.dir/dgemm-blas.cpp.o"

# External object files for target benchmark-blas
benchmark__blas_EXTERNAL_OBJECTS = \
"/home/gator/mmul-harness-instructional/build/CMakeFiles/benchmark.dir/benchmark.cpp.o"

benchmark-blas: CMakeFiles/benchmark-blas.dir/dgemm-blas.cpp.o
benchmark-blas: CMakeFiles/benchmark.dir/benchmark.cpp.o
benchmark-blas: CMakeFiles/benchmark-blas.dir/build.make
benchmark-blas: /usr/lib/x86_64-linux-gnu/libblas.so
benchmark-blas: CMakeFiles/benchmark-blas.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/gator/mmul-harness-instructional/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable benchmark-blas"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/benchmark-blas.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/benchmark-blas.dir/build: benchmark-blas
.PHONY : CMakeFiles/benchmark-blas.dir/build

CMakeFiles/benchmark-blas.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/benchmark-blas.dir/cmake_clean.cmake
.PHONY : CMakeFiles/benchmark-blas.dir/clean

CMakeFiles/benchmark-blas.dir/depend:
	cd /home/gator/mmul-harness-instructional/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/gator/mmul-harness-instructional /home/gator/mmul-harness-instructional /home/gator/mmul-harness-instructional/build /home/gator/mmul-harness-instructional/build /home/gator/mmul-harness-instructional/build/CMakeFiles/benchmark-blas.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/benchmark-blas.dir/depend

