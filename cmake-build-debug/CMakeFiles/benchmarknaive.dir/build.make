# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.6

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /home/jervelund/Downloads/clion-2016.2.3/bin/cmake/bin/cmake

# The command to remove a file.
RM = /home/jervelund/Downloads/clion-2016.2.3/bin/cmake/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/jervelund/ClionProjects/DM818

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/jervelund/ClionProjects/DM818/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/benchmarknaive.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/benchmarknaive.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/benchmarknaive.dir/flags.make

CMakeFiles/benchmarknaive.dir/src/benchmark.cpp.o: CMakeFiles/benchmarknaive.dir/flags.make
CMakeFiles/benchmarknaive.dir/src/benchmark.cpp.o: ../src/benchmark.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jervelund/ClionProjects/DM818/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/benchmarknaive.dir/src/benchmark.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/benchmarknaive.dir/src/benchmark.cpp.o -c /home/jervelund/ClionProjects/DM818/src/benchmark.cpp

CMakeFiles/benchmarknaive.dir/src/benchmark.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/benchmarknaive.dir/src/benchmark.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jervelund/ClionProjects/DM818/src/benchmark.cpp > CMakeFiles/benchmarknaive.dir/src/benchmark.cpp.i

CMakeFiles/benchmarknaive.dir/src/benchmark.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/benchmarknaive.dir/src/benchmark.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jervelund/ClionProjects/DM818/src/benchmark.cpp -o CMakeFiles/benchmarknaive.dir/src/benchmark.cpp.s

CMakeFiles/benchmarknaive.dir/src/benchmark.cpp.o.requires:

.PHONY : CMakeFiles/benchmarknaive.dir/src/benchmark.cpp.o.requires

CMakeFiles/benchmarknaive.dir/src/benchmark.cpp.o.provides: CMakeFiles/benchmarknaive.dir/src/benchmark.cpp.o.requires
	$(MAKE) -f CMakeFiles/benchmarknaive.dir/build.make CMakeFiles/benchmarknaive.dir/src/benchmark.cpp.o.provides.build
.PHONY : CMakeFiles/benchmarknaive.dir/src/benchmark.cpp.o.provides

CMakeFiles/benchmarknaive.dir/src/benchmark.cpp.o.provides.build: CMakeFiles/benchmarknaive.dir/src/benchmark.cpp.o


CMakeFiles/benchmarknaive.dir/src/dgemm-naive.cpp.o: CMakeFiles/benchmarknaive.dir/flags.make
CMakeFiles/benchmarknaive.dir/src/dgemm-naive.cpp.o: ../src/dgemm-naive.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jervelund/ClionProjects/DM818/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/benchmarknaive.dir/src/dgemm-naive.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/benchmarknaive.dir/src/dgemm-naive.cpp.o -c /home/jervelund/ClionProjects/DM818/src/dgemm-naive.cpp

CMakeFiles/benchmarknaive.dir/src/dgemm-naive.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/benchmarknaive.dir/src/dgemm-naive.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jervelund/ClionProjects/DM818/src/dgemm-naive.cpp > CMakeFiles/benchmarknaive.dir/src/dgemm-naive.cpp.i

CMakeFiles/benchmarknaive.dir/src/dgemm-naive.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/benchmarknaive.dir/src/dgemm-naive.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jervelund/ClionProjects/DM818/src/dgemm-naive.cpp -o CMakeFiles/benchmarknaive.dir/src/dgemm-naive.cpp.s

CMakeFiles/benchmarknaive.dir/src/dgemm-naive.cpp.o.requires:

.PHONY : CMakeFiles/benchmarknaive.dir/src/dgemm-naive.cpp.o.requires

CMakeFiles/benchmarknaive.dir/src/dgemm-naive.cpp.o.provides: CMakeFiles/benchmarknaive.dir/src/dgemm-naive.cpp.o.requires
	$(MAKE) -f CMakeFiles/benchmarknaive.dir/build.make CMakeFiles/benchmarknaive.dir/src/dgemm-naive.cpp.o.provides.build
.PHONY : CMakeFiles/benchmarknaive.dir/src/dgemm-naive.cpp.o.provides

CMakeFiles/benchmarknaive.dir/src/dgemm-naive.cpp.o.provides.build: CMakeFiles/benchmarknaive.dir/src/dgemm-naive.cpp.o


# Object files for target benchmarknaive
benchmarknaive_OBJECTS = \
"CMakeFiles/benchmarknaive.dir/src/benchmark.cpp.o" \
"CMakeFiles/benchmarknaive.dir/src/dgemm-naive.cpp.o"

# External object files for target benchmarknaive
benchmarknaive_EXTERNAL_OBJECTS =

benchmarknaive: CMakeFiles/benchmarknaive.dir/src/benchmark.cpp.o
benchmarknaive: CMakeFiles/benchmarknaive.dir/src/dgemm-naive.cpp.o
benchmarknaive: CMakeFiles/benchmarknaive.dir/build.make
benchmarknaive: CMakeFiles/benchmarknaive.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/jervelund/ClionProjects/DM818/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable benchmarknaive"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/benchmarknaive.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/benchmarknaive.dir/build: benchmarknaive

.PHONY : CMakeFiles/benchmarknaive.dir/build

CMakeFiles/benchmarknaive.dir/requires: CMakeFiles/benchmarknaive.dir/src/benchmark.cpp.o.requires
CMakeFiles/benchmarknaive.dir/requires: CMakeFiles/benchmarknaive.dir/src/dgemm-naive.cpp.o.requires

.PHONY : CMakeFiles/benchmarknaive.dir/requires

CMakeFiles/benchmarknaive.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/benchmarknaive.dir/cmake_clean.cmake
.PHONY : CMakeFiles/benchmarknaive.dir/clean

CMakeFiles/benchmarknaive.dir/depend:
	cd /home/jervelund/ClionProjects/DM818/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/jervelund/ClionProjects/DM818 /home/jervelund/ClionProjects/DM818 /home/jervelund/ClionProjects/DM818/cmake-build-debug /home/jervelund/ClionProjects/DM818/cmake-build-debug /home/jervelund/ClionProjects/DM818/cmake-build-debug/CMakeFiles/benchmarknaive.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/benchmarknaive.dir/depend

