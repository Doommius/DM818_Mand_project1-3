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
include CMakeFiles/benchmarkblocked.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/benchmarkblocked.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/benchmarkblocked.dir/flags.make

CMakeFiles/benchmarkblocked.dir/src/benchmark.cpp.o: CMakeFiles/benchmarkblocked.dir/flags.make
CMakeFiles/benchmarkblocked.dir/src/benchmark.cpp.o: ../src/benchmark.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jervelund/ClionProjects/DM818/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/benchmarkblocked.dir/src/benchmark.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/benchmarkblocked.dir/src/benchmark.cpp.o -c /home/jervelund/ClionProjects/DM818/src/benchmark.cpp

CMakeFiles/benchmarkblocked.dir/src/benchmark.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/benchmarkblocked.dir/src/benchmark.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jervelund/ClionProjects/DM818/src/benchmark.cpp > CMakeFiles/benchmarkblocked.dir/src/benchmark.cpp.i

CMakeFiles/benchmarkblocked.dir/src/benchmark.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/benchmarkblocked.dir/src/benchmark.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jervelund/ClionProjects/DM818/src/benchmark.cpp -o CMakeFiles/benchmarkblocked.dir/src/benchmark.cpp.s

CMakeFiles/benchmarkblocked.dir/src/benchmark.cpp.o.requires:

.PHONY : CMakeFiles/benchmarkblocked.dir/src/benchmark.cpp.o.requires

CMakeFiles/benchmarkblocked.dir/src/benchmark.cpp.o.provides: CMakeFiles/benchmarkblocked.dir/src/benchmark.cpp.o.requires
	$(MAKE) -f CMakeFiles/benchmarkblocked.dir/build.make CMakeFiles/benchmarkblocked.dir/src/benchmark.cpp.o.provides.build
.PHONY : CMakeFiles/benchmarkblocked.dir/src/benchmark.cpp.o.provides

CMakeFiles/benchmarkblocked.dir/src/benchmark.cpp.o.provides.build: CMakeFiles/benchmarkblocked.dir/src/benchmark.cpp.o


CMakeFiles/benchmarkblocked.dir/src/dgemm-blocked.cpp.o: CMakeFiles/benchmarkblocked.dir/flags.make
CMakeFiles/benchmarkblocked.dir/src/dgemm-blocked.cpp.o: ../src/dgemm-blocked.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jervelund/ClionProjects/DM818/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/benchmarkblocked.dir/src/dgemm-blocked.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/benchmarkblocked.dir/src/dgemm-blocked.cpp.o -c /home/jervelund/ClionProjects/DM818/src/dgemm-blocked.cpp

CMakeFiles/benchmarkblocked.dir/src/dgemm-blocked.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/benchmarkblocked.dir/src/dgemm-blocked.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jervelund/ClionProjects/DM818/src/dgemm-blocked.cpp > CMakeFiles/benchmarkblocked.dir/src/dgemm-blocked.cpp.i

CMakeFiles/benchmarkblocked.dir/src/dgemm-blocked.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/benchmarkblocked.dir/src/dgemm-blocked.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jervelund/ClionProjects/DM818/src/dgemm-blocked.cpp -o CMakeFiles/benchmarkblocked.dir/src/dgemm-blocked.cpp.s

CMakeFiles/benchmarkblocked.dir/src/dgemm-blocked.cpp.o.requires:

.PHONY : CMakeFiles/benchmarkblocked.dir/src/dgemm-blocked.cpp.o.requires

CMakeFiles/benchmarkblocked.dir/src/dgemm-blocked.cpp.o.provides: CMakeFiles/benchmarkblocked.dir/src/dgemm-blocked.cpp.o.requires
	$(MAKE) -f CMakeFiles/benchmarkblocked.dir/build.make CMakeFiles/benchmarkblocked.dir/src/dgemm-blocked.cpp.o.provides.build
.PHONY : CMakeFiles/benchmarkblocked.dir/src/dgemm-blocked.cpp.o.provides

CMakeFiles/benchmarkblocked.dir/src/dgemm-blocked.cpp.o.provides.build: CMakeFiles/benchmarkblocked.dir/src/dgemm-blocked.cpp.o


# Object files for target benchmarkblocked
benchmarkblocked_OBJECTS = \
"CMakeFiles/benchmarkblocked.dir/src/benchmark.cpp.o" \
"CMakeFiles/benchmarkblocked.dir/src/dgemm-blocked.cpp.o"

# External object files for target benchmarkblocked
benchmarkblocked_EXTERNAL_OBJECTS =

benchmarkblocked: CMakeFiles/benchmarkblocked.dir/src/benchmark.cpp.o
benchmarkblocked: CMakeFiles/benchmarkblocked.dir/src/dgemm-blocked.cpp.o
benchmarkblocked: CMakeFiles/benchmarkblocked.dir/build.make
benchmarkblocked: CMakeFiles/benchmarkblocked.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/jervelund/ClionProjects/DM818/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable benchmarkblocked"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/benchmarkblocked.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/benchmarkblocked.dir/build: benchmarkblocked

.PHONY : CMakeFiles/benchmarkblocked.dir/build

CMakeFiles/benchmarkblocked.dir/requires: CMakeFiles/benchmarkblocked.dir/src/benchmark.cpp.o.requires
CMakeFiles/benchmarkblocked.dir/requires: CMakeFiles/benchmarkblocked.dir/src/dgemm-blocked.cpp.o.requires

.PHONY : CMakeFiles/benchmarkblocked.dir/requires

CMakeFiles/benchmarkblocked.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/benchmarkblocked.dir/cmake_clean.cmake
.PHONY : CMakeFiles/benchmarkblocked.dir/clean

CMakeFiles/benchmarkblocked.dir/depend:
	cd /home/jervelund/ClionProjects/DM818/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/jervelund/ClionProjects/DM818 /home/jervelund/ClionProjects/DM818 /home/jervelund/ClionProjects/DM818/cmake-build-debug /home/jervelund/ClionProjects/DM818/cmake-build-debug /home/jervelund/ClionProjects/DM818/cmake-build-debug/CMakeFiles/benchmarkblocked.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/benchmarkblocked.dir/depend

