# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.24

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
CMAKE_COMMAND = /usr/local/bin/cmake

# The command to remove a file.
RM = /usr/local/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/yexw/maximalKPlex

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/yexw/maximalKPlex

# Include any dependencies generated for this target.
include src/graph/CMakeFiles/Graph.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include src/graph/CMakeFiles/Graph.dir/compiler_depend.make

# Include the progress variables for this target.
include src/graph/CMakeFiles/Graph.dir/progress.make

# Include the compile flags for this target's objects.
include src/graph/CMakeFiles/Graph.dir/flags.make

src/graph/CMakeFiles/Graph.dir/graph.cpp.o: src/graph/CMakeFiles/Graph.dir/flags.make
src/graph/CMakeFiles/Graph.dir/graph.cpp.o: src/graph/graph.cpp
src/graph/CMakeFiles/Graph.dir/graph.cpp.o: src/graph/CMakeFiles/Graph.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/yexw/maximalKPlex/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/graph/CMakeFiles/Graph.dir/graph.cpp.o"
	cd /home/yexw/maximalKPlex/src/graph && g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/graph/CMakeFiles/Graph.dir/graph.cpp.o -MF CMakeFiles/Graph.dir/graph.cpp.o.d -o CMakeFiles/Graph.dir/graph.cpp.o -c /home/yexw/maximalKPlex/src/graph/graph.cpp

src/graph/CMakeFiles/Graph.dir/graph.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Graph.dir/graph.cpp.i"
	cd /home/yexw/maximalKPlex/src/graph && g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/yexw/maximalKPlex/src/graph/graph.cpp > CMakeFiles/Graph.dir/graph.cpp.i

src/graph/CMakeFiles/Graph.dir/graph.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Graph.dir/graph.cpp.s"
	cd /home/yexw/maximalKPlex/src/graph && g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/yexw/maximalKPlex/src/graph/graph.cpp -o CMakeFiles/Graph.dir/graph.cpp.s

# Object files for target Graph
Graph_OBJECTS = \
"CMakeFiles/Graph.dir/graph.cpp.o"

# External object files for target Graph
Graph_EXTERNAL_OBJECTS =

src/graph/libGraph.a: src/graph/CMakeFiles/Graph.dir/graph.cpp.o
src/graph/libGraph.a: src/graph/CMakeFiles/Graph.dir/build.make
src/graph/libGraph.a: src/graph/CMakeFiles/Graph.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/yexw/maximalKPlex/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libGraph.a"
	cd /home/yexw/maximalKPlex/src/graph && $(CMAKE_COMMAND) -P CMakeFiles/Graph.dir/cmake_clean_target.cmake
	cd /home/yexw/maximalKPlex/src/graph && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/Graph.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/graph/CMakeFiles/Graph.dir/build: src/graph/libGraph.a
.PHONY : src/graph/CMakeFiles/Graph.dir/build

src/graph/CMakeFiles/Graph.dir/clean:
	cd /home/yexw/maximalKPlex/src/graph && $(CMAKE_COMMAND) -P CMakeFiles/Graph.dir/cmake_clean.cmake
.PHONY : src/graph/CMakeFiles/Graph.dir/clean

src/graph/CMakeFiles/Graph.dir/depend:
	cd /home/yexw/maximalKPlex && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/yexw/maximalKPlex /home/yexw/maximalKPlex/src/graph /home/yexw/maximalKPlex /home/yexw/maximalKPlex/src/graph /home/yexw/maximalKPlex/src/graph/CMakeFiles/Graph.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/graph/CMakeFiles/Graph.dir/depend

