# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/zzhfro/code/3DVision/3DReconstruciton/bundleadjustment

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/zzhfro/code/3DVision/3DReconstruciton/bundleadjustment/build

# Include any dependencies generated for this target.
include CMakeFiles/testba.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/testba.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/testba.dir/flags.make

CMakeFiles/testba.dir/batest.cpp.o: CMakeFiles/testba.dir/flags.make
CMakeFiles/testba.dir/batest.cpp.o: ../batest.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/zzhfro/code/3DVision/3DReconstruciton/bundleadjustment/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/testba.dir/batest.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/testba.dir/batest.cpp.o -c /home/zzhfro/code/3DVision/3DReconstruciton/bundleadjustment/batest.cpp

CMakeFiles/testba.dir/batest.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/testba.dir/batest.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/zzhfro/code/3DVision/3DReconstruciton/bundleadjustment/batest.cpp > CMakeFiles/testba.dir/batest.cpp.i

CMakeFiles/testba.dir/batest.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/testba.dir/batest.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/zzhfro/code/3DVision/3DReconstruciton/bundleadjustment/batest.cpp -o CMakeFiles/testba.dir/batest.cpp.s

# Object files for target testba
testba_OBJECTS = \
"CMakeFiles/testba.dir/batest.cpp.o"

# External object files for target testba
testba_EXTERNAL_OBJECTS =

testba: CMakeFiles/testba.dir/batest.cpp.o
testba: CMakeFiles/testba.dir/build.make
testba: /usr/local/lib/libopencv_dnn.so.3.4.16
testba: /usr/local/lib/libopencv_highgui.so.3.4.16
testba: /usr/local/lib/libopencv_ml.so.3.4.16
testba: /usr/local/lib/libopencv_objdetect.so.3.4.16
testba: /usr/local/lib/libopencv_shape.so.3.4.16
testba: /usr/local/lib/libopencv_stitching.so.3.4.16
testba: /usr/local/lib/libopencv_superres.so.3.4.16
testba: /usr/local/lib/libopencv_videostab.so.3.4.16
testba: /usr/local/lib/libopencv_viz.so.3.4.16
testba: /usr/local/lib/libopencv_calib3d.so.3.4.16
testba: /usr/local/lib/libopencv_features2d.so.3.4.16
testba: /usr/local/lib/libopencv_flann.so.3.4.16
testba: /usr/local/lib/libopencv_photo.so.3.4.16
testba: /usr/local/lib/libopencv_video.so.3.4.16
testba: /usr/local/lib/libopencv_videoio.so.3.4.16
testba: /usr/local/lib/libopencv_imgcodecs.so.3.4.16
testba: /usr/local/lib/libopencv_imgproc.so.3.4.16
testba: /usr/local/lib/libopencv_core.so.3.4.16
testba: CMakeFiles/testba.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/zzhfro/code/3DVision/3DReconstruciton/bundleadjustment/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable testba"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/testba.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/testba.dir/build: testba

.PHONY : CMakeFiles/testba.dir/build

CMakeFiles/testba.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/testba.dir/cmake_clean.cmake
.PHONY : CMakeFiles/testba.dir/clean

CMakeFiles/testba.dir/depend:
	cd /home/zzhfro/code/3DVision/3DReconstruciton/bundleadjustment/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/zzhfro/code/3DVision/3DReconstruciton/bundleadjustment /home/zzhfro/code/3DVision/3DReconstruciton/bundleadjustment /home/zzhfro/code/3DVision/3DReconstruciton/bundleadjustment/build /home/zzhfro/code/3DVision/3DReconstruciton/bundleadjustment/build /home/zzhfro/code/3DVision/3DReconstruciton/bundleadjustment/build/CMakeFiles/testba.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/testba.dir/depend
