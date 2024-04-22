#!/bin/bash

. .bashrc

cfg_file="config.cfg"
BUILD_TYPE=$(grep 'build_type' "$cfg_file" | cut -d '=' -f2 | tr -d ' ')
REBUILD=$(grep 'rebuild_cmake' "$cfg_file" | cut -d '=' -f2 | tr -d ' ')

cd resource
python3 sqls.py
cd ..

cd build/$BUILD_TYPE/generators
if [ -f "CMakeCache.txt" ]; then
  [ "$REBUILD" = true ] && rm CMakeCache.txt && cmake -DCMAKE_BUILD_TYPE=$BUILD_TYPE -DCMAKE_TOOLCHAIN_FILE=conan_toolchain.cmake ../../..
else
  cmake -DCMAKE_BUILD_TYPE=$BUILD_TYPE -DCMAKE_TOOLCHAIN_FILE=conan_toolchain.cmake ../../..
fi
cmake --build . --config $BUILD_TYPE
./main