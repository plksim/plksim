#!/bin/bash
set -e

# Read build type, default value is debug
BUILD_TYPE=${1:-debug}
BUILD_TYPES=('debug' 'release' 'test')
if [[ ! ${BUILD_TYPES[@]} =~ "${BUILD_TYPE}" ]]; then
  BUILD_TYPE='debug'
fi

if [ "${BUILD_TYPE}" = 'release' ]; then
    BUILD_DIR="release"
    CMAKE_BUILD_TYPE="Release"
else
    BUILD_DIR="debug"
    CMAKE_BUILD_TYPE="Debug"
fi

# Use specified cmake if necessary
CMAKE=cmake
${CMAKE} . -B ${BUILD_DIR} -D CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
${CMAKE} --build ${BUILD_DIR}

# Run tests
if [ "${BUILD_TYPE}" = 'test' ]; then
    ctest --test-dir ${BUILD_DIR} --verbose
fi
