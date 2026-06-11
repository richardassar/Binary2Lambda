# Binary2Lambda - C++ build.
#
# The C++ implementation is a single dependency-free C++17 translation unit;
# this file configures and builds it with the standard CMake workflow so the
# include paths, standard level and optimisation flags are set for you:
#
#   cmake -B build -DCMAKE_BUILD_TYPE=Release
#   cmake --build build
#   ctest --test-dir build --output-on-failure
#   ./build/lambda_bijection_cpp                 # self-test + demo
#
# The Rust port builds with Cargo (rust/Cargo.toml); the Python reference
# needs no build step.
cmake_minimum_required(VERSION 3.16)

project(Binary2Lambda
  VERSION 1.0.0
  DESCRIPTION "Total bijection between binary strings and closed lambda terms"
  LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Default to an optimised build for single-config generators when the caller
# does not pick one; the bignum is constant-factor heavy without it.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    Debug Release RelWithDebInfo MinSizeRel)
endif()

add_executable(lambda_bijection_cpp cpp/lambda_bijection.cpp)

enable_testing()
# The binary runs its full self-test on the default (no-argument) path and
# exits non-zero on any failure, so ctest needs only to run it.
add_test(NAME self_test COMMAND lambda_bijection_cpp)
add_test(NAME compress_vectors COMMAND lambda_bijection_cpp --compress-vectors)
add_test(NAME decode_vectors COMMAND lambda_bijection_cpp --vectors)
