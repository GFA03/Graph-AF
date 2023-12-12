cmake_minimum_required(VERSION 3.4)
project(Graph-AF)
add_subdirectory(pybind11)
pybind11_add_module(graf graf.cpp)