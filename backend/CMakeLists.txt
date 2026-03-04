cmake_minimum_required(VERSION 3.10)
project(rmobi CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include dirs
include_directories(${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/external)

# Single executable (no other main files!)
add_executable(rmobi
    src/main.cpp
    src/graph.cpp
    src/export_graph_json.cpp
    # add other non-main sources here:
    src/parser.cpp
    src/id_maps.cpp
    src/traversal.cpp
    src/config_loader.cpp
    src/export.cpp
    src/graph_utils.cpp
    src/analysis.cpp
)

# Optional warnings
if (MSVC)
  target_compile_options(rmobi PRIVATE /W4)
else()
  target_compile_options(rmobi PRIVATE -Wall -Wextra -Wpedantic)
endif()

