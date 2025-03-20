project(spg)

# There is no real support for .obj or .msh meshes in any of the demos yet, so disable for now
set(USE_TINYOBJLOADER TRUE)
set(USE_MSHIO FALSE)


add_library(${PROJECT_NAME})

file(GLOB SRCF
source/sim/simObject/*.cpp
source/sim/solver/*.cpp
source/sim/energy/*.cpp
source/geom/io/*.cpp
source/utils/*.cpp
)

file(GLOB INCF
include/spg/sim/simObject/*.h
include/spg/sim/energy/*.h
include/spg/sim/solver/*.h
include/spg/geom/*.h
include/spg/geom/io/*.h
include/spg/utils/*.h
include/spg/*.h
)

add_subdirectory(../external/eigen/ ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/eigen)
add_subdirectory(../external/TinyAD/ ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/TinyAD)

if(USE_TINYOBJLOADER)
  add_subdirectory(../external/tinyobjloader/ ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/tinyobjloader)
  add_compile_definitions(TINYOBJLOADER_ACTIVE)
endif()

if(USE_MSHIO)
  add_subdirectory(../external/MshIO/ ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/MshIO)
  add_compile_definitions(MSHIO_ACTIVE)
endif()

target_include_directories(
  ${PROJECT_NAME}
  PUBLIC
  ../external/TinyAD/include
  include
)

# Allow for big compilation objects creation
if (MSVC)
  target_compile_options(${PROJECT_NAME} PUBLIC /bigobj)
endif ()

target_sources(${PROJECT_NAME} PRIVATE ${SRCF} ${INCF})

target_link_libraries(${PROJECT_NAME} PUBLIC Eigen3::Eigen OpenMP::OpenMP_CXX)

if(USE_TINYOBJLOADER)
  target_link_libraries(${PROJECT_NAME} PUBLIC tinyobjloader)
endif()

if(USE_MSHIO)
  target_link_libraries(${PROJECT_NAME} PUBLIC mshio::mshio)
endif()