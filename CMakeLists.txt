cmake_minimum_required(VERSION 2.8.4)
project(Protein)

#ExternalProject_Add(GECODE gecode)

find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)

include_directories( ${OPENGL_INCLUDE_DIRS}  ${GLUT_INCLUDE_DIRS} )

#include_directories(BIU)
#include_directories(BIU/util)
#include_directories(BIU/qrng)
#include_directories(gecode/gecode)
#include_directories(CPSP)
#include_directories(CPSP/gecode)
include_directories(GL)


#add_subdirectory(BIU)
#add_subdirectory(BIU/util)
#add_subdirectory(BIU/qrng)
#add_subdirectory(gecode)
#add_subdirectory(CPSP)
#add_subdirectory(CPSP/gecode)
add_subdirectory(GL)


#link_directories(BIU/util)
#link_directories(BIU/qrng)
#link_directories(BIU)
#link_directories(gecode)
#link_directories(CPSP/gecode)
#link_directories(CPSP)
link_directories(GL)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(SOURCE_FILES main.cpp cLatticeModel.cpp cMonteCarloMover.cpp)
add_executable(Protein ${SOURCE_FILES})

#target_link_libraries(Protein GLFRAMEWORK CPSP GECODE_INTERFACE gecodeminimodel BIU BIU_UTILS BIU_QRNG ${OPENGL_LIBRARIES} ${GLUT_LIBRARY} )
target_link_libraries(Protein GLFRAMEWORK ${OPENGL_LIBRARIES} ${GLUT_LIBRARY} )
