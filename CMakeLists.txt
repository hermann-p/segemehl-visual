cmake_minimum_required(VERSION 2.6)
project(chipboard)

include(CheckCXXCompilerFlag) 
check_cxx_compiler_flag(-std=c++11 Cxx11_SUPPORTED)
if (NOT Cxx11_SUPPORTED)
  message(FATAL_ERROR "Compiler with C++11 support required!")
endif()

#find_package(Qt4 REQUIRED)

include_directories(${QT_INCLUDES} ${CMAKE_CURRENT_BINARY_DIR})

set(shared_sources genome.cpp readcontainer.cpp utils.cpp plotchromosome2.cpp vplot2.cpp linplot2.cpp)

set(CMAKE_CXX_FLAGS "-std=c++11 -g -fpermissive")

#qt4_automoc(${shared_sources})
add_library(shared_objects ${shared_sources})

#add_executable(chipboard main.cpp)
#target_link_libraries(chipboard ${QT_QTCORE_LIBRARY} ${QT_QTGUI_LIBRARY} shared_objects)

add_executable(chipboard-cli chipboard-cli.cpp)
target_link_libraries(chipboard-cli ${QT_QTCORE_LIBRARY} ${QT_QTGUI_LIBRARY} shared_objects)

#add_executable(testrun tests.cpp)
#target_link_libraries(testrun ${QT_QTCORE_LIBRARY} ${QT_QTGUI_LIBRARY} shared_objects)


#install(TARGETS chipboard RUNTIME DESTINATION bin)
