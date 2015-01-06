cmake_minimum_required(VERSION 2.6)
project(chipboard)

include(CheckCXXCompilerFlag) 
check_cxx_compiler_flag(-std=c++11 Cxx11_SUPPORTED)
if (NOT Cxx11_SUPPORTED)
  message(FATAL_ERROR "Compiler with C++11 support required!")
endif()

find_package(Qt4 REQUIRED)

include_directories(${QT_INCLUDES} ${CMAKE_CURRENT_BINARY_DIR})

set(chipboard_SRCS vplot.cpp linplot.cpp genome.cpp readcontainer.cpp utils.cpp)
set(CMAKE_CXX_FLAGS "-std=c++14 -g")

qt4_automoc(${chipboard_SRCS})
add_executable(chipboard main.cpp ${chipboard_SRCS})
target_link_libraries(chipboard ${QT_QTCORE_LIBRARY} ${QT_QTGUI_LIBRARY})

add_executable(chipboard-cli chipboard-cli.cpp ${chipboard_SRCS})
target_link_libraries(chipboard-cli ${QT_QTCORE_LIBRARY} ${QT_QTGUI_LIBRARY})

add_executable(testrun tests.cpp ${chipboard_SRCS})
target_link_libraries(testrun ${QT_QTCORE_LIBRARY} ${QT_QTGUI_LIBRARY})


#install(TARGETS chipboard RUNTIME DESTINATION bin)
