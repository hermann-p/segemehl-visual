cmake_minimum_required(VERSION 2.6)
project(chipboard)
find_package(Qt4 REQUIRED)

include_directories(${QT_INCLUDES} ${CMAKE_CURRENT_BINARY_DIR})

set(chipboard_SRCS chipboard.cpp main.cpp)

qt4_automoc(${chipboard_SRCS})
add_executable(chipboard ${chipboard_SRCS})
target_link_libraries(chipboard ${QT_QTCORE_LIBRARY} ${QT_QTGUI_LIBRARY})

install(TARGETS chipboard RUNTIME DESTINATION bin)
