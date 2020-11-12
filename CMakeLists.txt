cmake_minimum_required(VERSION 2.8)
project(tiebrush)

set (CMAKE_STATIC_LINKER_FLAGS "-Wl,--as-needed.-lcurl")

set(CMAKE_CXX_STANDARD 11)

include_directories("${CMAKE_SOURCE_DIR}/include/gclib/")
include_directories("${CMAKE_SOURCE_DIR}/include/htslib/")
link_directories("${CMAKE_SOURCE_DIR}/include/htslib")

include(ExternalProject)
ExternalProject_Add(htslib
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/include/htslib/
    BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND autoheader COMMAND autoconf COMMAND ./configure --without-libdeflate
    BUILD_COMMAND ${MAKE}
    INSTALL_COMMAND ""
)

add_executable(tiebrush src/tiebrush.cpp
                        src/commons.h
                        src/GSam.cpp
                        src/GSam.h
                        src/tmerge.cpp
                        src/tmerge.h
                        ${CMAKE_SOURCE_DIR}/include/gclib/GStr.h
                        ${CMAKE_SOURCE_DIR}/include/gclib/GStr.cpp
                        ${CMAKE_SOURCE_DIR}/include/gclib/GArgs.h
                        ${CMAKE_SOURCE_DIR}/include/gclib/GArgs.cpp
                        ${CMAKE_SOURCE_DIR}/include/gclib/GBase.h
                        ${CMAKE_SOURCE_DIR}/include/gclib/GBase.cpp)

add_dependencies(tiebrush htslib)
target_link_libraries(tiebrush -lhts -lz -lbz2 -llzma -pthread -lcrypto -lcurl)

add_executable(tiecov src/tiecov.cpp
                    src/commons.h
                    src/GSam.cpp
                    src/GSam.h
                    src/tmerge.cpp
                    src/tmerge.h
                    ${CMAKE_SOURCE_DIR}/include/gclib/GStr.h
                    ${CMAKE_SOURCE_DIR}/include/gclib/GStr.cpp
                    ${CMAKE_SOURCE_DIR}/include/gclib/GArgs.h
                    ${CMAKE_SOURCE_DIR}/include/gclib/GArgs.cpp
                    ${CMAKE_SOURCE_DIR}/include/gclib/GBase.h
                    ${CMAKE_SOURCE_DIR}/include/gclib/GBase.cpp)

add_dependencies(tiecov htslib)
target_link_libraries(tiecov -lhts -lz -lbz2 -llzma -pthread -lcrypto -lcurl)
