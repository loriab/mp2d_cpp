cmake_minimum_required(VERSION 3.5 FATAL_ERROR)  # 3.5 avoid invalid C compiler

project(mp2d
        VERSION 0.1
        LANGUAGES CXX)
set(mp2d_AUTHORS      "Chandler S. Greenwell and Gregory J. O. Beran")
set(mp2d_DESCRIPTION  "Dispersion correction to MP2")
set(mp2d_URL          "https://github.com/loriab/mp2d_cpp")
set(mp2d_LICENSE      "MIT")

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

################  Options: Overview and Not Otherwise Mentioned  ###############

#  <<<  CMake build overview  >>>
#
#    >>> ls
#    Dispersion.h LICENSE Params.txt ...
#    >>> cmake -H. -Bobjdir -DCMAKE_INSTALL_PREFIX=/path/to/install-mp2d ...
#    ...
#    -- Generating done
#    -- Build files have been written to: /current/dir/objdir
#    >>> cd objdir && make -j`getconf _NPROCESSORS_ONLN`
#    >>> make install

#  <<<  Required build dependencies that MP2D can't build itself  >>>
#
#    - CMake (e.g., `conda install cmake`)
#    - C++ compiler

############################  Options: Build How?  #############################
include(GNUInstallDirs)
include(psi4OptionsTools)
include(autocmake_safeguards)

option_with_default(CMAKE_BUILD_TYPE "Build type (Release or Debug)" Release)

#################################  Main Project  ################################

# <<<  Build >>>

configure_file(Coord_Num.cpp
               ${CMAKE_CURRENT_BINARY_DIR}/Coord_Num.cpp @ONLY)

set(sources_list Main.cpp
                 Dispersion.cpp
                 ${CMAKE_CURRENT_BINARY_DIR}/Coord_Num.cpp)

add_executable(mp2d ${sources_list})
target_include_directories(mp2d PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# <<<  Install  >>>

install(TARGETS mp2d
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

install(FILES GrimmeC6.txt Params.txt UCHF_daug_C6.txt
        DESTINATION ${CMAKE_INSTALL_DATADIR}/mp2d)
