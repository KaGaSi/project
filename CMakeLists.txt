cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
set(CMAKE_LEGACY_CYGWIN_WIN32 0) # something for cygwin

project(project)
enable_language (C Fortran)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# c compiler flags
set(CMAKE_C_FLAGS "-std=c11 -Wall -Werror") # always stop on errors
#set(CMAKE_C_FLAGS "-std=c11 -Wall") # only show errors
set(CMAKE_C_FLAGS_DEBUG "-ggdb -g") # for -DCMAKE_BUILD_TYPE=Debug
set(CMAKE_C_FLAGS_RELEASE "-O3") # for -DCMAKE_BUILD_TYPE=Release

set(LIBRARIES AnalysisTools
              Options)
foreach(library ${LIBRARIES})
  add_library(${library} ${library}.c)
endforeach(library)

set(COMMON Aggregates
           Aggregates-NotSameBeads
           Average
           BondLength
           Config
           DensityAggregates
           DensityMolecules
           DistrAgg
           GyrationAggregates
           GyrationMolecules
           JoinAggregates
#          JoinRuns
           PairCorrel
#          PairCorrelPerAgg
           SelectedVcf
           TransformVsf)
         foreach(utility ${COMMON})
  add_executable(${utility} Common/${utility}.c)
  target_link_libraries(${utility} m ${LIBRARIES})
endforeach(utility)

set(LINEAR EndToEnd)
         foreach(utility ${LINEAR})
  add_executable(${utility} Linear/${utility}.c)
  target_link_libraries(${utility} m ${LIBRARIES})
endforeach(utility)

# this somehow specifies fortran compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
set(F_UTILITIES traject-v2_5
                traject-v2_6)
foreach(utility ${F_UTILITIES})
  add_executable(${utility} Common/${utility}.f90)
endforeach(utility)
