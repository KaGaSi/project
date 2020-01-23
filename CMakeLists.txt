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

set(LIBRARIES Options
              AnalysisTools
              Errors)
foreach(library ${LIBRARIES})
  add_library(${library} ${library}.c)
endforeach(library)

set(UTILS AddToSystem
          Aggregates
          Aggregates-NotSameBeads
          AngleMolecules
          Average
          BondLength
          Config
          Config_from_xyz
          DensityAggregates
          DensityBox
          DensityMolecules
          DihedralMolecules
          DistrAgg
          GenSystem
          GenLayers
          GyrationAggregates
          GyrationMolecules
          Info
          JoinAggregates
          JoinRuns
          lmp_data
          lmp_vtf
          Orientation
          PairCorrel
          PotentialAggregates
          SelectedVcf
          TransformVsf
          VisualizeAgg)
        foreach(utility ${UTILS})
  add_executable(${utility} Utility/${utility}.c)
  target_link_libraries(${utility} ${LIBRARIES} m)
endforeach(utility)

#set(ONETIME PlaceSquares)
#        foreach(utility ${ONETIME})
#  add_executable(${utility} OneTime/${utility}.c)
#  target_link_libraries(${utility} m ${LIBRARIES})
#endforeach(utility)

# this somehow specifies fortran compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
set(F_UTILS traject-v2_5
            traject-v2_6)
foreach(utility ${F_UTILS})
  add_executable(${utility} Utility/${utility}.f90)
endforeach(utility)
