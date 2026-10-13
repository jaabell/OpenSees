#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_SysOfEqn
  PRIVATE
    DomainSolver.cpp
    LinearSOE.cpp
    LinearSOESolver.cpp
  PUBLIC
    DomainSolver.h
    LinearSOE.h
    LinearSOESolver.h
)

target_include_directories(OPS_SysOfEqn PUBLIC ${CMAKE_CURRENT_LIST_DIR})

add_subdirectory(bandGEN)
add_subdirectory(bandSPD)
add_subdirectory(diagonal)
add_subdirectory(fullGEN)
add_subdirectory(sparseGEN)
add_subdirectory(sparseSYM)
add_subdirectory(umfGEN)
add_subdirectory(sparsePython)

add_subdirectory(profileSPD)
#add_subdirectory(cg)
#add_subdirectory(petsc)
#add_subdirectory(mumps)
# MKL PARDISO (system Pardiso): OPS_HAVE_PARDISO is set in the top-level
# CMakeLists.txt when MKL is found on Windows or MKL_PARDISO_LINUX is ON.
if(OPS_HAVE_PARDISO)
  add_subdirectory(pardiso)
endif()
add_subdirectory(itpack)

