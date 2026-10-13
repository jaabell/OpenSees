#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

# Only the general pair is built. PARDISOSymLinSOE/PARDISOSymLinSolver are
# not compiled; PARDISOGenLinSOE covers the symmetric cases (-matrixType).
target_sources(OPS_SysOfEqn
    PRIVATE
        PARDISOGenLinSOE.cpp
        PARDISOGenLinSolver.cpp

    PUBLIC
        PARDISOGenLinSOE.h
        PARDISOGenLinSolver.h
)

target_include_directories(OPS_SysOfEqn PUBLIC ${CMAKE_CURRENT_LIST_DIR})

# mkl_pardiso.h, mkl_types.h, mkl_service.h
if(DEFINED MKL_ROOT)
    target_include_directories(OPS_SysOfEqn PUBLIC ${MKL_ROOT}/include)
endif()
