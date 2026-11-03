#==============================================================================
#
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Domain
  PRIVATE
    LadrunoContactDomain.cpp
    LadrunoContactSurface.cpp
    LadrunoContactCommands.cpp
  PUBLIC
    LadrunoContactDomain.h
    LadrunoContactSurface.h
    LadrunoContactProjection.h
    LadrunoContactNormalField.h
    LadrunoContactBucketSort.h
    LadrunoMortarKernel.h
    LadrunoFrictionKernel.h
    LadrunoEdgeKernel.h
    LadrunoContact2DKernel.h
)

target_include_directories(OPS_Domain PUBLIC ${CMAKE_CURRENT_LIST_DIR})
