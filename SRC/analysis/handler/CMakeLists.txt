#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Analysis
    PRIVATE
      ConstraintHandler.cpp 
      PlainHandler.cpp
      PenaltyConstraintHandler.cpp 
      LagrangeConstraintHandler.cpp
      TransformationConstraintHandler.cpp 
      AutoConstraintHandler.cpp
      LadrunoProjectionHandler.cpp
      LadrunoConstraintProjector.cpp
      LadrunoContactHandler.cpp
      LadrunoContactFE.cpp
    PUBLIC
      ConstraintHandler.h 
      PlainHandler.h
      PenaltyConstraintHandler.h 
      LagrangeConstraintHandler.h
      TransformationConstraintHandler.h 
      AutoConstraintHandler.h
      LadrunoProjectionConsumer.h
      LadrunoProjectionHandler.h
      LadrunoConstraintProjector.h
      LadrunoContactHandler.h
      LadrunoContactFE.h
)

target_include_directories(OPS_Analysis PUBLIC ${CMAKE_CURRENT_LIST_DIR})
