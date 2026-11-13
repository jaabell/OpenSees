#==============================================================================
#
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
# jaabell/OpenSees port of the Ladruno geometry-method layer: linear only
target_sources(OPS_Element
    PRIVATE
    SolidTransformationLinear.cpp
    PUBLIC
    SolidTransformation.h
    SolidTransformationLinear.h
)

target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})
