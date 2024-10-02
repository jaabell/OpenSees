#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Element
    PRIVATE
    FourNodeTetrahedron.cpp
    TenNodeTetrahedron.cpp
    TenNodeTetrahedronThermal.cpp
    SixNodeBoundryCondition.cpp
    PUBLIC
    FourNodeTetrahedron.h
    TenNodeTetrahedron.h
    TenNodeTetrahedronThermal.h
    SixNodeBoundryCondition.h
)

target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})
