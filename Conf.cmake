#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#     (c) Copyright 1999-2021 The Regents of the University of California
#                             All Rights Reserved
# (Copyright and Disclaimer @ http://www.berkeley.edu/OpenSees/copyright.html)
#
#==============================================================================


#==============================================================================
#                           Select Default Executable
#==============================================================================
set(OPS_FINAL_TARGET "OpenSees" CACHE STRING "OpenSees final target")
set_property(CACHE OPS_FINAL_TARGET PROPERTY STRINGS
    "G3;OpenSees;OpenSeesMP;OpenSeesSP;OpenSeesPy;OpenSeesPyMP")


#==============================================================================
#                            Basic Switches
#==============================================================================

option(OPS_Use_Dev_Directories
  "Include files in DEVELOPER directory"                   OFF)

option(FMK
  "Special FMK Code"                                       OFF)

set(OPS_Use_Graphics_Option
  None
  # Base
  # OpenGL
)

