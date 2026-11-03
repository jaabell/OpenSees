#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Recorder
    PRIVATE
      DamageRecorder.cpp
      DatastoreRecorder.cpp
      DriftRecorder.cpp
      ElementRecorder.cpp
      ElementRecorderRMS.cpp
      EnvelopeDriftRecorder.cpp
      EnvelopeElementRecorder.cpp
      EnvelopeNodeRecorder.cpp
      GSA_Recorder.cpp
      GmshRecorder.cpp
      MaxNodeDispRecorder.cpp
      NodeRecorder.cpp
      NodeRecorderRMS.cpp
      NormElementRecorder.cpp
      NormEnvelopeElementRecorder.cpp
      PatternRecorder.cpp
      PVDRecorder.cpp      
      Recorder.cpp
      RemoveRecorder.cpp
      VTK_Recorder.cpp
      EnergyBalanceRecorder.cpp
    PUBLIC
      DamageRecorder.h
      DatastoreRecorder.h
      DriftRecorder.h
      ElementRecorder.h
      ElementRecorderRMS.h
      EnvelopeDriftRecorder.h
      EnvelopeElementRecorder.h
      EnvelopeNodeRecorder.h
      GSA_Recorder.h
      GmshRecorder.h
      MaxNodeDispRecorder.h
      NodeRecorder.h
      NodeRecorderRMS.h
      NormElementRecorder.h
      NormEnvelopeElementRecorder.h
      PatternRecorder.h
      PVDRecorder.h      
      Recorder.h
      RemoveRecorder.h
      VTK_Recorder.h
      EnergyBalanceRecorder.h
      EnergyBalanceKernel.h
)


if(HDF5_FOUND)
   if (HDF5_VERSION VERSION_GREATER_EQUAL 1.12.0)
      target_sources(OPS_Recorder
        PRIVATE
          MPCORecorder.cpp
          VTKHDF_Recorder.cpp
        PUBLIC
         MPCORecorder.h
         VTKHDF_Recorder.h
      )
      # recorder ladruno (HDF5 .ladruno files)
      target_sources(OPS_Recorder
        PRIVATE
          LadrunoRecorder.cpp
          Ladruno_NodeResults.cpp
          Ladruno_ElementResults.cpp
          Ladruno_DomainResults.cpp
          Ladruno_Sinks.cpp
          LadrunoMonitorRecorder.cpp
          LadrunoMonitorSink.cpp
        PUBLIC
          LadrunoRecorder.h
          Ladruno_Types.h
          Ladruno_Hdf5.h
          Ladruno_ResultIO.h
          Ladruno_LaunchEnv.h
          Ladruno_NodeResults.h
          Ladruno_ElementResults.h
          Ladruno_DomainResults.h
          Ladruno_Sinks.h
          LadrunoMonitorRecorder.h
          LadrunoMonitorSink.h
      )
   endif()
endif()

target_sources(OPS_Paraview
    PRIVATE
      PVDRecorder.cpp
    PUBLIC
      PVDRecorder.h
)

target_sources(OPS_Graphics
    PRIVATE
      AlgorithmIncrements.cpp
      FilePlotter.cpp
    PUBLIC
      AlgorithmIncrements.h
      FilePlotter.h
)

add_subdirectory(response)

