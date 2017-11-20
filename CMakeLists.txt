cmake_minimum_required(VERSION 3.0)

project(KFParticleLibrary)

macro(add_target_property _target _prop _value)
   get_target_property(_oldprop "${_target}" ${_prop})
   if(NOT _oldprop)
      set_target_properties("${_target}" PROPERTIES ${_prop} "${_value}")
   else(NOT _oldprop)
      set_target_properties("${_target}" PROPERTIES ${_prop} "${_oldprop} ${_value}")
   endif(NOT _oldprop)
endmacro(add_target_property)

if (CMAKE_VERSION VERSION_LESS "3.1")
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set (CMAKE_CXX_FLAGS "--std=gnu++11 ${CMAKE_CXX_FLAGS}")
endif ()
else ()
    set (CMAKE_CXX_STANDARD 11)
endif ()

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

set(FIXTARGET FALSE CACHE BOOL "Compile for fix target geometry.")

find_package(ROOT REQUIRED COMPONENTS Core Hist Matrix Physics EG Gpad Graf Graf3d RIO MathCore)
find_package(Vc REQUIRED) 

Set(INCLUDE_DIRECTORIES
    ${PROJECT_SOURCE_DIR}/KFParticle
    ${PROJECT_SOURCE_DIR}/KFParticlePerformance
)

Include_Directories( ${CMAKE_SOURCE_DIR} ${ROOT_INCLUDE_DIRS} ${INCLUDE_DIRECTORIES} ${Vc_INCLUDE_DIR})
add_definitions(${ROOT_CXX_FLAGS})

include(${ROOT_USE_FILE})

set(LINK_DIRECTORIES
  ${ROOT_LIBRARY_DIR}
)
 
link_directories( ${LINK_DIRECTORIES})

set (SOURCES
   KFParticle/KFParticleTopoReconstructor.cxx
   KFParticle/KFVertex.cxx	
   KFParticle/KFPTrack.cxx
   KFParticle/KFPTrackVector.cxx
   KFParticle/KFPVertex.cxx
   KFParticle/KFParticlePVReconstructor.cxx
   KFParticle/KFParticleDatabase.cxx
   KFParticle/KFParticleBase.cxx
   KFParticle/KFParticleBaseSIMD.cxx
   KFParticle/KFParticle.cxx
   KFParticle/KFParticleSIMD.cxx
   KFParticle/KFParticleFinder.cxx
   KFParticlePerformance/KFMCVertex.cxx
   KFParticlePerformance/KFParticlePerformanceBase.cxx
   KFParticlePerformance/KFTopoPerformance.cxx
   KFParticlePerformance/KFMCParticle.cxx
)

set (HEADERS
   KFParticle/KFParticleBase.h
   KFParticle/KFParticle.h
   KFParticle/KFVertex.h
   KFParticlePerformance/KFPartEfficiencies.h
   KFParticlePerformance/KFMCParticle.h
)

if(FIXTARGET)
  ROOT_GENERATE_DICTIONARY(G__KF ${HEADERS} LINKDEF KFLinkDef.h OPTIONS "-DDO_TPCCATRACKER_EFF_PERFORMANCE" "-DNonhomogeneousField" "-DCBM" "-DUSE_TIMERS")
  add_library(KFParticle SHARED ${SOURCES} G__KF.cxx)
  target_link_libraries(KFParticle ${ROOT_LIBRARIES} ${Vc_LIBRARIES} )
  add_target_property(KFParticle COMPILE_FLAGS "-DDO_TPCCATRACKER_EFF_PERFORMANCE -DNonhomogeneousField -DCBM -DUSE_TIMERS")
else(FIXTARGET)
  ROOT_GENERATE_DICTIONARY(G__KF ${HEADERS} LINKDEF KFLinkDef.h OPTIONS "-DDO_TPCCATRACKER_EFF_PERFORMANCE" "-DHomogeneousField" "-DUSE_TIMERS")
  add_library(KFParticle SHARED ${SOURCES} G__KF.cxx)
  target_link_libraries(KFParticle ${ROOT_LIBRARIES} ${Vc_LIBRARIES} )
  add_target_property(KFParticle COMPILE_FLAGS "-DDO_TPCCATRACKER_EFF_PERFORMANCE -DHomogeneousField -DUSE_TIMERS")
endif(FIXTARGET)
