cmake_minimum_required(VERSION 2.8.9) 
#-----------------------------------------------------------------------------
set(MODULE_NAME ProstateSegmenter)
#-----------------------------------------------------------------------------
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})
#-----------------------------------------------------------------------------
set(MODULE_INCLUDE_DIRECTORIES)
set(MODULE_SRCS)
set(MODULE_TARGET_LIBRARIES ${ITK_LIBRARIES})
#-----------------------------------------------------------------------------
SEMMacroBuildCLI(
  NAME ${MODULE_NAME}
  TARGET_LIBRARIES ${MODULE_TARGET_LIBRARIES}
  INCLUDE_DIRECTORIES ${MODULE_INCLUDE_DIRECTORIES}
  ADDITIONAL_SRCS ${MODULE_SRCS}
  )
#-----------------------------------------------------------------------------
if(BUILD_TESTING)
  add_subdirectory(Testing)
endif()

