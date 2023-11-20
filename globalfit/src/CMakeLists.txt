
#GBMCMC library
include_directories ("${PROJECT_SOURCE_DIR}/tools/src/")              
include_directories ("${PROJECT_SOURCE_DIR}/data/src/")
include_directories ("${PROJECT_SOURCE_DIR}/ucb/src/")
include_directories ("${PROJECT_SOURCE_DIR}/noise/src/")
include_directories(SYSTEM ${GSL_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS})

add_executable(global_fit globalfit.c ucb_wrapper.h noise_wrapper.h ucb_wrapper.c noise_wrapper.c vgb_wrapper.c vgb_wrapper.h mbh_wrapper.c mbh_wrapper.h)
target_link_libraries(global_fit ${GSL_LIBRARIES} ${HDF5_LIBRARIES})
target_link_libraries(global_fit m)
target_link_libraries(global_fit tools)
target_link_libraries(global_fit data)
target_link_libraries(global_fit ucb)
target_link_libraries(global_fit noise)
target_link_libraries(global_fit MBH::mbh)
if(APPLE)
   target_link_libraries(global_fit omp)
endif()
install(TARGETS global_fit DESTINATION bin)

