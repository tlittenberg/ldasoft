
#GBMCMC library
include_directories ("${PROJECT_SOURCE_DIR}/utils/src/")
include_directories ("${PROJECT_SOURCE_DIR}/ucb/src/")
include_directories ("${PROJECT_SOURCE_DIR}/noise/src/")
include_directories(SYSTEM ${GSL_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS})

add_executable(global_fit globalfit.c glass_ucb_wrapper.h glass_noise_wrapper.h glass_ucb_wrapper.c glass_noise_wrapper.c glass_vgb_wrapper.c glass_vgb_wrapper.h glass_mbh_wrapper.c glass_mbh_wrapper.h)
target_link_libraries(global_fit ${GSL_LIBRARIES} ${HDF5_LIBRARIES})
target_link_libraries(global_fit m)
target_link_libraries(global_fit glass_utils)
target_link_libraries(global_fit glass_ucb)
target_link_libraries(global_fit glass_noise)
target_link_libraries(global_fit MBH::mbh)
if(APPLE)
   target_link_libraries(global_fit omp)
endif()
install(TARGETS global_fit DESTINATION bin)

