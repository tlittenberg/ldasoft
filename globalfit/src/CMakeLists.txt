
include_directories ("${PROJECT_SOURCE_DIR}/tools/src/")              
include_directories ("${PROJECT_SOURCE_DIR}/lisa/src/")
include_directories ("${PROJECT_SOURCE_DIR}/gbmcmc/src/")
include_directories ("${PROJECT_SOURCE_DIR}/noise/src/")
include_directories(SYSTEM ${GSL_INCLUDE_DIRS})

add_executable(global_fit GlobalFit.c GalacticBinaryWrapper.c NoiseWrapper.c)
target_link_libraries(global_fit ${GSL_LIBRARIES})
target_link_libraries(global_fit tools)
target_link_libraries(global_fit lisa)
target_link_libraries(global_fit gbmcmc)
target_link_libraries(global_fit noise)
target_link_libraries(global_fit omp)
target_link_libraries(global_fit hdf5)
install(TARGETS global_fit DESTINATION bin)

