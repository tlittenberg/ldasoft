
#GBMCMC library
include_directories ("${PROJECT_SOURCE_DIR}/tools/src/")              
include_directories ("${PROJECT_SOURCE_DIR}/lisa/src/")
include_directories ("${PROJECT_SOURCE_DIR}/util/src/")
include_directories ("${PROJECT_SOURCE_DIR}/gbmcmc/src/")
include_directories ("${PROJECT_SOURCE_DIR}/noise/src/")
include_directories(SYSTEM ${GSL_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS})

add_executable(global_fit GlobalFit.c GalacticBinaryWrapper.h NoiseWrapper.h GalacticBinaryWrapper.c NoiseWrapper.c VerificationBinaryWrapper.c VerificationBinaryWrapper.h MBHWrapper.c MBHWrapper.h)
target_link_libraries(global_fit ${GSL_LIBRARIES} ${HDF5_LIBRARIES})
target_link_libraries(global_fit m)
target_link_libraries(global_fit tools)
target_link_libraries(global_fit lisa)
target_link_libraries(global_fit util)
target_link_libraries(global_fit gbmcmc)
target_link_libraries(global_fit noise)
target_link_libraries(global_fit mpi)
target_link_libraries(global_fit MBH::mbh)
if(APPLE)
   target_link_libraries(global_fit omp)
endif()
install(TARGETS global_fit DESTINATION bin)

