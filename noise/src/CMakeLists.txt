# add library
add_library(noise STATIC Noise.c noise_mcmc.c)

include_directories ("${PROJECT_SOURCE_DIR}/tools/src/")
include_directories ("${PROJECT_SOURCE_DIR}/lisa/src/")
include_directories ("${PROJECT_SOURCE_DIR}/gbmcmc/src/")
include_directories(SYSTEM ${GSL_INCLUDE_DIRS})

install(TARGETS noise DESTINATION lib)
install(DIRECTORY "./" DESTINATION include FILES_MATCHING PATTERN "*.h")

add_executable(noise_mcmc noise_mcmc.c Noise.c)
target_link_libraries(noise_mcmc ${GSL_LIBRARIES})
target_link_libraries(noise_mcmc tools)
target_link_libraries(noise_mcmc lisa)
target_link_libraries(noise_mcmc gbmcmc)
target_link_libraries(noise_mcmc omp)
target_link_libraries(noise_mcmc hdf5)
install(TARGETS noise_mcmc DESTINATION bin)

