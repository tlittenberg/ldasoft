# add library
add_library(tools STATIC GMM_with_EM.c gaussian_mixture_model.c)

include_directories ("${PROJECT_SOURCE_DIR}/tools/src/")
include_directories(SYSTEM ${GSL_INCLUDE_DIRS})
target_link_libraries(tools ${GSL_LIBRARIES})

install(TARGETS tools DESTINATION lib)

add_executable(gaussian_mixture_model gaussian_mixture_model.c GMM_with_EM.c)
target_link_libraries(gaussian_mixture_model ${GSL_LIBRARIES})
install(TARGETS gaussian_mixture_model DESTINATION bin)

