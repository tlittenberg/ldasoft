# add library
add_library(tools STATIC gmm.c gaussian_mixture_model.c gmm.h)

include_directories ("${PROJECT_SOURCE_DIR}/tools/src/")
include_directories(SYSTEM ${GSL_INCLUDE_DIRS})

install(TARGETS tools DESTINATION lib)

add_executable(gaussian_mixture_model gaussian_mixture_model.c gmm.c gmm.h)
target_link_libraries(gaussian_mixture_model m)
target_link_libraries(gaussian_mixture_model ${GSL_LIBRARIES})
install(TARGETS gaussian_mixture_model DESTINATION bin)

