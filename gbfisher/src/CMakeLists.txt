# add library
add_executable(Confusion_Fit Confusion_Fit.c Subroutines.c arrays.c)

add_executable(Full_Residual Full_Residual.c Subroutines.c arrays.c)

add_executable(Bright_Remove Bright_Remove.c Subroutines.c arrays.c)

add_executable(Galaxy Galaxy.c Subroutines.c arrays.c)

add_executable(OrbitFile OrbitFile.c Subroutines.c arrays.c)

add_executable(Setup Setup.c)

add_executable(Fisher_Galaxy Fisher_Galaxy.c Subroutines.c arrays.c)


include_directories ("${PROJECT_SOURCE_DIR}/gbmcmc/src/")
include_directories ("${PROJECT_SOURCE_DIR}/gbfisher/src/")
include_directories(SYSTEM ${GSL_INCLUDE_DIRS})


target_link_libraries(Confusion_Fit GSL::gsl GSL::gslcblas)
target_link_libraries(Confusion_Fit gbmcmc)
target_link_libraries(Confusion_Fit m)

target_link_libraries(Full_Residual GSL::gsl GSL::gslcblas)
target_link_libraries(Full_Residual gbmcmc)
target_link_libraries(Full_Residual m)


target_link_libraries(Bright_Remove GSL::gsl GSL::gslcblas)
target_link_libraries(Bright_Remove gbmcmc)
target_link_libraries(Bright_Remove m)


target_link_libraries(Galaxy GSL::gsl GSL::gslcblas)
target_link_libraries(Galaxy gbmcmc)
target_link_libraries(Galaxy m)


target_link_libraries(OrbitFile GSL::gsl GSL::gslcblas)
target_link_libraries(OrbitFile gbmcmc)
target_link_libraries(OrbitFile m)


target_link_libraries(Fisher_Galaxy GSL::gsl GSL::gslcblas)
target_link_libraries(Fisher_Galaxy gbmcmc)
target_link_libraries(Fisher_Galaxy m)


target_link_libraries(Setup GSL::gsl GSL::gslcblas)
target_link_libraries(Setup gbmcmc)
target_link_libraries(Setup m)
