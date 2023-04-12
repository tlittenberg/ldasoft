# add library
add_library(gbmcmc STATIC GalacticBinaryFStatistic.c GalacticBinaryPrior.c GalacticBinaryWaveform.c
            GalacticBinaryCatalog.c GalacticBinaryIO.c GalacticBinaryMath.c GalacticBinaryProposal.c 
            GalacticBinaryData.c GalacticBinaryMCMC.c GalacticBinaryModel.c GalacticBinaryResidual.c
            GalacticBinary.h
            GalacticBinaryMCMC.h GalacticBinaryIO.h GalacticBinaryModel.h GalacticBinaryWaveform.h 
            GalacticBinaryMath.h GalacticBinaryData.h GalacticBinaryPrior.h GalacticBinaryProposal.h
            GalacticBinaryFStatistic.h GalacticBinaryCatalog.h)

include_directories ("${PROJECT_SOURCE_DIR}/tools/src/")              
include_directories ("${PROJECT_SOURCE_DIR}/lisa/src/")
include_directories ("${PROJECT_SOURCE_DIR}/util/src/")
include_directories ("${PROJECT_SOURCE_DIR}/gbmcmc/src/")
include_directories(SYSTEM ${GSL_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS})
target_link_libraries(gbmcmc tools m ${GSL_LIBRARIES} ${HDF5_LIBRARIES} util)

install(TARGETS gbmcmc DESTINATION lib)
install(DIRECTORY "./" DESTINATION include FILES_MATCHING PATTERN "*.h")

add_executable(gb_mcmc gb_mcmc.c)
                
target_link_libraries(gb_mcmc ${GSL_LIBRARIES})
target_link_libraries(gb_mcmc ${HDF5_LIBRARIES})
target_link_libraries(gb_mcmc m)
target_link_libraries(gb_mcmc tools)
target_link_libraries(gb_mcmc lisa)
target_link_libraries(gb_mcmc gbmcmc)
target_link_libraries(gb_mcmc util)
if(APPLE) 
   target_link_libraries(gb_mcmc omp)
endif()
install(TARGETS gb_mcmc DESTINATION bin)

add_executable(gb_catalog gb_catalog.c)
target_link_libraries(gb_catalog ${GSL_LIBRARIES})
target_link_libraries(gb_catalog ${HDF5_LIBRARIES})
target_link_libraries(gb_catalog m)
target_link_libraries(gb_catalog tools)
target_link_libraries(gb_catalog lisa)
target_link_libraries(gb_catalog gbmcmc)
target_link_libraries(gb_mcmc util)
if(APPLE)
   target_link_libraries(gb_catalog omp)
endif()
install(TARGETS gb_catalog DESTINATION bin)

add_executable(gb_match GalacticBinaryMatch.c)
target_link_libraries(gb_match ${GSL_LIBRARIES})
target_link_libraries(gb_match ${HDF5_LIBRARIES})
target_link_libraries(gb_match m)
target_link_libraries(gb_match tools)
target_link_libraries(gb_match lisa)
target_link_libraries(gb_match gbmcmc)
target_link_libraries(gb_mcmc util)
if(APPLE)
   target_link_libraries(gb_match omp)
endif()
install(TARGETS gb_match DESTINATION bin)

add_executable(vb_mcmc vb_mcmc.c)
                
target_link_libraries(vb_mcmc ${GSL_LIBRARIES})
target_link_libraries(vb_mcmc ${HDF5_LIBRARIES})
target_link_libraries(vb_mcmc m)
target_link_libraries(vb_mcmc tools)
target_link_libraries(vb_mcmc lisa)
target_link_libraries(gb_mcmc util)
if(APPLE)
   target_link_libraries(vb_mcmc omp)
endif()
target_link_libraries(vb_mcmc gbmcmc)
install(TARGETS vb_mcmc DESTINATION bin)
