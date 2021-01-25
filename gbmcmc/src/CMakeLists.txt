# add library
add_library(gbmcmc STATIC GalacticBinaryFStatistic.c GalacticBinaryMatch.c GalacticBinaryPrior.c GalacticBinaryWaveform.c gb_mcmc_chirpmass.c
            GalacticBinaryCatalog.c GalacticBinaryIO.c GalacticBinaryMath.c GalacticBinaryProposal.c 
            GalacticBinaryData.c GalacticBinaryMCMC.c GalacticBinaryModel.c GalacticBinaryResidual.c gb_mcmc_brans_dicke.c gb_mcmc.c)

include_directories ("${PROJECT_SOURCE_DIR}/tools/src/")              
include_directories ("${PROJECT_SOURCE_DIR}/lisa/src/")
include_directories ("${PROJECT_SOURCE_DIR}/gbmcmc/src/")
include_directories(SYSTEM ${GSL_INCLUDE_DIRS})
target_link_libraries(gbmcmc tools omp ${GSL_LIBRARIES})

install(TARGETS gbmcmc DESTINATION lib)
install(DIRECTORY "./" DESTINATION include FILES_MATCHING PATTERN "*.h")

add_executable(gb_mcmc gb_mcmc.c GalacticBinaryMCMC.c GalacticBinaryIO.c GalacticBinaryModel.c GalacticBinaryWaveform.c GalacticBinaryMath.c GalacticBinaryData.c GalacticBinaryPrior.c GalacticBinaryProposal.c GalacticBinaryFStatistic.c)
target_link_libraries(gb_mcmc ${GSL_LIBRARIES})
target_link_libraries(gb_mcmc tools)
target_link_libraries(gb_mcmc lisa)
target_link_libraries(gb_mcmc omp)
target_link_libraries(gb_mcmc hdf5)
install(TARGETS gb_mcmc DESTINATION bin)

add_executable(gb_catalog GalacticBinaryCatalog.c GalacticBinaryIO.c GalacticBinaryModel.c GalacticBinaryWaveform.c GalacticBinaryMath.c GalacticBinaryData.c GalacticBinaryPrior.c GalacticBinaryProposal.c GalacticBinaryFStatistic.c)
target_link_libraries(gb_catalog ${GSL_LIBRARIES})
target_link_libraries(gb_catalog tools)
target_link_libraries(gb_catalog lisa)
target_link_libraries(gb_catalog omp)
target_link_libraries(gb_catalog hdf5)
install(TARGETS gb_catalog DESTINATION bin)

add_executable(gb_residual GalacticBinaryResidual.c GalacticBinaryIO.c GalacticBinaryModel.c GalacticBinaryWaveform.c GalacticBinaryMath.c GalacticBinaryData.c GalacticBinaryPrior.c GalacticBinaryProposal.c GalacticBinaryFStatistic.c)
target_link_libraries(gb_residual ${GSL_LIBRARIES})
target_link_libraries(gb_residual tools)
target_link_libraries(gb_residual lisa)
target_link_libraries(gb_residual omp)
target_link_libraries(gb_residual hdf5)
install(TARGETS gb_residual DESTINATION bin)

add_executable(gb_mcmc_chirpmass gb_mcmc_chirpmass.c)
target_link_libraries(gb_mcmc_chirpmass ${GSL_LIBRARIES})
install(TARGETS gb_mcmc_chirpmass DESTINATION bin)

add_executable(gb_match GalacticBinaryMatch.c GalacticBinaryIO.c GalacticBinaryModel.c GalacticBinaryWaveform.c GalacticBinaryMath.c GalacticBinaryData.c GalacticBinaryPrior.c GalacticBinaryProposal.c GalacticBinaryFStatistic.c)
target_link_libraries(gb_match ${GSL_LIBRARIES})
target_link_libraries(gb_match tools)
target_link_libraries(gb_match lisa)
target_link_libraries(gb_match omp)
target_link_libraries(gb_match hdf5)
install(TARGETS gb_match DESTINATION bin)

add_executable(gb_mcmc_brans_dicke gb_mcmc_brans_dicke.c)
target_link_libraries(gb_mcmc_brans_dicke ${GSL_LIBRARIES})
install(TARGETS gb_mcmc_brans_dicke DESTINATION bin)
