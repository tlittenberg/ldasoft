# add library
add_library(gbmcmc STATIC BayesLine.c GalacticBinaryFStatistic.c GalacticBinaryMatch.c GalacticBinaryPrior.c GalacticBinaryWaveform.c gb_mcmc_chirpmass.c
            GalacticBinaryCatalog.c GalacticBinaryIO.c GalacticBinaryMath.c GalacticBinaryProposal.c LISA.c
            GalacticBinaryData.c GalacticBinaryMCMC.c GalacticBinaryModel.c GalacticBinaryResidual.c gb_mcmc_brans_dicke.c)
              
include_directories ("${PROJECT_SOURCE_DIR}/gbmcmc/src/")
include_directories(SYSTEM ${GSL_INCLUDE_DIRS})
target_link_libraries(gbmcmc ${GSL_LIBRARIES})

install(TARGETS gbmcmc DESTINATION lib)


add_executable(gb_mcmc GalacticBinaryMCMC.c LISA.c GalacticBinaryIO.c GalacticBinaryModel.c GalacticBinaryWaveform.c GalacticBinaryMath.c GalacticBinaryData.c GalacticBinaryPrior.c GalacticBinaryProposal.c GalacticBinaryFStatistic.c)
target_link_libraries(gb_mcmc ${GSL_LIBRARIES})
install(TARGETS gb_mcmc DESTINATION bin)

add_executable(gb_catalog GalacticBinaryCatalog.c LISA.c GalacticBinaryIO.c GalacticBinaryModel.c GalacticBinaryWaveform.c GalacticBinaryMath.c GalacticBinaryData.c GalacticBinaryPrior.c GalacticBinaryProposal.c GalacticBinaryFStatistic.c)
target_link_libraries(gb_catalog ${GSL_LIBRARIES})
install(TARGETS gb_catalog DESTINATION bin)

add_executable(gb_residual GalacticBinaryResidual.c LISA.c GalacticBinaryIO.c GalacticBinaryModel.c GalacticBinaryWaveform.c GalacticBinaryMath.c GalacticBinaryData.c GalacticBinaryPrior.c GalacticBinaryProposal.c GalacticBinaryFStatistic.c)
target_link_libraries(gb_residual ${GSL_LIBRARIES})
install(TARGETS gb_residual DESTINATION bin)

add_executable(gb_mcmc_chirpmass gb_mcmc_chirpmass.c)
target_link_libraries(gb_mcmc_chirpmass ${GSL_LIBRARIES})
install(TARGETS gb_mcmc_chirpmass DESTINATION bin)

add_executable(gb_match GalacticBinaryMatch.c LISA.c GalacticBinaryIO.c GalacticBinaryModel.c GalacticBinaryWaveform.c GalacticBinaryMath.c GalacticBinaryData.c GalacticBinaryPrior.c GalacticBinaryProposal.c GalacticBinaryFStatistic.c)
target_link_libraries(gb_match ${GSL_LIBRARIES})
install(TARGETS gb_match DESTINATION bin)

add_executable(gb_mcmc_brans_dicke gb_mcmc_brans_dicke.c)
target_link_libraries(gb_mcmc_brans_dicke ${GSL_LIBRARIES})
install(TARGETS gb_mcmc_brans_dicke DESTINATION bin)
