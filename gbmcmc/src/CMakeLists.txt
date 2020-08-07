# add library
add_library(gbmcmc STATIC BayesLine.c GalacticBinaryFStatistic.c GalacticBinaryMatch.c GalacticBinaryPrior.c GalacticBinaryWaveform.c gb_mcmc_chirpmass.c
            GalacticBinaryCatalog.c GalacticBinaryIO.c GalacticBinaryMath.c GalacticBinaryProposal.c LISA.c
            GalacticBinaryData.c GalacticBinaryMCMC.c GalacticBinaryModel.c GalacticBinaryResidual.c gb_mcmc_brans_dicke.c)
              
include_directories ("${PROJECT_SOURCE_DIR}/gbmcmc/src/")
include_directories(SYSTEM ${GSL_INCLUDE_DIRS})
target_link_libraries(gbmcmc ${GSL_LIBRARIES})


install(TARGETS gbmcmc DESTINATION lib)
