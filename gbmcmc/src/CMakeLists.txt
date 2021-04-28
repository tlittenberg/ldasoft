# add library
add_library(gbmcmc STATIC GalacticBinaryFStatistic.c GalacticBinaryMatch.c GalacticBinaryPrior.c GalacticBinaryWaveform.c
            GalacticBinaryCatalog.c GalacticBinaryIO.c GalacticBinaryMath.c GalacticBinaryProposal.c 
            GalacticBinaryData.c GalacticBinaryMCMC.c GalacticBinaryModel.c GalacticBinaryResidual.c
            GalacticBinary.h
            GalacticBinaryMCMC.h GalacticBinaryIO.h GalacticBinaryModel.h GalacticBinaryWaveform.h 
            GalacticBinaryMath.h GalacticBinaryData.h GalacticBinaryPrior.h GalacticBinaryProposal.h
            GalacticBinaryFStatistic.h GalacticBinaryCatalog.h)

include_directories ("${PROJECT_SOURCE_DIR}/tools/src/")              
include_directories ("${PROJECT_SOURCE_DIR}/lisa/src/")
include_directories ("${PROJECT_SOURCE_DIR}/gbmcmc/src/")
include_directories(SYSTEM ${GSL_INCLUDE_DIRS})
target_link_libraries(gbmcmc tools omp ${GSL_LIBRARIES})

install(TARGETS gbmcmc DESTINATION lib)
install(DIRECTORY "./" DESTINATION include FILES_MATCHING PATTERN "*.h")

add_executable(gb_mcmc gb_mcmc.c GalacticBinary.h 
                GalacticBinaryMCMC.h GalacticBinaryIO.h GalacticBinaryModel.h GalacticBinaryWaveform.h GalacticBinaryMath.h GalacticBinaryData.h GalacticBinaryPrior.h GalacticBinaryProposal.h GalacticBinaryFStatistic.h GalacticBinaryCatalog.h
                GalacticBinaryMCMC.c GalacticBinaryIO.c GalacticBinaryModel.c GalacticBinaryWaveform.c GalacticBinaryMath.c GalacticBinaryData.c GalacticBinaryPrior.c GalacticBinaryProposal.c GalacticBinaryFStatistic.c GalacticBinaryCatalog.c)
                
target_link_libraries(gb_mcmc ${GSL_LIBRARIES})
target_link_libraries(gb_mcmc tools)
target_link_libraries(gb_mcmc lisa)
target_link_libraries(gb_mcmc omp)
target_link_libraries(gb_mcmc hdf5)
install(TARGETS gb_mcmc DESTINATION bin)

add_executable(gb_catalog gb_catalog.c
                GalacticBinaryCatalog.h GalacticBinaryIO.h GalacticBinaryModel.h GalacticBinaryWaveform.h GalacticBinaryMath.h GalacticBinaryData.h GalacticBinaryPrior.h GalacticBinaryProposal.h GalacticBinaryFStatistic.h 
                GalacticBinaryCatalog.c GalacticBinaryIO.c GalacticBinaryModel.c GalacticBinaryWaveform.c GalacticBinaryMath.c GalacticBinaryData.c GalacticBinaryPrior.c GalacticBinaryProposal.c GalacticBinaryFStatistic.c)
target_link_libraries(gb_catalog ${GSL_LIBRARIES})
target_link_libraries(gb_catalog tools)
target_link_libraries(gb_catalog lisa)
target_link_libraries(gb_catalog omp)
target_link_libraries(gb_catalog hdf5)
install(TARGETS gb_catalog DESTINATION bin)

add_executable(gb_residual GalacticBinaryResidual.c GalacticBinaryCatalog.c GalacticBinaryIO.c GalacticBinaryModel.c GalacticBinaryWaveform.c GalacticBinaryMath.c GalacticBinaryData.c GalacticBinaryPrior.c GalacticBinaryProposal.c GalacticBinaryFStatistic.c)
target_link_libraries(gb_residual ${GSL_LIBRARIES})
target_link_libraries(gb_residual tools)
target_link_libraries(gb_residual lisa)
target_link_libraries(gb_residual omp)
target_link_libraries(gb_residual hdf5)
install(TARGETS gb_residual DESTINATION bin)

add_executable(gb_match GalacticBinaryMatch.c GalacticBinaryCatalog.c GalacticBinaryIO.c GalacticBinaryModel.c GalacticBinaryWaveform.c GalacticBinaryMath.c GalacticBinaryData.c GalacticBinaryPrior.c GalacticBinaryProposal.c GalacticBinaryFStatistic.c)
target_link_libraries(gb_match ${GSL_LIBRARIES})
target_link_libraries(gb_match tools)
target_link_libraries(gb_match lisa)
target_link_libraries(gb_match omp)
target_link_libraries(gb_match hdf5)
install(TARGETS gb_match DESTINATION bin)

