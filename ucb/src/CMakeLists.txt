# add library
add_library(ucb STATIC ucb_fstatistic.c ucb_prior.c ucb_waveform.c
            ucb_catalog.c ucb_io.c ucb_proposal.c 
            ucb_data.c ucb_sampler.c ucb_model.c ucb_residual.c
            ucb_sampler.h ucb_io.h ucb_model.h ucb_waveform.h 
            ucb_data.h ucb_prior.h ucb_proposal.h
            ucb_fstatistic.h ucb_catalog.h ucb.h)

include_directories ("${PROJECT_SOURCE_DIR}/tools/src/")              
include_directories ("${PROJECT_SOURCE_DIR}/data/src/")
include_directories ("${PROJECT_SOURCE_DIR}/ucb/src/")
include_directories(SYSTEM ${GSL_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS})

install(TARGETS ucb DESTINATION lib)
install(DIRECTORY "./" DESTINATION include FILES_MATCHING PATTERN "*.h")

