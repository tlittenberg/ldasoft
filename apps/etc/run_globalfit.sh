globalfit="${GLASS_PREFIX}/global_fit"
data="${HOME}/ldc/sangria/data/LDC2_sangria_training_v2.h5"
vgb="${GLASS_REPO}/apps/etc/hmcnc.dat"
outdir="${GLASS_TEST}/globalfit"

#path/to/search_sources.dat (contains starting point for MBH sampler)
mbh="${HOME}/ldasoft/apps/etc/"

Tobs=2621440
samples=64
padding=8

#Tobs=7864320
#padding=32

cmd="${globalfit} --h5-data ${data} --h5-no-mbh --sangria --chains 12 --threads 12 --duration ${Tobs} --rundir ${outdir} --known-sources ${vgb} --h5-no-mbh --padding ${padding} --samples ${samples}"


#--mbh-search-path ${mbh}"


echo $cmd
export OMP_NUM_THREADS=12

mpirun -np 4 --oversubscribe $cmd
