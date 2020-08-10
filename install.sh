pushd tools/src
make clean
make install
popd

pushd gbmcmc/src
make clean
make install
popd

pushd gbfisher/src
make clean
make install
popd

