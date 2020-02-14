#!/bin/bash
source scl_source enable devtoolset-8
source loadLSST.bash
setup -t sims_w_2019_42 lsst_sims
setup -t DC2production throughputs
setup -t DC2production sims_skybrightness_data
pip install nose
pip install coveralls
pip install pylint
eups declare sims_truthcatalog -r ${TRAVIS_BUILD_DIR} -t current
setup sims_truthcatalog
cd ${TRAVIS_BUILD_DIR}
scons
nosetests -s --with-coverage --cover-package=desc.sims_truthcatalog
