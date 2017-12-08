#$ -S /bin/bash
#$ -cwd
#$ -j y

cd ~/diffsplice_0.1.2beta

ResultsDir=~/ResultsSim${1}

mkdir -p $ResultsDir

./diffsplice sim_settings.cfg datafile_sim${1}.cfg ${ResultsDir} > ${ResultsDir}/diffsplicev012run_sim${1}.log
