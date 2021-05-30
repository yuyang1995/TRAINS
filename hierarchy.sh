#!/bin/bash
# primary parameter file
Param=$1
# secondary parameter file
ParamPT=$(grep '^ParamPT' ${Param} | awk '{print $2}')
FileDir=$(grep '^FileDir' ${Param} | awk '{print $2}')
# FlagMethod=1
sed -i 's/\(FlagMethod\s*\)[0-9]/\10/' ${FileDir}/${ParamPT}
# FlagPrior=1
sed -i 's/\(FlagPrior\s*\)[0-9]/\11/' ${FileDir}/${ParamPT}

# Minimum and Maximum number of GW sources
MinNs=1
MaxNs=$2
# Number of cores used in parallel sampling
Ncores=20
for Ns in $(seq ${MinNs} 1 ${MaxNs}); do
    # SourceNumber=Ns
    echo "The number of GW sources: "${Ns}
    echo ""
    sed -i 's/\(SourceNumber\s*\)[0-9]\+/\1'${Ns}'/' ${FileDir}/${ParamPT}
    # MaxNumberSaves = 5000
    NSaves=5000
    sed -i 's/\(MaxNumberSaves\s*\)[0-9]\+/\1'${NSaves}'/' ${FileDir}/options/OPTIONSPT

    # do diffusive nested sampling
    mpirun -np $Ncores ./trains ${Param} >output_${Ns}.txt
    # get the size of the posterior sample
    Nps=$(wc -l ${FileDir}/data/posterior_sample_pt.txt | awk '{print $1}')
    # Done creating levels?
    Done=$(grep Done output_${Ns}.txt)

    while [ $Nps -le 100 -o -z "$Done" ]; do
        # prepare the restart file
        cp ${FileDir}/data/restart_pt_dnest.txt_${NSaves} ${FileDir}/data/restart_pt_dnest.txt
        # restore the range file
        if [ $Ns -gt 1 ]; then
            cp ${FileDir}/data/data_$(expr ${Ns} - 1)/range_pt.txt ${FileDir}/data
        fi
        # increase MaxNumberSaves
        NSaves=$(expr ${NSaves} + 5000)
        sed -i 's/\(MaxNumberSaves\s*\)[0-9]\+/\1'${NSaves}'/' ${FileDir}/options/OPTIONSPT

        # restart the sampling
        mpirun -np $Ncores ./trains -r ${Param} >>output_${Ns}.txt
        # update the size of the posterior sample
        Nps=$(wc -l ${FileDir}/data/posterior_sample_pt.txt | awk '{print $1}')
        # Done creating levels?
        Done=$(grep Done output_${Ns}.txt)
    done

    # evidence, information and effective sample size
    grep Compress output_${Ns}.txt | head -1
    grep -A2 Z output_${Ns}.txt | tail -3
    echo ""
    # time used
    grep Time output_${Ns}.txt
    echo ""

    # save the data
    mkdir -p ${FileDir}/data/data_${Ns}
    cp ${FileDir}/data/*_pt.txt ${FileDir}/data/data_${Ns}
    cp ${FileDir}/data/DNEST_OPTIONS ${FileDir}/data/data_${Ns}
    cp ${FileDir}/data/restart_pt_dnest.txt_${NSaves} ${FileDir}/data/data_${Ns}/restart_pt_dnest.txt
    rm ${FileDir}/data/restart_pt_dnest.txt*
    echo ""
done
