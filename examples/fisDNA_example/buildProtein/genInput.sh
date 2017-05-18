#########################################################################
# Author: Charlse.Zhang
# Created Time: Fri 04 Apr 2014 09:34:14 AM CDT
# File Name: run.sh
# Description: 
#########################################################################
#!/bin/bash -l

scriptPath="/home/mt42/opt/script"
proteinName="fis"

case "$1" in 
1)  #python $scriptPath/PDBToCoordinates_includeDNA.py $proteinName ${proteinName}.coord
    python $scriptPath/PDBToCoordinates.py $proteinName ${proteinName}.coord
    ;;
2)  #python $scriptPath/CoordinatesToWorkLammpsDataFile_includeDNA.py ${proteinName}.coord data.${proteinName} -b
    python $scriptPath/CoordinatesToWorkLammpsDataFile.py ${proteinName}.coord data.${proteinName} -b
    ;;
3)  python $scriptPath/frag_mem_tools/Pdb2Gro.py $proteinName ${proteinName}_chainA.gro A
    python $scriptPath/frag_mem_tools/Pdb2Gro.py $proteinName ${proteinName}_chainB.gro B
    ;;
4)
    python $scriptPath/GenSswight.py ./secondary_structure/chainA_stride.txt ./secondary_structure/ssweight_chainA
    python $scriptPath/GenSswight.py ./secondary_structure/chainB_stride.txt ./secondary_structure/ssweight_chainB
    cat ./secondary_structure/ssweight_chainA \
        ./secondary_structure/ssweight_chainB > ssweight
    ;;
esac
