#!/bin/bash
# ##################################################################
# run this script in its own directory to cp all the python codes 
# and data files into $HOME/opt/script                     
# ##################################################################

#Note: the files in the parameters directory still have to be copied to 
#the specific simulation directory to make it run.

export s=$HOME/opt/script

mkdir -p $HOME/opt/script

dir=`pwd`

#cp files in create_project_tools/
cp -v * $s/

cd ..

#cp files in results_analysis_tools/
if [ ! -d results_analysis_tools ]; then
	echo -e "\n\nCan't find results_analysis_tools directory! Go figure...\n\n"
	exit
fi

cp -v results_analysis_tools/* $s/

echo -e "\nFinished copying all the files and data!"
echo "ls your $s directory to double-check!"
cd $dir
