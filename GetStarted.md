

# Building LAMMPS with AWSEM #


## What do you need ##
  * **October 2012** version of LAMMPS: http://awsemmd.googlecode.com/files/lammps-9Oct12.tar.gz
  * FFTW 2.1.5: http://www.fftw.org/fftw-2.1.5.tar.gz
  * AWSEM code (see Source->Checkout)


## Step by step guide ##
  * Download and unzip FFTW 2.1.5
  * Install FFTW 2.1.5 somewhere locally using `./conﬁgure --preﬁx=/your_home_directory/local; make; make install`
  * Download and unzip LAMMPS source code
  * Copy all _.cpp_ and _.h_ file from awsemmd directory to LAMMPS\_DIR/src/
  * Go to LAMMPSDIR/src/STUBS/ and execute make command to compile a dummy MPI library
  * Copy awsemmd/Makefiles/Makefile.serial file to LAMMPS\_DIR/src/MAKE/ and adjust the paths according to FFTW installation above
  * Go back to LAMMPS\_DIR/src/ directory and execute `make clean-all; make serial`
  * If everything goes well, you must find lmp\_serial executable in src directory.
See the list of known issues [here](Possible_Issues#.md).

For more information on LAMMPS compilation and inclusion of optional packages see the following link
http://lammps.sandia.gov/doc/Section_start.html#start_2


# Generating project #


The easiest way to generate a LAMMPS-AWSEM project is to do that using a pdb file which contains the structure of interest. You will need to have python and biopython installed on your computer.
  * First of all go to awsemmd/create\_project\_tool/ and execute _install\_tools.sh_ script. This will copy all necessary scripts in /your\_home\_directory/opt/script/.
  * Add /your\_home\_directory/opt/script/ to PATH variable (call _echo $PATH_ to make sure that you did it right)
  * Copy lmp\_serial and the pdb file to a same directory
  * From that directory call _PdbCoords2Lammps.sh pdb\_id project\_name_

After this you should find the following files in your simulation (current) directory
  1. _data._ file which contains atom and bond description
  1. _.seq_ file with protein sequance
  1. _.coord_ with all atom backbone description and which is not needed for LAMMPS-AWSEM simulations
  1. _.in_ which contains an input script

You can modify the _.in_ according your needs. To learn more about the commands you will find there or perhaps some other ones you may use refer to LAMMPS documentation.
http://lammps.sandia.gov/doc/Section_commands.html


# Running simulations #

To run a simulation you will also need several parameter files. Most of those files you can find in awsemmd/parameters/ directory. Bellow is a short description of the files located there
  * **fix\_backbone\_coeff.data** - main parameter file for AWSEM potential. Each term has its own section within _fix\_backbone\_coeff.data_ entitled with a term name enclosed in `[]` brackets. To switch any of the terms off put a "-" after `[]` brackets (or alternat its title in any other way). `[Epsilon]` and `[ABC]` sections must be always on. `[Epsilon]` sets general energy scale of AWSEM potential. `[ABC]` coefficients are used to calculate the positions of implicit backbone atoms from CA and O coordinates.
  * **gamma.dat** - contains gamma parameters for direct or protein/water mediated potential (`[Water]`)
  * **burial\_gamma.dat** - contains gammas for burial potential
  * **para\_one, para\_HB, anti\_one, anti\_HB, anti\_NHB** - coefficients for beta hydrogen bonding potential (`[Dssp_Hdrgn]`)
  * **amh-go.gamma** - example of gamma file used by AMH-Go term
  * **uniform.gamma, alpha.gamma** - sample gamma files for Fragment Memory term

You may also need _ssweight_ file if `[SSWeight]` section in _fix\_backbone\_coeff.data_ in on. This file is typically used to apply a secondary structure bias based on prediction. This file should contain two columns of float numbers; one line for each residue indicating a probability of the residue forming alpha-helix or beta-sheet. _GenSswight.py_ script from awsemmd/create\_project\_tools/ can be used to convert a prediction from JPRED online server (http://www.compbio.dundee.ac.uk/www-jpred/) to the desired format. The input file for _GenSswight.py_ should contain only one line with a prediction string.

After obtaining all necessary parameter files you can use _lmp\_serial < your\_project.in_ command to run a LAMMPS-AWSEM-MD simulation.