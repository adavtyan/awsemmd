# This is a TCL script that was adapted by Nick Schafer from code found on the internet
# It can be used to load a pdb file with color information in the b-factor field of a
# PDB and update that information frame by frame during the playback in VMD.
# This was originally designed to be used with Bobby Kim's local Q calculation, but
# should be suitable for animating any color information in the b-factor field, so long
# as the min, max and midpoint are adjusted properly (according to the range of the data).
###########################################################
# To use this script, use:                                # 
# vmd -e animatecolor.tcl -args movie                     # 
# where "movie" is the prefix to your .pdb and .psf files.# 
###########################################################

# A function that sets the user color equal to the value in the B-factor field
proc animatecolor { fname } {
    # load the pdb
    mol load pdb ${fname}.pdb 
    # add the psf file so that the secondary structure can be properly visualized
    mol addfile ${fname}.psf type {psf} first 0 last -1 step 1 waitfor 1 0

    # select all the atoms in the file that was just loaded
    set all [atomselect top all]
    # start with the first frame
    set frame 0
    # set the input model name for reading
    set in [open ${fname}.pdb r]
    # make a beta variable for storing the b-factor information later
    set beta {}
    # go through the whole file
    while { [gets $in line] != -1 } {
	# look at what is in the first field of the given line
	# if it is "END", set the user field to the value of beta,
	# reinitialize beta and increment the frame
	# if it is "ATOM" or "HETA", get the b-factor value (in cols 61-66)
	# and append it to the beta variable
	switch -- [string range $line 0 3] {
	    END {
		$all frame $frame
		$all set user $beta
		set beta {}
		incr frame
	    }
	    ATOM -
	    HETA {
		lappend beta [string range $line 61 66]
	    }
	}
    }
    # set the coloring method to "User"
    mol modcolor 0 top "User"
    # go to the first frame and ...
    animate goto 0

    color scale min 0.0
    color scale max 1.0
    color midpoint 0.5
    mol scaleminmax top 0 0.33 0.66
}

# make the default load style newcartoon
mol default style newcartoon
# use a red-white-green coloring scheme
color scale method RWG
# set the min, max and midpoint so that 0-.33 is red, .33-.66 is white and .66 to 1 is green
color scale min 0.33
color scale max 1.0
color midpoint 0.5
# call the function above with the first argument, 
# e.g. "movie" if you have movie.pdb and movie.psf in your current directory
# animatecolor $argv


