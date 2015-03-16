

# Secondary Structure Bias #

## Introduction ##

This document explains how to apply secondary structure (SS) bias in simulations with AWSEM, and how to turn SS bias off if you don't need it.

## Ramachandran plot parameters ##

```
[Rama]
2.0
5
 1.3149  15.398 0.15   1.74 0.65 -2.138
1.32016 49.0521 0.25  1.265 0.45  0.318
 1.0264 49.0954 0.65 -1.041 0.25  -0.78
    2.0   419.0  1.0  0.995  1.0  0.820
    2.0  15.398  1.0   2.25  1.0  -2.16

[Rama_P]
2
2.17 105.52 1.0 1.153 0.15  -2.4
2.15 109.09 1.0  0.95 0.15 0.218
```

## SSWeight ##
```
[SSWeight]
0 0 0 1 1 0
0 0 0 0 0 0
```

## File format for _ssweight_ ##
To apply SS bias you need a _ssweight_ file in a specific format, which must correspond to Ramachandran plot definition and to '[SSWeight](SSWeight.md)' section in the parameter file. In the most common case described above this file must contain two columns with float numbers (typically 0.0 or 1.0), with the number of lines equal to the number of residues in the structure of interest. Values of 1.0 in the first column set helical bias to the Ramachandaran plots of corresponding residues and 1.0s in the second column bias towards beta-sheets. Note, that technically both helical and beta-sheet biases can be on simultaneously.

In general case _ssweight_ can have up to 12 columns. The first six will set a per residue weighting factor to the terms of the general Ramachandran plot, and the next six columns will do the same for Proline plot.

## Obtaining SS prediction and _ssweight_ file ##

First you need to obtain an SS prediction using any appropriate software or an online server. For instance, you can use Jpred prediction server (http://www.compbio.dundee.ac.uk/www-jpred/), which takes a protein sequence as an input.

After you run Jpred server with desired sequence you should get to the Result page, which will offer you the prediction results in various of different formats. Choose _Simple HTML_ output, which should look like this.

```
SISSRVKSKRIQLGLNQAELAQKVGTTQQSIEQLENGKTKRPRFLPELASALGVSVDWLLNGTSDSNVR
-HHHHHHHHHHH----HHHHHHHH---HHHHHHHH------HHHHHHHHHHH------H----------
```

Copy the prediction (second) line to a separate file and use _GenSswight.py_ script located in awsemmd/create\_project\_tools/ directory to convert it to _ssweight_ file in the desired format.

## Turning SS Bias off ##
```
[Rama]
2.0
3
 1.3149  15.398 0.15   1.74 0.65 -2.138
1.32016 49.0521 0.25  1.265 0.45  0.318
 1.0264 49.0954 0.65 -1.041 0.25  -0.78
    2.0   419.0  1.0  0.995  1.0  0.820
    2.0  15.398  1.0   2.25  1.0  -2.16
```
When compared with the default fix\_backbone\_coeff.data, the "5" in the second row of the `[Rama]` block has been changed from "5" to "3". This has the effect of ignoring the last two rows in the Rama block, which correspond to the extra biases that would be applied if a secondary structure prediction was used.

It is important to note that simply "turning off" the `[SSWeight]` block  using "`[SSweight]`-" without making the above change to `[Rama]` WILL NOT turn off the bias as you might expect. Instead, the alpha AND beta bias will be applied to all residues. Both biases will also be applied to all residues if you change the 1s to 0s in the `[SSWeight]` block. Therefore it is not recommended to ever use `[SSWeight]-` or to delete the `[SSWeight]` block from the fix\_backbone\_coeff.data file or to have 0s in a column of the `[SSWeight]` block that correspond to rows of the `[Rama]` block which you don't want to apply to all residues. Incidentally, if you make the above change to the `[Rama]` block (change 5 to 3) and turn off `[SSWeight]` using `[SSWeight]-` then the extra alpha and beta biases will not be applied, but in that case using `[SSWeight]-` is unnecessary.

As an alternative to the above solution, you may choose to simply edit your ssweight file to contain two columns of all 0.0. This will turn off the bias for all residues, as you might expect.