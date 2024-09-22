# Fit to the $^{113}$ Cd $\beta^{-}$ spectrum to extract the effective value of $g_{A}$.
## Author: Toby Dixon (toby.dixon@universite-paris-saclay.fr)
## Date  : June 1 2023

## Step 0: Setup the enviroment

Use:
	source /software/LoadBat.sh
This will load BAT, root etc.	

## Step 1: Convolution of theory with MC

The first step is a macro `ConvolveMC.cxx`
This will create the input file for BAT, the options are precised on the command line if compile with:

     g++ macros/ConvolveMC.cxx -o ConvolveMC `root-config  --cflags --glibs`

This will create an output file of `{OUTPATH}/{NAME}_{RESOTYPE}_{MODEL}.root` used for the gA fit.


Note that the input file containing the resolution, bias etc must have been updated to include the range of efficiency curve parameters.




## Step 2: gA fit

Inside:

Next you run the fit in gA_fit subdirectory, compile with make and then run:
     ./rungA_fit --help
to see options.
