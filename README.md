Fit to the ^{113}Cd \beta^{-} spectrum to extract the effective value of g_{A}.
Author: Toby Dixon (toby.dixon@universite-paris-saclay.fr)
Date  : June 1 2023



The first step is a macro ConvolveMC.cxx
This will create the input file for BAT, the options are precised on the command line if compile with:
g++ macros/ConvolveMC.cxx -o ConvolveMC `root-config  --cflags --glibs`
This will create an output file of {OUTPATH}/{NAME}_{RESOTYPE}_{MODEL}.root used for the gA fit.

Next you run the fit in gA_fit subdirectory, compile with make and then run ./rungA_fit --help to see options.

Note that the input file containing the resolution, bias etc must have been updated to include the range of efficiency curve parameters.




