#ifndef DSELECTOR_HELPERFUNCS
#define DSELECTOR_HELPERFUNCS
// PUT THINGS HERE THAT ARE NEEDED BY THE DSELECTOR PROGRAM + OTHER PROGRAMS.
// FOR EXAMPLE FOR SIDEBAND SUBTRACTION WE NEED TO DEFINE THE REGIONS WHICH WILL BE USED TO CALCULATE THE SIDEBAND WEIGHTS BUT ALSO
// NEED TO THE REGIONS IN THE MAKEGRAPHS.C CODE TO DRAW THE LINES


//Variables for determining sideband subtraction regions.
double pi0Mean=0.135881;
double etaMean=0.548625;
double pi0Std=0.0076;
double etaStd=0.0191;
double pi0Sig=3;
double pi0Skip=1;
double etaSig=3;
double etaSkip=1;

//// Option 1 - simply from the phase space size 
//double pi0SB=2;
//double etaSB=2;
//double etaSBweight= -1.0*etaSig/etaSB;
//double pi0SBweight= -1.0*pi0Sig/pi0SB;

// Option 2 from fitting Mpi0 and Meta independently and extracting weights from fits 
double pi0SB=1;
double etaSB=1;
double etaSBweight= -2.92;
double pi0SBweight= -3.52;
#endif
