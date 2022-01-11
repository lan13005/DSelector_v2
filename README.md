## DSelector + Accompanying Software for GlueX Analysis
DSelector is a program that is written in C++ using the object oriented program and library, ROOT. ROOT is the main program for analyzing large physics datasets. This DSelector is used to analyze terabytes of GlueX data in the gamma p -> eta pi0 p -> 4gamma p reaction channel, which is an especially interesting channel to look for exotic hybrid mesons. The main purpose is to perform selections and creating histograms to extract physical insight. Supporting scripts to convert histograms in root files to images and to prune root trees. 

### Analysis of double regge beam asymmetries
The double-Regge exchange process is an important to understand in order to search for exotic hybrid mesons in this reaction channel. A signature of exotic hybrid mesons is an asymmetry in the angular distribution of the decay of the etapi0 system. The double-regge exchange process is expected to contribute to the resulting signal and could produce such asymmetries. fitAsymmetryPlots.C takes in a flat tree (post-DSelector) and makes measurements of the beam asymmetry at the top-vertex of the double regge exchange process in the photoproduction of etapi0. Measurements of the beam asymmetry gives access to the naturality of the exchange and allows one to build models describing the potential asymmetry. 

### Analysis of b1->omegapi0->5gamma leakage
Previous studies have suggested that the main contributor to combinatorial background (from other channels) are from the b1->5g process. It is possible to determine the expected leakage after event selections are applied since the photoprodution cross sections for the b1 have already been measured at GlueX. Remaining b1 leakage into the etapi0 dataset can be determined and limits can be set.

