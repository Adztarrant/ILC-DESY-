# ILC-DESY- Finding efficency of reconstucted tracks

This software will find the efficency of reconstructed and marlin tracks for muons in the ILCsoft marlin framework.

#runing 

src/RECDiagnostics.cc and include/RECDiagnostics.h must first be built and installed. 

mysteer.xml is given as an example as to the structure your stearing file needs. The software requires a ILC input file that contains MCParticle, Tracks, ReconstructedParticles and LCRelation. Most self run simulations this will the be file endding in REC.slcio.

run in command line Marlin mysteer.xml

This will produce a root file. Opening this root file will allow access to all the histograms and efficency plots produced by the RECDiagnostics programe. 

#code 

To change the particle type function varible mcpdg should be changed to the relvant monte carlo particle code. In defult this is set to be a muon (13). 

If momentum above 20GeV is used then the max momentum in the efficency graph will need to be changed. 