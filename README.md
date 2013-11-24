SingleTprime-analysis
=====================

This README file draw how to use this repository in order to do:

1. Run the extractor + analysis in a given dataset
2. How to access the code for the analysis
3. How to use the default macro to plot the Tprime mass

First you should set up your working area. For that please foolow the instructions here:

https://github.com/IPNL-CMS/PatTopProduction/blob/master/README.md

Then you should clone the repository and compile it in your area. For this you should do:

src> mkdir Extractors
src> cd Extractors
Extractors> git clone https://github.com/IPNL-CMS/SingleTprime-analysis.git
Extractors> cd ..
src> scram b -j8

Running on MC for the Signal:

In order to run the extractor + analysis on the signal you should do:

cd test/
cmsRun SingleTprime_analysis_cfg.py

With the flag OnlyAnalysis set to False. You can find this flag in line 4 of SingleTprime_analysis_cfg.py. This will create an extracted_mc.root file in the test directory.

After this you can run only the analysis. For this set OnlyAnalysis to True and redo 

cmsRun SingleTprime_analysis_cfg.py

This step will create an analyzed_mc.root file with only the information relevant to the analysis, and no information coming from the extractor.

In analyzed_mc.root you will fin one tree with two branches: stp and cuts. In the cuts branchs only the info of events passing or not each cut of the analysis is contained. Each leaf has the number of the cut with values o or 1 meaning that an event not passed or passed the respective cut.

The stp branch contains information about the analysis and the variables of interest. The most important is the TLorentzVector Reconstructed_Tprime which is the lorentz vector of the reconstructed T prime according to the analysis procedure after all cuts. For example, from this vector you can recover the mass that corresponds to the mass of the 5 jets system identified as coming from the decay of the T prime. In the same logic you can find a Reconstructed_W, Reconstructed_Higgs and Reconstructed_Top TLorentzVectors.

If the task involves a change in the analysis procedure, you can find the analysis code in plugins/SingleTprime_analysis.cc and plugins/SingleTprime_analysis.h

After having changed the code you can compile your new code with "scram b -j 8" command. Then you can rerun "cmsRun SingleTprime_analysis_cfg.py" in the analysis only mode beacuse you don't need to rerun the extractor.

In order to obtain automatic info from the analyzed MC from the signal you can go into the bin/ folder. In this folder you can find the GeneralAnalyzer.cc root macro. If you set lines 59 and 60 to

  int ParcialTestMax=NumberOfProcesses;
  int ParcialTestMin=NumberOfProcesses-1;

you will have printed in the screen the total number of events processed and the events passing each cut. Additionally, a plot of the mass of the 5 jets system.

Running on MC for the backgrounds:

To run over the backgrounds you have to use crab, as all datasets are not present in in2p3 T2. Then, for the moment, the only possibility is to run in the extractor + analysis mode.

To run over the backgrounds, 
