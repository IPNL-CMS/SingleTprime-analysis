SingleTprime-analysis
=====================

This README file draw how to use this repository in order to do:

1. Run the extractor + analysis in a given dataset
2. How to access the code for the analysis

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

To run over the backgrounds you have to use crab, as all datasets are not present in in2p3 T2. For the moment, the only possibility is to run in the extractor + analysis mode.

For this task you will find in the test/ folder the script createAndRunMCCrab.py. In between lines 20 to 31 you can add new datasets to run on the extracto. All this datasets are already pat tuples. Additionally, you can comment out a dataset, for example if you have already run over it, putting a # at the beginning of the corresponding line.

With all datasets you want to process with the extractor added, just type from the test/ folder:

./createAndRunMCCrab.py --create-cfg

to create the crab configuration files,

./createAndRunMCCrab.py --run

to launch all the crab jobs,

./createAndRunMCCrab.py --status

to get the status of all jobs, and

./createAndRunMCCrab.py --get

to get the output of all crab jobs.

With this done you should find the extractor results in the in2p3 T2 under the following location:

/dpm/in2p3.fr/home/cms/data/store/user/jruizalv/jruizalv/Extracted_MC/

Search a subdirectory named under the date you launched the production and a directory under with the name you inserted createAndRunMCCrab.py for the dataset. For example for TTbar, if you launched a crab submission the 3 of december you should find the results of the extractor in the following folder:

/dpm/in2p3.fr/home/cms/data/store/user/jruizalv/jruizalv/Extracted_MC/3Dec13/TTJets/

Remember the commands rfdir and frcp commands to list and copy from the T2, respectively.

The extractor files are not very heavy, so you copy them into /gridgroup/cms/ area to have a local access. There is a script to do this task, test/CopyExtractorResulst.py

