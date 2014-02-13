import commands
import os, time, sys

if len(sys.argv)<2 or len(sys.argv)>2:
    print "Usage: python RunMany_OnlyAnalysis.py Queue"
    sys.exit(2)

FirstPath="/afs/cern.ch/work/b/beaucero/public/JoseSamples/"
SecondPath="/afs/cern.ch/work/j/jruizalv/private/Analysis/CMSSW_5_3_9_patch3/src/Extractors/Samples/"
#DirFirstPath=["QCD_PT_170_300"]
DirFirstPath=["QCD_PT_170_300","QCD_PT_300_470","QCD_PT_470_600","TTJets"]
DirSecondPath=["QCD_PT_600_800","QCD_PT_800_1000","Tbar-s","Tbar-t","Tbar-tw","T-s","T-t","T-tw"]

dir="Workplace_" + str(time.strptime(time.ctime(time.time())).tm_hour) + "_"  + str(time.strptime(time.ctime(time.time())).tm_min) + "_"  + str(time.strptime(time.ctime(time.time())).tm_mday) + "_" + str(time.strptime(time.ctime(time.time())).tm_mon) + "_" + str(time.strptime(time.ctime(time.time())).tm_year)
commands.getoutput("mkdir " + dir)
PWD=commands.getoutput("pwd")

for i in DirFirstPath:
    a=commands.getoutput("ls " + FirstPath + i +"/")    
    r=a.split('\n')
    k=0
    for j in r:
        k+=1
        curdir=dir + "/Job_" + i + "_" + str(k)
        FILEIN=FirstPath + i + "/" + j
        FILEOUT="analyzed_mc_" + i + "_" + str(k) +".root"
        print "Running for file " + FILEIN
        commands.getoutput("mkdir " + curdir)
        commands.getoutput('cat Template_SingleTprime_analysis_cfg.py | sed -e"s#FILEIN#' + FILEIN + '#g" | sed -e"s#FILEOUT#' + FILEOUT + '#g" > ' + curdir + '/run_cfg.py')
        commands.getoutput("cp triggers_jets.xml " + curdir + "/")
        commands.getoutput("cp Extractor_MTT_ScaleFactors.py " + curdir + "/")
        commands.getoutput("cp Muon_ID_iso_Efficiencies_Run_2012ABCD_53X.pkl " + curdir + "/")
        commands.getoutput("cp Electron_scale_factors.json " + curdir + "/")
        commands.getoutput("export LS_SUBCWD=`pwd`/" + curdir)
        os.chdir(curdir)
        commands.getoutput("bsub -q "+ sys.argv[1] +" -o /tmp/Job_out < ../../batchScript.sh")
        os.chdir(PWD + "/")

for i in DirSecondPath:
    a=commands.getoutput("ls " + SecondPath + i +"/")    
    r=a.split('\n')
    k=0
    for j in r:
        k+=1
        curdir=dir + "/Job_" + i + "_" + str(k)
        FILEIN=SecondPath + i + "/" + j
        FILEOUT="analyzed_mc_" + i + "_" + str(k) +".root"
        print "Running for file " + FILEIN
        commands.getoutput("mkdir " + curdir)
        commands.getoutput('cat Template_SingleTprime_analysis_cfg.py | sed -e"s#FILEIN#' + FILEIN + '#g" | sed -e"s#FILEOUT#' + FILEOUT + '#g" > ' + curdir + '/run_cfg.py')
        commands.getoutput("cp triggers_jets.xml " + curdir + "/")        
        commands.getoutput("cp Extractor_MTT_ScaleFactors.py " + curdir + "/")
        commands.getoutput("cp Muon_ID_iso_Efficiencies_Run_2012ABCD_53X.pkl " + curdir + "/")
        commands.getoutput("cp Electron_scale_factors.json " + curdir + "/")
        commands.getoutput("export LS_SUBCWD=`pwd`/" + curdir)
        os.chdir(curdir)
        commands.getoutput("bsub -q "+ sys.argv[1] +" -o /tmp/Job_out < ../../batchScript.sh")
        os.chdir(PWD + "/")
