import commands
import os, time, sys

if len(sys.argv)<4 or len(sys.argv)>4:
    print "Usage: python FileMover.py Workplace Version FileCheck"
    sys.exit(2)

#Logical switch for different tasks
#Only Checking if all folders have a root file (Not moving files nor merging root files)
RootFileCheck=sys.argv[3]

Version=sys.argv[2]+"/"

Samples=["QCD_PT_170_300","QCD_PT_300_470","QCD_PT_470_600","TTJets","QCD_PT_600_800","QCD_PT_800_1000","Tbar-s","Tbar-tw","Tbar-t","T-s","T-tw","T-t"]
Path="/afs/cern.ch/work/j/jruizalv/private/Analysis/CMSSW_5_3_9_patch3/src/Extractors/Analyzed_Files/WithTrigger/" + Version
MergedSamplesPath="/afs/cern.ch/work/j/jruizalv/private/Analysis/CMSSW_5_3_9_patch3/src/Extractors/PatExtractor/bin/WithTrigger/" + Version

Folder=sys.argv[1]
a=commands.getoutput("ls " + Folder)    
SubFolders=a.split('\n')
NotRootFile=[]

if RootFileCheck=="False":
    for i in Samples:
        commands.getoutput("mkdir " + Path)
        commands.getoutput("mkdir " + MergedSamplesPath)
        print "Creating folder " + Path + i
        commands.getoutput("mkdir " + Path + i)

for i in SubFolders:
    Indic=0
    b=commands.getoutput("ls " + Folder + i)
    Files=b.split('\n')
    #print Files
    for j in Files:
        if ".root" in j:
            Indic+=1
            #print "Found root file in: " + i
            if RootFileCheck=="False":
                #print "Moving: " + Folder + i + "/" + j
                for k in Samples:
                    if k in i:
                        #print "Corresponding to " + Path + k + "/"
                        print "Doing: " + "mv " + Folder + i + "/" + j + " " + Path + k + "/"
                        commands.getoutput("mv " + Folder + i + "/" + j + " " + Path + k + "/")
    #print Indic
    if Indic==0: NotRootFile.append(i)

print "List of subfolders without a root file: ", NotRootFile

PWD=commands.getoutput("pwd")
if RootFileCheck=="True":
    if len(NotRootFile)!=0:
        for i in NotRootFile:
            print "Entering ",  Folder + i, " and reruning"
            os.chdir(Folder + i)
            commands.getoutput("cmsRun run_cfg.py")
            os.chdir(PWD)            

if RootFileCheck=="False":
    for i in Samples:
        print "Beginning merging of root files..."
        commands.getoutput(" hadd " + MergedSamplesPath + i + "_Full_analyzed.root " + Path + i + "/*")
