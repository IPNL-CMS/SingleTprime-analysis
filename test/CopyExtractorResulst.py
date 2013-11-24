#! /usr/bin/env python

import os, shutil, subprocess
from optparse import OptionParser
import commands

parser = OptionParser()
parser.add_option("-p", "--path", dest="path", type="string", help="where to Extractor results where store by crab")
parser.add_option("-o", "--obje", dest="obje", type="string", help="Name of folder to copy")
#parser.add_option("-m", "--merge", dest="merge", type="string", help="Name of merged file")

#Usage example: ./CopyExtractorResulst.py --path /dpm/in2p3.fr/home/cms/data/store/user/jruizalv/jruizalv/Extracted_MC/17Nov13/T-t/ --obje T-t/ --merge extracted_mc_merged.root

(options, args) = parser.parse_args()

#if options.path is None or not os.path.isdir(options.path) or options.obje is None or not os.path.isdir(options.path):
if options.path is None or options.obje is None:
    parser.error("you must specify a valid path")

FullString=commands.getoutput('rfdir ' + options.path)
AllLines=FullString.split('\n')
FilesNames=[]

for i in xrange(len(AllLines)):
    FilesNames.append(AllLines[i].rpartition(' ')[-1])

commands.getoutput('mkdir /gridgroup/cms/jruizalv/Extracted_MC/' + options.obje)

for CuFile in FilesNames:
    print 'Copying ' + CuFile
    OutPut=commands.getoutput('rfcp ' + options.path + CuFile + ' /gridgroup/cms/jruizalv/Extracted_MC/' + options.obje)
    
#if len(FilesNames)<50:
#    FullStringToMerge="hadd /gridgroup/cms/jruizalv/Extracted_MC/" + options.obje + options.merge + " " 
#
#    for CuFile in FilesNames:
#        FullStringToMerge+=("/gridgroup/cms/jruizalv/Extracted_MC/" + options.obje + CuFile + " ")
#
#    commands.getoutput(FullStringToMerge)
#
#else:
#    for i in xrange(len(FilesNames)):
#        if i==0: FullStringToMerge="hadd /gridgroup/cms/jruizalv/Extracted_MC/" + options.obje + "extracted_mc_merged_1.root "
#        if i==50:
#            commands.getoutput(FullStringToMerge)
#            FullStringToMerge="hadd /gridgroup/cms/jruizalv/Extracted_MC/" + options.obje + "extracted_mc_merged_2.root "
#        if i==100:
#            commands.getoutput(FullStringToMerge)
#            FullStringToMerge="hadd /gridgroup/cms/jruizalv/Extracted_MC/" + options.obje + "extracted_mc_merged_3.root "
#        if i==150:
#            commands.getoutput(FullStringToMerge)
#            FullStringToMerge="hadd /gridgroup/cms/jruizalv/Extracted_MC/" + options.obje + "extracted_mc_merged_4.root "
#        if i==200:
#            commands.getoutput(FullStringToMerge)
#            FullStringToMerge="hadd /gridgroup/cms/jruizalv/Extracted_MC/" + options.obje + "extracted_mc_merged_5.root "
#        if i==250:
#            commands.getoutput(FullStringToMerge)
#            FullStringToMerge="hadd /gridgroup/cms/jruizalv/Extracted_MC/" + options.obje + "extracted_mc_merged_6.root "
#        if i==300:
#            commands.getoutput(FullStringToMerge)
#            FullStringToMerge="hadd /gridgroup/cms/jruizalv/Extracted_MC/" + options.obje + "extracted_mc_merged_7.root "
#        if i==350:
#            commands.getoutput(FullStringToMerge)
#            FullStringToMerge="hadd /gridgroup/cms/jruizalv/Extracted_MC/" + options.obje + "extracted_mc_merged_8.root "
#        if i==400:
#            commands.getoutput(FullStringToMerge)
#            FullStringToMerge="hadd /gridgroup/cms/jruizalv/Extracted_MC/" + options.obje + "extracted_mc_merged_9.root "
#        if i==450:
#            commands.getoutput(FullStringToMerge)
#            FullStringToMerge="hadd /gridgroup/cms/jruizalv/Extracted_MC/" + options.obje + "extracted_mc_merged_10.root "
#        FullStringToMerge+=("/gridgroup/cms/jruizalv/Extracted_MC/" + options.obje + FilesNames[i] + " ")
#        if i==500: commands.getoutput(FullStringToMerge)
#        
#for CuFile in FilesNames:
#    commands.getoutput("rm " + "/gridgroup/cms/jruizalv/Extracted_MC/" + options.obje + CuFile)
#

print 'All Done!'
