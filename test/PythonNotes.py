import commands
a=commands.getoutput('ls /gridgroup/cms/jruizalv/Extracted_MC/Wjets_VBF/')
r=a.split('\n')
for i in r: print '"' + i + '", ',
print len(r)
