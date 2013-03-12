#!/usr/bin/python

import sys
import os
import subprocess
import time

max_jobs = 13
root_path = "../"
data_path = "data/"
results_path = "results/"
log_path = "logs/"

pars = ['','','','','','','','','','','']

def parseMetafile(dataset,pars):
    with open(dataset+".meta") as f:
	readingHelicityPars = False
	for line in f:
	    if(line.strip() == ''): continue # Ignore empty lines
	    if(line.find(':') == -1):
		print "ERROR: Cannot parse line " + str(linenum) + ": '"+ line.strip('\n') +"'"
		quit()

	    tag, val = line.split(':',1)
	    val = val.strip()
	    
	    if tag == 'hp':   readingHelicityPars = True
	    elif tag == 'sm': readingHelicityPars = False

	    if readingHelicityPars: continue
	    
	    if   tag == 'ap':  pars[0] = val
	    elif tag == 'apa': pars[1] = val
	    elif tag == 'a0':  pars[2] = val
	    elif tag == 'ata': pars[3] = val
	    elif tag == 'phiw':pars[4] = val
	    elif tag == 'rp':  pars[5] = val
	    elif tag == 'r0':  pars[6] = val
	    elif tag == 'rt':  pars[7] = val
	    elif tag == 'sp':  pars[8] = val
	    elif tag == 's0':  pars[9] = val
	    elif tag == 'st':  pars[10] = val

def createCommand(dataset,pars):
    command = 'nice bin/Release/DSRhoFit ' + dataset + ' ' + dataset.replace(data_path,results_path) + \
	      '.res ' + ' '.join(pars) + ' 1 0' +' > ' + dataset.replace(data_path,log_path) + '.log ' + '2>&1 &'
    return command

def numberOfJobs():
    stdout = subprocess.check_output('ps aux|grep bin/Release/DSRhoFit |wc -l',shell=True)
    jobs = int(stdout) - 2
    return jobs



if len(sys.argv) != 3 :
    print "ERROR: Wrong number of arguments!"
    print "Usage: fitMultipleDatasets.py beg_dataset_no end_dataset_no"
    quit()
elif not (sys.argv[1].isdigit() and sys.argv[2].isdigit()):
    print "ERROR: Argument(s) not positive integer(s)!"
    quit()

beg_dataset = int(sys.argv[1])
end_dataset = int(sys.argv[2])

os.chdir(root_path)
datasets = []

for i in range(beg_dataset,end_dataset+1):
    datasets.append(data_path+"dataset_"+str(i))

for dataset in datasets:
    parseMetafile(dataset,pars)
    command = createCommand(dataset,pars)
    
    while(numberOfJobs() >= max_jobs):
	time.sleep(10)

    subprocess.call(command,shell=True)
