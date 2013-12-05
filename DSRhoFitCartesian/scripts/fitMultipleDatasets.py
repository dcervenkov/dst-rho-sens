#!/usr/bin/python3

import sys
import os
import subprocess
import time
import cmath

# What do_fit values mean what
optFit = 1
optGenerate = 4

# Change these settings to configure the script
###############################################
max_jobs = 22
root_path = "../"
data_path = "data_gen5/"
results_path = "results/"
log_path = "logs/"
do_fit = optFit
do_plot = 0
enable_checks = 1 # Disable checks which require user input when running with, e.g., nohup
###############################################


def parseMetafile(dataset,pars):
    with open(dataset+".meta") as f:
        readingHelicityPars = False
        for line in f:
            if(line.strip() == ''): continue # Ignore empty lines
            if(line.find(':') == -1):
                print("ERROR: Cannot parse line " + str(linenum) + ": '"+ line.strip('\n') +"'")
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

def convertToSemiCartesian(pars):
    rp = float(pars[5])
    r0 = float(pars[6])
    rt = float(pars[7])
    sp = float(pars[8])
    s0 = float(pars[9])
    st = float(pars[10])
    crp = cmath.rect(rp,sp)
    cr0 = cmath.rect(r0,s0)
    crt = cmath.rect(rt,st)
    xp = crp.real
    yp = crp.imag
    x0 = cr0.real
    y0 = cr0.imag
    xt = crt.real
    yt = crt.imag
    pars[5] = "%.3f" % round(xp,3)
    pars[6] = "%.3f" % round(x0,3)
    pars[7] = "%.3f" % round(xt,3)
    pars[8] = "%.3f" % round(yp,3)
    pars[9] = "%.3f" % round(y0,3)
    pars[10]= "%.3f" % round(yt,3)
    
def convertToCartesian(pars):
    phiw = float(pars[4])
    rp = float(pars[5])
    r0 = float(pars[6])
    rt = float(pars[7])
    sp = float(pars[8])
    s0 = float(pars[9])
    st = float(pars[10])

    crp = cmath.rect(rp,-phiw+sp)
    cr0 = cmath.rect(r0,-phiw+s0)
    crt = cmath.rect(rt,-phiw+st)
    
    xp = crp.real
    yp = crp.imag
    x0 = cr0.real
    y0 = cr0.imag
    xt = crt.real
    yt = crt.imag
    
    crpb = cmath.rect(rp,+phiw+sp)
    cr0b = cmath.rect(r0,+phiw+s0)
    crtb = cmath.rect(rt,+phiw+st)
    
    xpb = crpb.real
    ypb = crpb.imag
    x0b = cr0b.real
    y0b = cr0b.imag
    xtb = crtb.real
    ytb = crtb.imag

    pars[4] = "%.3f" % round(xp,3)
    pars[5] = "%.3f" % round(x0,3)
    pars[6] = "%.3f" % round(xt,3)
    pars[7] = "%.3f" % round(yp,3)
    pars[8] = "%.3f" % round(y0,3)
    pars[9] = "%.3f" % round(yt,3)
    pars[10]= "%.3f" % round(xpb,3)
    pars[11]= "%.3f" % round(x0b,3)
    pars[12]= "%.3f" % round(xtb,3)
    pars[13]= "%.3f" % round(ypb,3)
    pars[14]= "%.3f" % round(y0b,3)
    pars[15]= "%.3f" % round(ytb,3)

def createCommand(dataset,pars):
    command = 'nice bin/Release/DSRhoFitCartesian ' + dataset + ' ' + dataset.replace(data_path,results_path) + \
              '.res ' + ' '.join(pars) + ' ' + str(do_fit) + ' ' + str(do_plot) +' > ' + dataset.replace(data_path,log_path) + '.log ' + '2>&1 &'
    return command

def numberOfJobs():
    stdout = subprocess.check_output('ps aux|grep bin/Release/DSRhoFitCartesian |wc -l',shell=True)
    jobs = int(stdout) - 2
    return jobs

def checkIfEmpty(path):
    if os.listdir(path) != [] :
        answer = input(path + " is not empty, continue? [y/N]: ")
        if answer != "y": quit(1)

def initChecks():
    if do_fit == optGenerate :
        checkIfEmpty(log_path)
        checkIfEmpty(data_path) #TODO Should check only for data files
        for dataset in datasets:
            if not os.path.isfile(dataset + ".meta"):
                print("ERROR: Not all metafiles exist in " + data_path)
                quit(2)
    elif do_fit == optFit :
        checkIfEmpty(results_path)
        checkIfEmpty(log_path)
        for dataset in datasets:
            if not os.path.isfile(dataset):
                print("ERROR: Not all datasets exist in " + data_path)
                quit(2)
    
    print("   max_jobs = %i" % max_jobs)
    print("   data_path = %s" % data_path)
    print("   results_path = %s" % results_path)
    print("   log_path = %s" % log_path)
    print("   do_fit = %s" % do_fit)
    print("   do_plot = %s" % do_plot)
    answer = input("Submit the jobs with this configuration? [y/N]: ")
    if answer != "y": quit(2)
    
    
if len(sys.argv) != 3 :
    print("ERROR: Wrong number of arguments!")
    print("Usage: fitMultipleDatasets.py beg_dataset_no end_dataset_no")
    quit(85)
elif not (sys.argv[1].isdigit() and sys.argv[2].isdigit()):
    print("ERROR: Argument(s) not positive integer(s)!")
    quit(86)

beg_dataset = int(sys.argv[1])
end_dataset = int(sys.argv[2])

os.chdir(root_path)
datasets = []
pars = [None]*16

for i in range(beg_dataset,end_dataset+1):
    datasets.append(data_path+"dataset_"+str(i))

if enable_checks: initChecks()

job_counter = 0
for dataset in datasets:
    parseMetafile(dataset,pars)
    convertToCartesian(pars)
    command = createCommand(dataset,pars)
    
    while(numberOfJobs() >= max_jobs):
        time.sleep(10)

    subprocess.call(command,shell=True)
    job_counter += 1
    print("Submitted job %i/%i, %i to go." % (job_counter, len(datasets), len(datasets)-job_counter))

print("All jobs have been submitted.")
