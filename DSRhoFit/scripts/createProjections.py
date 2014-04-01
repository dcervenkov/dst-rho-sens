#!/usr/bin/python3

import sys
import os
import subprocess
import time

# Change these settings to configure the script
###############################################
root_path = "../"
data_path = "data/"
projection_type = "evt"
results_path = "results_evt3/"
log_path = "logs/"
plots_path = "plots/"
do_fit = 0
do_plot = 1
plot_fit_results = 1 # Whether to plot fit results or generator values
enable_checks = 1 # Disable checks which require user input when running with, e.g., nohup
###############################################

pars = ['','','','','','','','','','','']

def parseResultFile(dataset,pars):
    init_value = 5 + plot_fit_results
    with open(dataset.replace(data_path,results_path) + ".res") as f:
        for line in f:
            if(line.strip() == ''): continue # Ignore empty lines

            pars2=[]
            vals = line.split(' ')
            for i in range(13):
                if i in (3,4): continue
                pars2.append(vals[init_value+i*4])

            for i in range(0,11):
                pars[i] = pars2[i]

def createCommand(dataset,pars):
    command = 'nice bin/Release/DSRhoFit ' + dataset + ' ' + dataset.replace(data_path,results_path) + \
              '.res ' + ' '.join(pars) + ' ' + str(do_fit) + ' ' + str(do_plot) +' > ' + dataset.replace(data_path,log_path) + '.plotlog ' + '2>&1 '
    return command

def checkIfEmpty(path):
    if os.listdir(path) != [] :
        answer = input(path + " is not empty, continue? [y/N]: ")
        if answer != "y": quit(1)

if len(sys.argv) != 3 :
    print("ERROR: Wrong number of arguments!")
    print("Usage: " + os.path.basename(__file__) + " beg_dataset_no end_dataset_no")
    quit(85)
elif not (sys.argv[1].isdigit() and sys.argv[2].isdigit()):
    print("ERROR: Argument(s) not positive integer(s)!")
    quit(86)

beg_dataset = int(sys.argv[1])
end_dataset = int(sys.argv[2])

os.chdir(root_path)
datasets = []

if enable_checks: checkIfEmpty(plots_path)

for i in range(beg_dataset,end_dataset+1):
    datasets.append(data_path+"dataset_"+str(i))

job_counter = 0
for dataset in datasets:
    parseResultFile(dataset,pars)
    command = createCommand(dataset,pars)
    job_counter += 1
    print("Submitted job %i/%i, %i to go." % (job_counter, len(datasets), len(datasets)-job_counter))

    subprocess.call(command,shell=True)
    projection_file = dataset.replace(data_path,plots_path).replace("dataset_","projections_"+projection_type) + ".root"
    os.rename(plots_path + "projections.root", projection_file)

print("All jobs have been submitted.")
