import subprocess

start_job = 45
start_val = 0.0
num_jobs = 17
step = 0.1

for i in range(num_jobs):
    subprocess.call(["../bin/Release/DSRhoFit ../data/dataset_" + str(i+start_job) + " ../results/result 0.107 1.42 0.941 0.31 1.7937 0.1 0.1 0.1 0 " +str(i*step+start_val)+ " 0 0 1"],shell=True)    
    subprocess.call(["mv ../plots/nll_phiw.png ../plots_nll/nll_phiw_" +str(i*step+start_val)+ ".png"],shell=True)
    subprocess.call(["mv ../data/dataset_from_pdf ../data/dataset_from_pdf_r_0.1_s0_" +str(i*step+start_val)],shell=True)

