import math
import sys
import subprocess

PI = 3.1415

ap_min = 0.15	#0.1
ap_max = 0.35	#0.4
ap_step = 0.1
ap = ap_min
ap_steps = math.ceil((ap_max-ap_min)/ap_step + 0.00001)

apa_min = 0.1	    #0
apa_max = 1.1	    #2*PI
apa_step = 0.15
apa = apa_min
apa_steps = math.ceil((apa_max-apa_min)/apa_step + 0.00001)

a0_min = 0.85	#0.8
a0_max = 0.95	#1
a0_step = 0.1
a0 = a0_min
a0_steps = math.ceil((a0_max-a0_min)/a0_step + 0.00001)

ata_min = 2.3	    #0
ata_max = 3.3	    #2*PI
ata_step = 0.15
ata = ata_min
ata_steps = math.ceil((ata_max-ata_min)/ata_step + 0.00001)

numFits = ap_steps*apa_steps*a0_steps*ata_steps
print "Will perform %i fits" % numFits
print "ap : from %0.2f to %0.2f in %i steps of %0.2f" % (ap_min, ap_max, ap_steps, ap_step)
print "apa: from %0.2f to %0.2f in %i steps of %0.2f (%0.0f deg)" % (apa_min, apa_max, apa_steps, apa_step, apa_step*180)
print "a0 : from %0.2f to %0.2f in %i steps of %0.2f" % (a0_min, a0_max, a0_steps, a0_step)
print "ata: from %0.2f to %0.2f in %i steps of %0.2f (%0.0f deg)" % (ata_min, ata_max, ata_steps, ata_step, ata_step*180)
print "Should I proceed? (y,n)"
answer = raw_input()

if answer != "y":
    sys.exit(1)

while ap <= ap_max:
    while apa <= apa_max:
	while a0 <= a0_max:
	    while ata <= ata_max:
		output = str(ap)+"-"+str(apa)+"-"+str(a0)+"-"+str(ata)
		#subprocess.call(["echo","bsub -q g ./DSRhoFit","data/dataset_trans_org_1","results/"+output,str(ap),str(apa),str(a0),str(ata)])
                subprocess.call(["bsub -q g ./DSRhoFit data/dataset_trans_org_1 results/"+output+" "+str(ap)+" "+str(apa)+" "+str(a0)+" "+str(ata)],shell=True)
		ata += ata_step
	    a0 += a0_step
	    ata = ata_min
	apa += apa_step
	a0 = a0_min
    ap += ap_step
    apa = apa_min
