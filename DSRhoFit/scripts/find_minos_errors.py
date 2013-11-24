#!/usr/bin/python3

import sys

def parseLogFile(logfile,errors,error_names):
    with open(logfile) as f:
        for line in f:
            if(line.find('Minos:') == -1): continue
            if(line.find(':',34) == -1): continue

            line = line.replace('Minos: ','')

            isLower = False
            isUpper = False

            if line.find('Lower') != -1: 
                isLower = True
                line = line.replace('Lower error for parameter ','')
            if line.find('Upper') != -1: 
                isUpper = True
                line = line.replace('Upper error for parameter ','')

            if isLower and isUpper:
                print("ERROR: both upper and lower")
                quit(2)
            if not isLower and not isUpper:
                print("ERROR: neither upper nor lower")
                quit(2)

            tag, val = line.split(':',1)
            tag = tag.strip()
            val = val.strip(' \n-')
            
            for i in range(len(error_names)):
                if tag == error_names[i]:
                    if isLower: errors[i*2] = val
                    else: errors[i*2+1] = val


error_names = ['phiw','rp','r0','rt','sp','s0','st']
errors = [0] * (len(error_names)*2)
outputfile = sys.argv[1]

with open(outputfile, 'w') as f:
#    for name in error_names:
#        f.write('%s- %s+ ' % (name, name))
#    f.write('\n')
    for file in sys.argv[2:] :
        parseLogFile(file,errors,error_names)
        f.write(' '.join(errors)+'\n')
#print(errors)
