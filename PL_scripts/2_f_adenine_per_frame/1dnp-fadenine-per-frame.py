import os
import numpy as np
import re
np.set_printoptions(threshold=np.inf)
numbers = re.compile(r'(\d+)')
def numericalSort(value):
    parts = numbers.split(value)
    parts[1::2] = map(int, parts[1::2])
    return parts
path = '/home/bjr29/new-A2-PL-nextstep/1dnp-310K/logs' #change
os.chdir(path)

donorsS1 = ["N1", "C2", "O2", "N3", "C4", "C4A", "C5A", "C6", "C7", "C7M", "C8", "C8M", "C9", "C9A", "C10", "N10", "H14", "H15", "H19", "H18", "O4", "N5", "H30", "H31", "C1S", "H20", "H12", "H13", "H21", "H16", "H17", "C2S"]
donorsS2 = ["N1", "C2", "O2", "N3", "C4", "C4A", "C5A", "C6", "C7", "C7M", "C8", "C8M", "C9", "C9A", "C10", "N10", "H14", "H15", "H19", "H18", "O4", "N5", "H30", "H31", "C1S", "H20", "H12", "H13", "H21", "H16", "H17", "C2S"]
donorfracS1 = np.array([0.02344, 0.02357, 0.0119, 0.01363, 0.01372, 0.02928, 0.05162, 0.0842, 0.21643, 0.01125, 0.06349, 0.01089, 0.06977, 0.21086, 0.08412, 0.01295, 0.01629, 0.01803, 0.00895, 0.00678, 0.00599, 0.00417, 0.00218, 0.00173, 0.00128, 0.00108, 0.001, 0.00052, 0.00039, 0.00003, 0.00001, -0.00002])
donorfracS2 = np.array([0.00131, 0.00468, 0.00298, 0.00036, 0.00141, 0.00474, 0.22206, 0.19494, 0.01395, 0.00496, 0.16236, 0.00189, 0.22196, 0.03042, 0.01612, 0.0119, 0.00709, -0.00043, 0.00324, 0.01945, 0.00015, 0.03986, 0.00605, 0.00741, -0.00337, 0.00268, 0.00012, 0.00685, 0.01227, 0.00062, 0.00199, 0.0002])
acceptors = ["C4", "O4", "C4P", "C2", "C5", "O4P", "O2", "C5P", "C2P", "N3", "C6", "O2P", "H10", "N1", "N3P", "H25", "C1R", "N1P", "C5M", "C5N", "H13", "C6P", "H11", "C1S", "H9", "H21", "H23", "H12", "H24", "H22"]
acceptorfrac = np.array([0.25593, 0.12197, 0.11132, 0.10099, 0.07434, 0.06442, 0.05872, 0.05345, 0.03261, 0.02112, 0.02103, 0.01859, 0.01232, 0.01064, 0.00901, 0.00842, 0.00834, 0.00649, 0.00526, -0.00251, 0.00221, 0.00213, -0.00159, 0.0011, 0.00108, 0.00098, 0.00088, -0.00074, 0.00053, 0.00043])

donorslistS1 = dict( (key, value) for (key, value) in zip(donorsS1, donorfracS1) )
donorslistS2 = dict( (key, value) for (key, value) in zip(donorsS2, donorfracS2) )
acceptorslist = dict( (key, value) for (key, value) in zip(acceptors, acceptorfrac) )

first=int(0) #indices start at 0, don't change.
lastplusone=int(6251) #this should be 1 greater than number of frames 
numpaths_int=lastplusone-first
numpaths=float(numpaths_int)
runningweightarrayS1 = np.zeros(numpaths_int)
runningweightarrayS2 = np.zeros(numpaths_int)
pathways = []
fadenineS1 = []
fadenineS2 = []
estfrac = []
###start acceptor for loop
adenineatoms = ("H9\n", "H10\n", "N71\n", "N61\n", "H32\n", "C81\n", "N91\n", "C51\n", "C61\n", "C41\n", "N11\n", "N31\n", "C21\n", "H11\n")
for acceptor in acceptors:
    aindex=acceptors.index(acceptor)
    filename = '1dnp_%s' %acceptor
    var1='egrep "glob|Trajectory|FAD|1:" <%s.log> aux1.txt' %filename
    os.system(var1)
    var2 = '''awk '{print $1 " " $3 " " $5}' aux1.txt > aux2.txt'''
    os.system(var2)
    os.system('egrep -v "resname|within|net" <aux2.txt> aux3.txt')
    os.system('csplit --digits=2 --elide-empty-files --quiet --prefix=outfile aux3.txt "/pattern/+1" "{*}"')
    for item in os.listdir(path): 
        if item.startswith("outfile"):
            var4 = 'grep -v pattern <%s> %s.dut' %(item, item)
            os.system(var4)
            var41 = 'tail -n +2 %s.dut > %s.dat' %(item, item)
            os.system(var41)
    outfilescounts = []
    legend = []
    for item in os.listdir(path):
        if item.startswith("outfile"):
            if item.endswith(".dat"):
                legend.append(item)
    slegend=sorted(legend, key=numericalSort)
    for nitem in slegend:
        dindex=slegend.index(nitem)
        with open(nitem, 'a+') as f:
            framelegend = []
            binarray = []
            weightarrayS1 = []
            weightarrayS2 = []
            varinf = 'csplit --digits=4 --elide-empty-files --quiet --prefix=frame %s "/Trajectory/" "{*}"' %(nitem)
            os.system(varinf)
            for fitem in os.listdir(path):
                if fitem.startswith("frame"):
                    framelegend.append(fitem)
            sframelegend=sorted(framelegend, key=numericalSort)
            for lel in sframelegend[first:lastplusone]:
                with open(lel, 'a+') as ff:
                    for line in ff:
                        if any(s in line for s in adenineatoms):
                            frameval=1
                        else:
                            frameval=0
                binarray.append(frameval)
                frameval=0
            weightarrayS1 = np.array(binarray)*donorfracS1[dindex]*acceptorfrac[aindex] #weight the array here before switching to a new donor.
            weightarrayS2 = np.array(binarray)*donorfracS2[dindex]*acceptorfrac[aindex]
            runningweightarrayS1 = runningweightarrayS1 + weightarrayS1
            runningweightarrayS2 = runningweightarrayS2 + weightarrayS2
            os.system('rm frame*') #clear out the frames for the next of the 32 outfiles. As we change the outfile, we switch to a new Donor.
    os.system('find . -type f ! -name "*.log" -delete')

sumS1 = 0
sumS2 = 0
sumS1 = np.mean(runningweightarrayS1)
sumS2 = np.mean(runningweightarrayS2)
print "number of paths is %d" %numpaths
print '********************'
print 'output bulk S1 fadenine'
print sumS1
print '********************'
print 'output bulk S2 fadenine'
print sumS2
print '********************'
print 'output S1 fadenine for each frame'
print runningweightarrayS1
print '********************'
print 'output S2 fadenine for each frame'
print runningweightarrayS2
