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
donorfracS1 = np.array([-0.17952, -0.11294, 0.11138, 0.17978, -0.04544, -0.22115, -0.06983, 0.15051, -0.18502, -0.05564, 0.07419, 0.02578, 0.21570, -0.07621, 0.23639, 0.13716, -0.39710, 0.40502, 0.37660, -0.33851, 0.09831, -0.08143, 0.01235, 0.01905, 0.01517, 0.07386, -0.01219, -0.02325, -0.02903, 0.00672, 0.01163, -0.00858])
donorfracS2 = np.array([0.02433, 0.06247, -0.06286, -0.05354, -0.03486, -0.02802, 0.23592, 0.48174, -0.06222, -0.02343, -0.35726, -0.22048, 0.31350, -0.29772, 0.03222, 0.11668, 0.32729, -0.14080, -0.10230, 0.31203, 0.00433, 0.24002, -0.00257, 0.08272, -0.11830, -0.04357, 0.02127, -0.17812, 0.21137, 0.03562, 0.07559, 0.01044])
acceptors = ["C4", "O4", "C4P", "C2", "C5", "O4P", "O2", "C5P", "C2P", "N3", "C6", "O2P", "H10", "N1", "N3P", "H25", "C1R", "N1P", "C5M", "C5N", "H13", "C6P", "H11", "C1S", "H9", "H21", "H23", "H12", "H24", "H22"]
acceptorfrac = np.array([-0.09801, -0.18604, 0.36294, 0.01224, -0.13189, -0.12859, 0.01090, 0.35921, -0.22202, -0.08401, 0.20405, 0.06618, 0.24853, -0.00663, -0.14914, -0.17992, -0.24923, 0.08454, -0.07485, 0.07713, -0.08168, 0.49777, -0.08938, -0.03688, -0.04271, -0.03536, -0.07353, -0.07829, 0.04466, 0.01065])
#acceptors = ["C4", "O4"]
#acceptorfrac = np.array([0.25593, 0.12197])
#acceptors = ["C4"]
#acceptorfrac = np.array([0.25593])

donorslistS1 = dict( (key, value) for (key, value) in zip(donorsS1, donorfracS1) )
donorslistS2 = dict( (key, value) for (key, value) in zip(donorsS2, donorfracS2) )
acceptorslist = dict( (key, value) for (key, value) in zip(acceptors, acceptorfrac) )

first=int(0) #indices start at 0
lastplusone=int(6251) #one over target, as this value is the first index you don't need: as in, it will be excluded from the loop. As such, if I wanted to loop the first 1000 frames, I'd do first=0 lastplusone=1000. conveniently, the index of the frames starts at 0, so first will always match the frame index you want.
numpaths_int=lastplusone-first
numpaths=float(numpaths_int)
runningweightarrayS1 = np.zeros(numpaths_int)
runningweightarrayS2 = np.zeros(numpaths_int)
pathways = []
fadenineS1 = []
fadenineS2 = []
estfrac = []
###start acceptor for loop
for acceptor in acceptors:
    aindex=acceptors.index(acceptor)
    filename = '1dnp_%s' %acceptor
    var1='egrep "PATH 1|pattern" <%s.log> aux1.txt' %filename
    os.system(var1)
    var2 = '''awk '{print $5}' aux1.txt > aux2.txt'''
    os.system(var2)
    os.system('csplit --digits=2 --elide-empty-files --quiet --prefix=outfile aux2.txt "/pattern/+1" "{*}"')      #csplit might not be on every computer. Blocks the log into smaller logs. each log is PW from a new donor to the acceptor being studied.
    for item in os.listdir(path): 
        if item.startswith("outfile"):
            var4 = 'grep -v pattern <%s> %s.dat' %(item, item)
            os.system(var4)
    outfilescounts = []
    legend = []
    for item in os.listdir(path):
        if item.startswith("outfile"):
            if item.endswith(".dat"):
                legend.append(item)
    slegend=sorted(legend, key=numericalSort)
    for nitem in slegend:
        dindex=slegend.index(nitem)
        with open(nitem, 'r') as f:
            couplings = f.readlines()
            decimalcoups = np.array(couplings, dtype=float)
            weightarrayS1 = decimalcoups*donorfracS1[dindex]*acceptorfrac[aindex] #weight the array here before switching to a new donor.
            weightarrayS2 = decimalcoups*donorfracS2[dindex]*acceptorfrac[aindex]
#this is where the tallying is done. At this point, we have 2 arrays, each numframes long, containing the collective fadenine for each frame. As such, a single frame's fadenine won't be complete until we loop over all 960 Donor-Acceptor pairs. To get a single fadenine, we will later average the fadenine of each frame by summing the final array, and then dividing by the number of elements in the array (which, as a sanity check, we can show also equals numframes).
            runningweightarrayS1 = runningweightarrayS1 + weightarrayS1
            runningweightarrayS2 = runningweightarrayS2 + weightarrayS2
    os.system('find . -type f ! -name "*.log" -delete')

meanweightcoupS1 = np.mean(runningweightarrayS1)
meanweightcoupS2 = np.mean(runningweightarrayS2)

print "number of paths is %d" %numpaths
print '********************'
print 'S1 mean weighted coupling'
print meanweightcoupS1
print '********************'
print 'S2 mean weighted coupling'
print meanweightcoupS2
print '********************'
print 'S1 weighted coupling for each frame'
print runningweightarrayS1
print '********************'
print 'S2 weighted coupling for each frame'
print runningweightarrayS2
print '********************'

