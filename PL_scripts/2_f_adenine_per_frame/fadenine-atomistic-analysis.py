import os
import numpy as np
import re
import sys
np.set_printoptions(threshold=np.nan)
numbers = re.compile(r'(\d+)')
def numericalSort(value):
    parts = numbers.split(value)
    parts[1::2] = map(int, parts[1::2])
    return parts
prot_system = sys.argv[1]
prot = prot_system[0]+prot_system[1]+prot_system[2]+prot_system[3]
outfile_modifier = '-f-atomistic-vs-frame-matrix'
outfile_s1 = prot_system + outfile_modifier + '-s1.log'
outfile_s2 = prot_system + outfile_modifier + '-s2.log'
path = '/home/bjr29/new-A2-PL-nextstep/' + prot_system + '/logs'
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

first=int(0) #indices start at 0, this doesn't mean I can start at some value after 0 sadly.
lastplusone=int(6251) 
numframes=lastplusone-first
numpaths=float(numframes)
pathways = []
fadenineS1 = []
fadenineS2 = []
estfrac = []
###start acceptor for loop


adenineatoms = ("H9\n", "H10\n", "N71\n", "N61\n", "H32\n", "C81\n", "N91\n", "C51\n", "C61\n", "C41\n", "N11\n", "N31\n", "C21\n", "H11\n")
num_adenine_atoms = len(adenineatoms)
running_weighted_matrix_adenine_atomistic_involvement_S1 = np.array([[0.0 for cols in range(num_adenine_atoms)] for rows in range(numframes)])
running_weighted_matrix_adenine_atomistic_involvement_S2 = np.array([[0.0 for cols in range(num_adenine_atoms)] for rows in range(numframes)])
for acceptor in acceptors:
    aindex=acceptors.index(acceptor)
    filename = prot + '_%s' %acceptor
    var1='egrep "glob|Trajectory|FAD|1:" <%s.log> aux1.txt' %filename
    os.system(var1)
    var2 = '''awk '{print $1 " " $3 " " $5}' aux1.txt > aux2.txt'''
    os.system(var2)
    os.system('egrep -v "resname|within|net" <aux2.txt> aux3.txt') #aux3.txt contains frame and FAD atoms involved
    os.system('csplit --digits=2 --elide-empty-files --quiet --prefix=outfile aux3.txt "/pattern/+1" "{*}"')
    for item in os.listdir(path): 
        if item.startswith("outfile"):
            var4 = 'grep -v pattern <%s> %s.dut' %(item, item)
            os.system(var4)
            var41 = 'tail -n +2 %s.dut > %s.dat' %(item, item)
            os.system(var41)
    outfilescounts = []
    legend = []
    for item in os.listdir(path): #sort the filenames (the outfiles, containing the singleDonor-singleAcceptor-allFrames blocks) in slegend so the ordering matches donorsS1 (which is equal to donorsS2)
        if item.startswith("outfile"):
            if item.endswith(".dat"):
                legend.append(item)
    slegend=sorted(legend, key=numericalSort)
    for donor in slegend:#now take a singleAcceptor-singleDonor-allFrames block one at a time and in the order as listed in donorsS1
        dindex=slegend.index(donor)
        with open(donor, 'a+') as f:#now take the tidy DA block and split it into all the frames with csplit. The directory will have N_acceptor logs, 32 outfiles, and N_frames files titled framexxxx
            framelegend = []
            binarray = []
            weightarrayS1 = []
            weightarrayS2 = []
            varinf = 'csplit --digits=4 --elide-empty-files --quiet --prefix=frame %s "/Trajectory/" "{*}"' %(donor)#split the current DA block (the outfilexx.dat specified by the donor loop) into frames
            os.system(varinf)
            for fitem in os.listdir(path):
                if fitem.startswith("frame"):
                    framelegend.append(fitem)
            sframelegend=sorted(framelegend, key=numericalSort)#sorted framexxxx filenames are in sframelegend
            frame_iter=0
            for framefile in sframelegend[first:lastplusone]:#go through the framexxxx files one at a time, in order of increasing frame number
                with open(framefile, 'a+') as ff:
                    for line in ff:
                        ade_atom_iter=0#set to 0 here so that I can loop through all adenine atoms for each line in the path of a single D-A-frame combination.
                        for atomname in adenineatoms:
                            if (atomname in line):
                                running_weighted_matrix_adenine_atomistic_involvement_S1[frame_iter][ade_atom_iter] += donorfracS1[dindex]*acceptorfrac[aindex]
                                running_weighted_matrix_adenine_atomistic_involvement_S2[frame_iter][ade_atom_iter] += donorfracS2[dindex]*acceptorfrac[aindex]
                            ade_atom_iter += 1
                frame_iter += 1
                binarray.append(frameval)#by the end of this for loop, you have appended a 0 or 1 to the binarray, which will end up being numframes long.
                frameval=0
            weightarrayS1 = np.array(binarray)*donorfracS1[dindex]*acceptorfrac[aindex] #weight the array here before switching to a new donor.
            weightarrayS2 = np.array(binarray)*donorfracS2[dindex]*acceptorfrac[aindex]
            runningweightarrayS1 = runningweightarrayS1 + weightarrayS1
            runningweightarrayS2 = runningweightarrayS2 + weightarrayS2
            os.system('rm frame*') #clear out the frames of the given outfile as we move on to the next of the 32 outfiles. As we change the outfile, we switch to a new Donor.
    os.system('find . -type f ! -name "*.log" -delete')#once we're done with an Acceptor, clear all the temp files, leaving just the logs, and repeat for the next log.


running_weighted_matrix_adenine_atomistic_involvement_S1.tofile(outfile_s1,sep=",",format="%1.15f")
running_weighted_matrix_adenine_atomistic_involvement_S2.tofile(outfile_s2,sep=",",format="%1.15f")

#sumS1 = 0
#sumS2 = 0
#sumS1 = np.mean(runningweightarrayS1)
#sumS2 = np.mean(runningweightarrayS2)
print "number of frames is %d" %numpaths
print '********************'
print 'output S1 adenine-atomistic fadenine for each frame (ROWS = frames, COLS = adenine atoms in the order: H9, H10, N71, N61, H32, C81, N91, C51, C61, C41, N11, N31, C21, H11)'
print 'SEE ' + outfile_s1 + ' for f_atomistic matrix of dimensions (num_frames, num_adenine_atoms) = ' + repr(running_weighted_matrix_adenine_atomistic_involvement_S1.shape)
#print running_weighted_matrix_adenine_atomistic_involvement_S1
print '********************'
print 'output S2 adenine-atomistic fadenine for each frame (ROWS = frames, COLS = adenine atoms in the order: H9, H10, N71, N61, H32, C81, N91, C51, C61, C41, N11, N31, C21, H11)'
print 'SEE ' + outfile_s2 + ' for f_atomistic matrix of dimensions (num_frames, num_adenine_atoms) = ' + repr(running_weighted_matrix_adenine_atomistic_involvement_S2.shape)
#print running_weighted_matrix_adenine_atomistic_involvement_S2
