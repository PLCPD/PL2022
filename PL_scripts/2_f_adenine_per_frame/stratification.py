import os
import numpy as np
import re
import statistics as stat
import matplotlib.pyplot as plt
import operator
import sys
import csv
np.set_printoptions(threshold=np.inf)
numbers = re.compile(r'(\d+)')
def numericalSort(value):
    parts = numbers.split(value)
    parts[1::2] = map(int, parts[1::2])
    return parts
for prot_temp in ['1dnp-293K', '1dnp-310K', '1iqr-333K', '1qnf-293K', '1qnf-306K', '2e0i-353K']:
    ##### READ IN MESSY ORIGINAL f_full LOG FILE #####
    #infile_f_full = sys.argv[1]
    infile_f_full = prot_temp + '-fadenine-per-frame.log'
    with open(infile_f_full, 'r') as fin:
        data = fin.read().splitlines(True)
    with open('testa.txt', 'w') as fout:
        fout.writelines(data[9:])
    os.system('''grep -v "*" testa.txt > testb.txt''')
    os.system('''csplit --elide-empty-files --quiet testb.txt "/output/"''')
    os.system('''grep -v "output" xx01 > xx02''')
    os.system('''mv xx00 f_s1.txt''')
    os.system('''mv xx02 f_s2.txt''')

    with open('f_s1.txt') as f:
        f_s1 = f.readlines()
    f_s1 = [x.strip() for x in f_s1]
    f_s1 = [x.strip(' []') for x in f_s1]
    full_f_list_s1 = [float(y) for x in f_s1 for y in x.split('  ')]

    with open('f_s2.txt') as f:
        f_s2 = f.readlines()
    f_s2 = [x.strip() for x in f_s2]
    f_s2 = [x.strip(' []') for x in f_s2]
    full_f_list_s2 = [float(y) for x in f_s2 for y in x.split('  ')]
    os.system('''rm testa.txt testb.txt f_s1.txt f_s2.txt xx01''')
    ##### END READ IN MESSY ORIGINAL f_full LOG FILE #####
    readin_f_atomistic_s1 = str(prot_temp) + '-f-atomistic-vs-frame-matrix-s1.log'
    f_atomistic_vector_s1 = np.fromfile(readin_f_atomistic_s1,sep=",")
    f_atomistic_matrix_s1 = np.reshape(f_atomistic_vector_s1,(len(f_atomistic_vector_s1)//14,14))
    readin_f_atomistic_s2 = str(prot_temp) + '-f-atomistic-vs-frame-matrix-s2.log'
    f_atomistic_vector_s2 = np.fromfile(readin_f_atomistic_s2,sep=",")
    f_atomistic_matrix_s2 = np.reshape(f_atomistic_vector_s2,(len(f_atomistic_vector_s2)//14,14))
    enumerated_full_f_list_s1 = list(enumerate(full_f_list_s1))
    enumerated_full_f_list_s2 = list(enumerate(full_f_list_s2))
    enumerated_full_f_list_s1.sort(key=lambda x:x[1], reverse=True)
    enumerated_full_f_list_s2.sort(key=lambda x:x[1], reverse=True)

    #Start descending down the list in order of decreasing f value, creating a list for every singificant jump in f value. 
    #How to define a significant jump? Two criteria, only one need be satisfied: 
    #1) a jump in magnitude greater than a threshold (that I will define as around the size of the smallest spacing separating two clear strata); 
    fixed_incr_thresh = 0.01
    preset_min_list_size_orig = 250 #* scale for systems with different simulation lengths
    preset_min_list_size = preset_min_list_size_orig
#    iqr_multiplier=.2272#because analysis was done on 500-1068 ns for 1IQR, vs 500-3000 ns for all others
#    if prot_temp in ['1iqr-333K']:
#        preset_min_list_size=round(preset_min_list_size_orig*iqr_multiplier)
#    else:
#        preset_min_list_size=preset_min_list_size_orig

    #2) NOT USED:if more than a preset number (30?) of data are added to the same strata (so in descending order there are 30 data that do not cause a jump that    #triggers criterion 1), then a jump that's larger than some pre-determined multiple (could be less than 1, could be greater) of the std deviation of that stratum.
    stddev_thresh_mult = 2.5 # I took out std dev element as it was yielding uncontrollable behavior. Turns out, cranking down on the incr thresh left some great strata! 0.01 or even 0.005 are excellent.
    ##### STRATIFY THE f_full, YIELDING BOTH STRATA AND FRAME TAGS FOR ATOMISTIC f ANALYSIS #####
    f_full_vs_stratum_f_avg_dict_s1 = {}#for plotting strata (visual verification of stratification)
    frame_vs_stratum_f_avg_dict_s1 = {}#for tagging frames (used then to sort the atomistic f values into strata in the next loop)
    temp_list_f = []
    temp_list_frame = []
    datum_count = 0
    datum_jump = 0.0
    stddev_temp = 10000
    for frame_idx, datum in enumerated_full_f_list_s1:
        if datum_count == 0:#new temp_list
            temp_list_f.append(datum)
            temp_list_frame.append(frame_idx)
        else:#temp_list_f already has entry/entries
            #datum_jump = abs(datum - temp_list_f[-1])#OLD WAY. FORM1. GOOD WITH SMALL INCR, BUT PRONE TO DRIFT, AS ADJACENT DATA SLOWLY DECREASE.
            datum_jump = abs(datum - temp_list_f[0])#NEW WAY. ALWAYS REFERENCE TO FIRST DATUM.
            #stddev_temp_old = stat.stdev(temp_list_f)
            if(len(temp_list_f)>30):
                #stddev_temp = stat.stdev(temp_list_f)#get stddev of temp_list_f before inclusion of datum
                thresh_used_this_iter = min(fixed_incr_thresh, stddev_thresh_mult*stddev_temp)#pick the minimum of the two metrics for the threshold
            else: 
                thresh_used_this_iter = fixed_incr_thresh
            if(datum_jump <= thresh_used_this_iter):#datum jump was lower than smallest threshold
                temp_list_f.append(datum)
                temp_list_frame.append(frame_idx)
                datum_count = datum_count+1
                continue #go on to test next datum.
            else:#datum was further than the smallest threshold from the last point on the ongoing list. Create new list with datum and terminate old list.
                if(len(temp_list_f) >= preset_min_list_size):#store the list if terminating list is larger than preset minimum
                    avg_to_store = stat.mean(temp_list_f)
                    f_full_vs_stratum_f_avg_dict_s1[str(avg_to_store)] = temp_list_f.copy()#store identified stratum with a label of its frame-averaged f_adenine
                    frame_vs_stratum_f_avg_dict_s1[str(avg_to_store)] = temp_list_frame.copy()
                temp_list_f = []#with list now stored if it was big enough, reset temp_list_f and populate 1st element with datum.
                temp_list_frame = []
                temp_list_f.append(datum)
                temp_list_frame.append(frame_idx)
        datum_count = datum_count+1
    #print(f_full_vs_stratum_f_avg_dict_s1.items())

    f_full_vs_stratum_f_avg_dict_s2 = {}#for plotting strata (visual verification of stratification)
    frame_vs_stratum_f_avg_dict_s2 = {}#for tagging frames (used then to sort the atomistic f values into strata in the next loop)
    temp_list_f = []
    temp_list_frame = []
    datum_count = 0
    datum_jump = 0.0
    stddev_temp = 10000
    for frame_idx, datum in enumerated_full_f_list_s2:
        if datum_count == 0:#new temp_list
            temp_list_f.append(datum)
            temp_list_frame.append(frame_idx)
        else:#temp_list_f already has entry/entries
            #datum_jump = abs(datum - temp_list_f[-1])#OLD WAY. FORM1. GOOD WITH SMALL INCR, BUT PRONE TO DRIFT, AS ADJACENT DATA SLOWLY DECREASE.
            datum_jump = abs(datum - temp_list_f[0])#NEW WAY. ALWAYS REFERENCE TO FIRST DATUM.
            #stddev_temp_old = stat.stdev(temp_list_f)
            if(len(temp_list_f)>30):
                #stddev_temp = stat.stdev(temp_list_f)#get stddev of temp_list_f before inclusion of datum
                thresh_used_this_iter = min(fixed_incr_thresh, stddev_thresh_mult*stddev_temp)#pick the minimum of the two metrics for the threshold
            else: 
                thresh_used_this_iter = fixed_incr_thresh
            if(datum_jump <= thresh_used_this_iter):#datum jump was lower than smallest threshold
                temp_list_f.append(datum)
                temp_list_frame.append(frame_idx)
                datum_count = datum_count+1
                continue #go on to test next datum.
            else:#datum was further than the smallest threshold from the last point on the ongoing list. Create new list with datum and terminate old list.
                if(len(temp_list_f) >= preset_min_list_size):#store the list if terminating list is larger than preset minimum
                    avg_to_store = stat.mean(temp_list_f)
                    f_full_vs_stratum_f_avg_dict_s2[str(avg_to_store)] = temp_list_f.copy()#store identified stratum with a label of its frame-averaged f_adenine
                    frame_vs_stratum_f_avg_dict_s2[str(avg_to_store)] = temp_list_frame.copy()
                temp_list_f = []#with list now stored if it was big enough, reset temp_list_f and populate 1st element with datum.
                temp_list_frame = []
                temp_list_f.append(datum)
                temp_list_frame.append(frame_idx)
        datum_count = datum_count+1


    ##### SAVE THE STRATA (f_full) AS IMAGES FOR VISUAL INSPECTION #####
    filename_stratified_f_full_plot_s1 = readin_f_atomistic_s1.replace('.log', '').replace('atomistic-vs-frame-matrix-','full-strata-')
    for k, v in f_full_vs_stratum_f_avg_dict_s1.items():#python3 uses items, python2 uses iteritems 
        plt.ylim([0,1])
        plt.xlabel('frame index (within each stratum)')
        plt.ylabel('f_adenine_s1')
        xs = range(0, len(v))
        ys = v
        plt.plot(xs, ys, '-.')
    plt.savefig(filename_stratified_f_full_plot_s1+'.png')
    plt.cla()
    filename_stratified_f_full_plot_s2 = readin_f_atomistic_s2.replace('.log', '').replace('atomistic-vs-frame-matrix-','full-strata-')
    for k2, v2 in f_full_vs_stratum_f_avg_dict_s2.items():#python3 uses items, python2 uses iteritems 
        plt.ylim([0,1])
        plt.xlabel('frame index (within each stratum)')
        plt.ylabel('f_adenine_s2')
        xs = range(0, len(v2))
        ys = v2
        plt.plot(xs, ys, '-.')
    plt.savefig(filename_stratified_f_full_plot_s2+'.png')
    plt.cla()
    ##### SAVE THE ATOMISTIC f FOR EACH STRATUM IN INDIVIDUAL CSV FILES TO BE READ AS MATRICES INTO EXCEL #####
    filename_base_stratified_f_atomistic_s1 = readin_f_atomistic_s1.replace('.log', '').replace('vs-frame-matrix-','') + '-stratum-'
    stratum_count = 0

    for f_avg, frame_list in frame_vs_stratum_f_avg_dict_s1.items():
        mean_list = np.zeros(14)
        temp_mat = np.zeros((len(frame_list)+1,14))#The last row of temp_mat contains column- (f_atomistic strata-)averages
        frame_count = 0
        for frame in frame_list:#inside a given stratum's frame list now
            for i in range(0,14):
                temp_mat[frame_count][i] = f_atomistic_matrix_s1[frame][i]
                mean_list[i] = mean_list[i] + temp_mat[frame_count][i]
            frame_count = frame_count + 1
        for i in range(0,14):
            temp_mat[frame_count][i] = mean_list[i]/frame_count 
        f_avg_flt=float(f_avg)
        listofstrings = [filename_base_stratified_f_atomistic_s1, str(stratum_count), str("favg-"), str(round(f_avg_flt,6)), ".csv"]
        fname_s1 = "".join(listofstrings)
        np.savetxt(fname_s1, temp_mat, delimiter=",")
        stratum_count = stratum_count + 1

    filename_base_stratified_f_atomistic_s2 = readin_f_atomistic_s2.replace('.log', '').replace('vs-frame-matrix-','') + '-stratum-'
    stratum_count = 0
    for f_avg_s2, frame_list_s2 in frame_vs_stratum_f_avg_dict_s2.items():
        mean_list = np.zeros(14)
        temp_mat = np.zeros((len(frame_list_s2)+1,14))#The last row of temp_mat contains column- (f_atomistic strata-)averages
        frame_count = 0
        for frame in frame_list_s2:#inside a given stratum's frame list now
            for i in range(0,14):
                temp_mat[frame_count][i] = f_atomistic_matrix_s2[frame][i]
                mean_list[i] = mean_list[i] + temp_mat[frame_count][i]
            frame_count = frame_count + 1
        for i in range(0,14):
            temp_mat[frame_count][i] = mean_list[i]/frame_count 
        f_avg_s2_flt=float(f_avg_s2)
        np.savetxt(filename_base_stratified_f_atomistic_s2 + str(stratum_count) + 'favg-' + str(round(f_avg_s2_flt,6)) + '.csv', temp_mat, delimiter=",")
        stratum_count = stratum_count + 1
