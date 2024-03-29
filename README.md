# Python scripts for comprehensive ET analysis of protein complexes  
This repository contains the tools to perform a comprehensive analysis of electron transfer (ET) in protein complex systems within the framework of the tunneling pathways model. As written, these scripts are for the photolyase-CPD (PL-CPD) complex.

Any questions about the code can be directed to Ben Rousseau (b.g.rousseau@gmail.com) or Dr. Agostino Migliore (agostino.migliore@unipd.it).

Citation: Rousseau, B. J. G.;  Migliore, A.;  Stanley, R. J.; Beratan, D. N. Adenine Fine-Tunes DNA Photolyase’s Repair Mechanism. J. Phys. Chem. B 2023, 127, 2941-2954.


Link: https://pubs.acs.org/doi/full/10.1021/acs.jpcb.3c00566

Here we present the scripts necessary to evaluate the extent of adenine’s involvement in the forward ET repair process from flavin adenine dinucleotide (FADH–) to cyclobutane pyrimidine dimer (CPT) in the E. *coli* PL-CPD complex simulated at 310 Kelvin. 

# Environment
Unless otherwise stated, python2.7 scripts were run with python2.7.16, sys 2.7.16, python3 scripts were run with python3.7.3, sys 3.7.3, statistics and operator from Python3.7 library, matplotlib.pyplot from Python3 library, and all python scripts used numpy 1.16.2, re 2.2.1, csv 1.0. with the Debian GNU/Linux 10 (buster) OS, in the Linux 4.19.0-20-amd64 distribution with x86-64 architecture.
Pathways analysis and distance calculation were performed with VMD/1.9.2 and the Pathways 1.2 plugin [I.A. Balabin, X. Hu, D.N. Beratan, J. Comput. Chem. 2012, 33, 906-910. DOI:10.1002/jcc.22927].

# Pathways analysis:
First, pathways analysis  is performed on a simulation trajectory to extract the strongest electronic coupling paths between each donor atom on FAD and each acceptor atom on CPD. 

To generate the input files used for pathways analysis, run 
```
python3 new-PL-generate-pw-tcls.py 
```
In the case of S. tokodaii (PDBID 2E0I), which contains a second FAD molecule as its antenna cofactor, the above script will not work. Instead, run 
```
python3 new-2e0i-PL-generate-pw-tcls.py
``` 
Please note that the resid of the CPD and FADH– must be updated if you are not using the structures we provide in our previous work [B.J.G. Rousseau, S. Shafei, A. Migliore, R.J. Stanley, D.N. Beratan, J. Am. Chem. Soc. 2018, 140, 8, 2853–2861. DOI: 10.1021/jacs.7b11926].

After generating the pathways input files, run the pathways analyses with
```
sbatch 1dnp-pw.sh # change to whatever job controller your computer cluster uses
``` 
The resultant log files will be processed in the subsequent sections.

# f_adenine analysis:
## *f_adenine_calculation*
On the pathways log files, extract the value of f per frame (stored in *-fadenine-per-frame.log) with (note python 2.7 is used)
```
python2.7 1dnp-fadenine-per-frame.py
```
The value of f per frame for each individual adenine atom is obtained by running
```
python2.7 fadenine-atomistic-analysis.py proteintemp > log.txt
```
Where proteintemp is the identifier of your system (refer to fadenine-atomistic-analysis.py to see how you should name your files and to change the directory to your pathway analysis log files).
The resultant fatomistic values will be stored in f-atomistic-vs-frame-matrix-s#.log (#=1 or 2).

## *f_adenine_stratification*
The following Python code is used to identify the MD snapshots belonging to a same *f* stratum. Beginning with the snapshot with the largest *f* value, we adopted *∆f=0.01* as the threshold of acceptance for the next snapshots. A stratum was identified only if it contained a minimum of 250 MD snapshots (this number was scaled down for the TT system, proportionally to the shorter MD production run examined). The first discarded snapshot was used to start a next stratum search. The protocol used avoided any ambiguity in assigning a snapshot to different strata.

To perform the stratification of the f_adenine data, run
```
python3 stratification.py
```
#### Descriptions of the outputs of stratification.py:
*-f-atomistic-s#-strata.csv
Contains matrices organized in descending order of the \<f\><sub>stratum</sub>, one matrix for each stratum, tagged by the stratum’s average f.  In a given stratum’s matrix, the columns are the adenine atoms, and the rows are the frames identified to belong to that stratum.

*-f-atomistic-s#-strata-averages.csv contains the average f<sub>atomistic</sub> of each adenine atom for each stratum identified (average f value of each stratum provided in first column), and is therefore a matrix where the columns are the adenine atoms, and the rows are the strata identified. 

*-f-full-strata-s#.png are images of the strata. 


# Additional scripts:
calculate_N6C7N10_O_distances.tcl calculates the distances between the set of oxygen atoms on CPD {O4,O4'} and the set of atoms {N6,C7,N10} on FAD (see main paper figure 1) over the MD trajectory. It was run in vmd/1.9.2 with 
```
vmd -e calculate_N6C7N10_O_distances.tcl -dispdev text
```
Electronic couplings (eq 7a in citation) (from which you can calculate mean-square tunneling pathway products and values of the coherence parameter C) for states S1 and S2 were calculated using the 1dnp-310K-getweightedcouplings-cc.py script, which was run as-included on the output log files of the pathways analyses:
```
python2.7 1dnp-310K-getweightedcouplings-cc.py
```
