foreach prot {1dnp-293K 1dnp-310K 1qnf-293K 1qnf-306K 1iqr-333K 2e0i-353K} { 
    cd /home/bjr29/new-A2-PL-PW/$prot
    mol load parm7 ${prot}.top
    animate read dcd ${prot}-trimmed-wrapped.dcd beg 1042 waitfor all
    set N6ind [atomselect top "resname FAD and name N61 and same residue as within 5 of resname CPD"]
    set C7ind [atomselect top "resname FAD and name C7 and same residue as within 5 of resname CPD"]
    set N10ind [atomselect top "resname FAD and name N10 and same residue as within 5 of resname CPD"]
    set O4ind [atomselect top "resname CPD and name O4"]
    set O4pind [atomselect top "resname CPD and name O4'"]
    set filename "${prot}-A2-distances-post500ns.txt"
    set fileId [open $filename "a"]
    foreach fadtag "$N6ind $C7ind $N10ind" {
        foreach cpdtag "$O4ind $O4pind" {
            set ent1 "$prot: [$fadtag get name] - [$cpdtag get name] distances"
            puts $fileId $ent1 
            puts $fileId [measure bond "[$fadtag get index] [$cpdtag get index]" frame all]
        }
    }
    close $fileId
    mol delete top
}
exit
