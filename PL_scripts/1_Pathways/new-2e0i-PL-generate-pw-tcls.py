
acceptors = ["C4","O4","C4'","C2","C5","O4'","O2","C5'","C2'","N3","C6","O2'","H10","N1","N3'","H25","C1R","N1'","C5M","C5N","H13","C6'","H11","'C1*'","H9","H21","H23","H12","H24","H22"]
for acceptor in acceptors:
    query = '''
    mol load parm7 2e0i.top  dcd  2e0i-353K-trimmed-wrapped.dcd

    package require pathways

    pathtraj -d "resid 429  and name N1" -a "resname CPD and name {atom}" -b "resid 429 431 or {{within 8 of resid 429}} and {{within 8 of resname CPD}}" -withh "1" -p "1" -exph "1.146" -expts "1.146"
    pathtraj -d "resid 429  and name C2" -a "resname CPD and name {atom}" -b "resid 429 431 or {{within 8 of resid 429}} and {{within 8 of resname CPD}}" -withh "1" -p "1" -exph "1.146" -expts "1.146"
    pathtraj -d "resid 429  and name O2" -a "resname CPD and name {atom}" -b "resid 429 431 or {{within 8 of resid 429}} and {{within 8 of resname CPD}}" -withh "1" -p "1" -exph "1.146" -expts "1.146"
    pathtraj -d "resid 429  and name N3" -a "resname CPD and name {atom}" -b "resid 429 431 or {{within 8 of resid 429}} and {{within 8 of resname CPD}}" -withh "1" -p "1" -exph "1.146" -expts "1.146"
    pathtraj -d "resid 429  and name C4" -a "resname CPD and name {atom}" -b "resid 429 431 or {{within 8 of resid 429}} and {{within 8 of resname CPD}}" -withh "1" -p "1" -exph "1.146" -expts "1.146"
    pathtraj -d "resid 429  and name C4A" -a "resname CPD and name {atom}" -b "resid 429 431 or {{within 8 of resid 429}} and {{within 8 of resname CPD}}" -withh "1" -p "1" -exph "1.146" -expts "1.146"
    pathtraj -d "resid 429  and name C5A" -a "resname CPD and name {atom}" -b "resid 429 431 or {{within 8 of resid 429}} and {{within 8 of resname CPD}}" -withh "1" -p "1" -exph "1.146" -expts "1.146"
    pathtraj -d "resid 429  and name C6" -a "resname CPD and name {atom}" -b "resid 429 431 or {{within 8 of resid 429}} and {{within 8 of resname CPD}}" -withh "1" -p "1" -exph "1.146" -expts "1.146"
    pathtraj -d "resid 429  and name C7" -a "resname CPD and name {atom}" -b "resid 429 431 or {{within 8 of resid 429}} and {{within 8 of resname CPD}}" -withh "1" -p "1" -exph "1.146" -expts "1.146"
    pathtraj -d "resid 429  and name C7M" -a "resname CPD and name {atom}" -b "resid 429 431 or {{within 8 of resid 429}} and {{within 8 of resname CPD}}" -withh "1" -p "1" -exph "1.146" -expts "1.146"
    pathtraj -d "resid 429  and name C8" -a "resname CPD and name {atom}" -b "resid 429 431 or {{within 8 of resid 429}} and {{within 8 of resname CPD}}" -withh "1" -p "1" -exph "1.146" -expts "1.146"
    pathtraj -d "resid 429  and name C8M" -a "resname CPD and name {atom}" -b "resid 429 431 or {{within 8 of resid 429}} and {{within 8 of resname CPD}}" -withh "1" -p "1" -exph "1.146" -expts "1.146"
    pathtraj -d "resid 429  and name C9" -a "resname CPD and name {atom}" -b "resid 429 431 or {{within 8 of resid 429}} and {{within 8 of resname CPD}}" -withh "1" -p "1" -exph "1.146" -expts "1.146"
    pathtraj -d "resid 429  and name C9A" -a "resname CPD and name {atom}" -b "resid 429 431 or {{within 8 of resid 429}} and {{within 8 of resname CPD}}" -withh "1" -p "1" -exph "1.146" -expts "1.146"
    pathtraj -d "resid 429  and name C10" -a "resname CPD and name {atom}" -b "resid 429 431 or {{within 8 of resid 429}} and {{within 8 of resname CPD}}" -withh "1" -p "1" -exph "1.146" -expts "1.146"
    pathtraj -d "resid 429  and name N10" -a "resname CPD and name {atom}" -b "resid 429 431 or {{within 8 of resid 429}} and {{within 8 of resname CPD}}" -withh "1" -p "1" -exph "1.146" -expts "1.146"
    pathtraj -d "resid 429  and name H14" -a "resname CPD and name {atom}" -b "resid 429 431 or {{within 8 of resid 429}} and {{within 8 of resname CPD}}" -withh "1" -p "1" -exph "1.146" -expts "1.146"
    pathtraj -d "resid 429  and name H15" -a "resname CPD and name {atom}" -b "resid 429 431 or {{within 8 of resid 429}} and {{within 8 of resname CPD}}" -withh "1" -p "1" -exph "1.146" -expts "1.146"
    pathtraj -d "resid 429  and name H19" -a "resname CPD and name {atom}" -b "resid 429 431 or {{within 8 of resid 429}} and {{within 8 of resname CPD}}" -withh "1" -p "1" -exph "1.146" -expts "1.146"
    pathtraj -d "resid 429  and name H18" -a "resname CPD and name {atom}" -b "resid 429 431 or {{within 8 of resid 429}} and {{within 8 of resname CPD}}" -withh "1" -p "1" -exph "1.146" -expts "1.146"
    pathtraj -d "resid 429  and name O4" -a "resname CPD and name {atom}" -b "resid 429 431 or {{within 8 of resid 429}} and {{within 8 of resname CPD}}" -withh "1" -p "1" -exph "1.146" -expts "1.146"
    pathtraj -d "resid 429  and name N5" -a "resname CPD and name {atom}" -b "resid 429 431 or {{within 8 of resid 429}} and {{within 8 of resname CPD}}" -withh "1" -p "1" -exph "1.146" -expts "1.146"
    pathtraj -d "resid 429  and name H30" -a "resname CPD and name {atom}" -b "resid 429 431 or {{within 8 of resid 429}} and {{within 8 of resname CPD}}" -withh "1" -p "1" -exph "1.146" -expts "1.146"
    pathtraj -d "resid 429  and name H31" -a "resname CPD and name {atom}" -b "resid 429 431 or {{within 8 of resid 429}} and {{within 8 of resname CPD}}" -withh "1" -p "1" -exph "1.146" -expts "1.146"
    pathtraj -d "resid 429  and name 'C1*' " -a "resname CPD and name {atom}" -b "resid 429 431 or {{within 8 of resid 429}} and {{within 8 of resname CPD}}" -withh "1" -p "1" -exph "1.146" -expts "1.146"
    pathtraj -d "resid 429  and name H20" -a "resname CPD and name {atom}" -b "resid 429 431 or {{within 8 of resid 429}} and {{within 8 of resname CPD}}" -withh "1" -p "1" -exph "1.146" -expts "1.146"
    pathtraj -d "resid 429  and name H12" -a "resname CPD and name {atom}" -b "resid 429 431 or {{within 8 of resid 429}} and {{within 8 of resname CPD}}" -withh "1" -p "1" -exph "1.146" -expts "1.146"
    pathtraj -d "resid 429  and name H13" -a "resname CPD and name {atom}" -b "resid 429 431 or {{within 8 of resid 429}} and {{within 8 of resname CPD}}" -withh "1" -p "1" -exph "1.146" -expts "1.146"
    pathtraj -d "resid 429  and name H21" -a "resname CPD and name {atom}" -b "resid 429 431 or {{within 8 of resid 429}} and {{within 8 of resname CPD}}" -withh "1" -p "1" -exph "1.146" -expts "1.146"
    pathtraj -d "resid 429  and name H16" -a "resname CPD and name {atom}" -b "resid 429 431 or {{within 8 of resid 429}} and {{within 8 of resname CPD}}" -withh "1" -p "1" -exph "1.146" -expts "1.146"
    pathtraj -d "resid 429  and name H17" -a "resname CPD and name {atom}" -b "resid 429 431 or {{within 8 of resid 429}} and {{within 8 of resname CPD}}" -withh "1" -p "1" -exph "1.146" -expts "1.146"
    pathtraj -d "resid 429  and name 'C2*' " -a "resname CPD and name {atom}" -b "resid 429 431 or {{within 8 of resid 429}} and {{within 8 of resname CPD}}" -withh "1" -p "1" -exph "1.146" -expts "1.146"
    '''
    with open("%s.tcl" %acceptor, "w") as f:
        f.write(query.format( atom = acceptor))
