#!/bin/bash -e
    sys=U1mer
    for pH in 08.00 09.00 10.00
    do
    cat > plumed_rew_${pH}.dat << EOF
    
    # Activate MOLINFO functionalities
    MOLINFO STRUCTURE=data/U1mer.pdb
    chi:  TORSION ATOMS=8,9,19,28
    puck: PUCKERING ATOMS=6,8,9,13,11
        
    # Activate well-tempered metadynamics in Chi and Zx
    metad: METAD ARG=chi,puck.Zx PACE=100000000 HEIGHT=0.0 BIASFACTOR=4 SIGMA=0.35,0.35 FILE=data/HILLS_pH${pH} GRID_MIN=-pi,-pi GRID_MAX=pi,pi RESTART=YES # <- this is the new stuff! 
    
    as: REWEIGHT_BIAS ARG=metad.bias

    hhchi:  HISTOGRAM ARG=chi STRIDE=1 GRID_MIN=-pi GRID_MAX=pi GRID_BIN=100 BANDWIDTH=0.1 LOGWEIGHTS=as 
    hhpuck: HISTOGRAM ARG=puck.Zx STRIDE=1 GRID_MIN=-pi GRID_MAX=pi GRID_BIN=100 BANDWIDTH=0.1 LOGWEIGHTS=as 

    # Convert histograms h(s) to free energies F(s) = -kBT * log(h(s))
    ffchi: CONVERT_TO_FES GRID=hhchi 
    ffpuck: CONVERT_TO_FES GRID=hhpuck 
    # Print out the free energies F(s) to file once the entire trajectory is processed
    DUMPGRID GRID=ffchi FILE=ffchi_${pH}.dat 
    DUMPGRID GRID=ffpuck FILE=ffpuck_${pH}.dat 

# Make sure you print out the 3 variables in the specified order at every step
PRINT ARG=chi,puck.Zx,metad.bias FILE=COLVAR_REWEIGHT_${pH} STRIDE=1 # <- also change this one! 
EOF

plumed driver --mf_xtc data/${sys}_pH${pH}.xtc --plumed plumed_rew_${pH}.dat --kt 2.5
done