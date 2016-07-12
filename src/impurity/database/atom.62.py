n=[5,6,7] # Atomic valences for oca
nc=[5,6,7] # Atomic valences for ctqmc
l=3       # Orbital angular momentum of the shel
J=0.8355  # Slatter integrals F2=J*11.9219653179191 from the atomic physics program
cx=0.112 # spin-orbit coupling from the atomic physics program

para=1    # Runs in paramagnetic mode

Ex=[0.5,0.5,0.5]  # Energy window treated exactly - relevant only in magnetic state
Ep=[3.0,3.0,3.0]  # Energy window treated approximately

qOCA=1    # OCA diagrams are computes in addition to NCA diagrams
Eoca=1.0  # Energy window for OCA diagrams
qatom=0   # Prints full atomic energis rather than E-E_{ls}
