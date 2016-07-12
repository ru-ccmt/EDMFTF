n=[3,4,5,6] # Atomic valences for oca
nc=[3,4,5,6] # Atomic valences for ctqmc
l=3         # Orbital angular momentum of the shel
J=0.604     # Slatter integrals F2=J*11.9219653179191 from the atomic physics program
cx=0.2096   # spin-orbit coupling from the atomic physics program

para=1      # Runs in paramagnetic mode

Ex=[0.5,0.5,0.5] # Energy window treated exactly - relevant only in magnetic states
Ep=[4.0,2.0,2.0] # Energy window treated approximately

qOCA=1     # OCA diagrams are computes in addition to NCA diagrams
Eoca=0.2   # Energy window for OCA diagrams
qatom=0    # Prints full atomic energis rather than E-E_{ls}
