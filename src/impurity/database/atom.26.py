n=[5,6,7] # Atomic valences
nc=[5,6,7] # Atomic valences                                                                                                                   
l=2           # Orbital angular momentum of the shel
J=0.79         # Slatter integrals F2=J*11.9219653179191 from the atomic physics program, must check this for Fe as I just used Haule's value for Mn here
cx=0.0 # 0.052        # spin-orbit coupling from the atomic physics program
qOCA=1        # OCA diagrams are computes in addition to NCA diagrams
Eoca=1.       # If the atomi energy of any state, involved in the diagram, is higher that Eoca from the ground state, the diagram ms is neglected
mOCA=1e-3     # If maxtrix element of OCA diagram is smaller, it is neglected 
Ncentral=[6]  # OCA diagrams are selected such that central occupancy is in Ncentral
Ex=[0.5,0.5,0.5]  # Energy window treated exactly - relevant only in magnetic state
Ep=[3.0,3.0,3.0]  # Energy window treated approximately
qatom=0

