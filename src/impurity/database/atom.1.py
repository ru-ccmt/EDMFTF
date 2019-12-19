n=[0,1,2] # Atomic valences
nc=[0,1,2] # Atomic valences for ctqmc
l=0           # Orbital angular momentum of the shel
J=1e-6        # Slater integrals F2=J*11.9219653179191 from the atomic physics program  (2009 Nov 19, J=0.3 originally, changed by chuckyee)
cx=0.0        # spin-orbit coupling from the atomic physics program
qOCA=0        # OCA diagrams are computes in addition to NCA diagrams
Eoca=1.       # If the atomi energy of any state, involved in the diagram, is higher that Eoca from the ground state, the diagram ms is neglected
mOCA=1e-3     # If maxtrix element of OCA diagram is smaller, it is neglected 
Ncentral=[1]  # OCA diagrams are selected such that central occupancy is in Ncentral
