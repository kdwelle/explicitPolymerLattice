# explicitPolymerLattice
c++ code for running ion conductivity simulations using an explicit polymer on a lattice

The code is broken up into a few different classes:
- lattice: stores all the information about the simulation box. positions, energies, distances between ions etc.
- polymer: uses a self-avoiding-random walk to lay down a polymer and a variety of moves to sample configurational space.

