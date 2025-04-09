# Conductance-SMJS
Julia based scripts to calculate conductance and current for Single-molecule junctions

To run this file you will need to have the SMJ_Constants.jl and TransmissionFunctions.jl files in your path or located in the working directory. 

SMJ_Constants.jl need not be edited and has some effective masses calculated using the complex extension of k-vectors via Tomfohr and Sankey. DOI: https://doi.org/10.1103/PhysRevB.65.245105

TransmissionFunctions.jl will require QuadGk, Cubature and Printf packages to be loaded in your environment. The user can edit line 58 of this to change the name of Conductance_OUTPUT.txt to a system specific name, i.e. octanedithiol_conductance.txt. Additionally, one can change the output for the I-V on line 91.

The main.jl fucntion will prompt the user for the effective mass of the electron (hole), the jucntions barrier width in Angstroms, the Fermi level which we approximate as the HOMO of the Electrode-molecule-Electrode system for hole transport junctions with electron donanting contacts, the transport level which is described in this paper "PCCP paper" and the user specified voltage. Additionally, the user can change the working directory path (working_path) to their current directory or one of their choosing.
