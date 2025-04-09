module SMJ_Constants

# Physical constants
const e = 1.0 # charge of an electron in AU
const hbar = 1.0 # reduced Planck's constant in AU
const h = 2 * π # Planck's constant in AU
const m_e = 1.0 # mass of an electron in AU
const I_SI = 6.6236182375082e-3 # convert current in AU to current in Ampere

# Conversion factors
const Angstrom = 1.0 / 0.529177210903 # 1 Angstrom in AU
const eV = 1.0 / 27.21138624598853 # 1 eV in AU
const joule = 1.0 / 4.359744722207185e-18 # 1 joule in AU
const kB = 1.380649e-23 * joule # Boltzmann constant in AU/K

end # module SMJ_Constants
########## EFFECTIVE MASSES COMPUTED BY ROY GROUP via DOI: https://doi.org/10.1103/PhysRevB.65.245105   ########
#0.13 effective mass of the electron in BENZENE--TETRACENE-&-PENTACENE in AU; 0.11 E-Mass-NAPTHALENE; 
#0.12 E-Mass-ANTHRACENE, 0.26 effective mass of electron in ALKANE in AU, 0.013 effective mass in ALKENE in AU
