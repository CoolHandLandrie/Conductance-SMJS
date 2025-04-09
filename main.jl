include("SMJ_Constants.jl")
include("TransmissionFunctions.jl")

using .SMJ_Constants       # Import SMJ_Constants module
using .TransmissionFunctions # Import TransmissionFunctions module

function main()
    # User inputs provided here
    working_path = "/Users/path/to/working/directory/or/output/location"
    
    println("Enter effective mass of the electron (in terms of m_e): ")
    m_eff_input = parse(Float64, readline()) 
    m_eff = m_eff_input * SMJ_Constants.m_e
    
    println("Enter barrier width (in Angstroms): ")
    a_input_angstroms = parse(Float64, readline())
    a = a_input_angstroms * SMJ_Constants.Angstrom
    
    println("Enter Fermi level energy (in AU): ")
    E_Fermi = parse(Float64, readline())
    
    println("Enter transport energy level (in AU): ")
    E_Transport = parse(Float64, readline())
    
    println("Enter applied bias voltage (in volts): ")
    Voltage = parse(Float64, readline())
    
    # Derived values
    Temp = 298.15  # Room temperature in Kelvin
    E_kin = 0.5 * SMJ_Constants.kB * Temp # Total energy of the electron in AU (KE)
    PhiB = abs(E_Transport - E_Fermi)     # Height of the rectangular barrier in AU

    # Calculate transmission and current values
    Transmission_zerobias = transmission_zero_bias(m_eff, E_kin, PhiB, a)
    Transmission_lowbias, errorA = transmission_low_bias(m_eff, E_kin, PhiB, Voltage, a)
    Current_lowbias, errorB = calculate_current(m_eff, E_kin, PhiB, Voltage, a)

    # Print results to output files
    print_results(
        m_eff,
        a,
        E_Fermi,
        E_Transport,
        PhiB,
        Transmission_zerobias,
        Transmission_lowbias,
        errorA,
        Current_lowbias,
        errorB,
        Voltage,          # Ensure Voltage is passed here
        working_path
    )


println("Calculating I-V curve...")
    
# Call calculate_iv_curve to compute and save I-V data
calculate_iv_curve(
    m_eff,
    E_kin,
    PhiB,
    a,
    working_path;
    Voltage_increment=0.001,
    range_start=-500,
    range_end=501
)
end

main()
