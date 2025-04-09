module TransmissionFunctions

using QuadGK   # for 1D integration
using Cubature # for 2D integration
using Printf   # for @printf macro
include("SMJ_Constants.jl")
using .SMJ_Constants

# Export functions so they can be accessed directly
export transmission_zero_bias, transmission_low_bias, calculate_current, print_results, calculate_iv_curve

# Function to calculate transmission at zero bias voltage
function transmission_zero_bias(m_eff, E_kin, PhiB, a)
    k1 = sqrt(2.0 * m_eff * E_kin / (SMJ_Constants.hbar^2))
    kappa = sqrt(2.0 * m_eff * (PhiB - E_kin) / (SMJ_Constants.hbar^2))
    return 4.0 * k1^2 * kappa^2 / ((k1^2 + kappa^2)^2 * sinh(kappa * a)^2 + 4.0 * k1^2 * kappa^2)
end

# Function to calculate transmission at low bias voltage using integration
function transmission_low_bias(m_eff, E_kin, PhiB, Voltage, a)
    function t2(x)
        k1 = sqrt(2.0 * m_eff * E_kin / (SMJ_Constants.hbar^2))
        kappa_argument = 2.0 * m_eff * ((PhiB + Voltage * SMJ_Constants.eV * x / a) - E_kin) / (SMJ_Constants.hbar^2)
        kappa = sqrt(max(0.0, kappa_argument))
        k2_argument = 2.0 * m_eff * (E_kin - Voltage * SMJ_Constants.eV) / (SMJ_Constants.hbar^2)
        k2 = sqrt(max(0.0, k2_argument))
        return 4.0 * k1 * kappa^2 * k2 / ((k1^2 + kappa^2) * (k2^2 + kappa^2) * sinh(kappa * a)^2 + (k1 * kappa + k2 * kappa)^2)
    end
    result, errorA = quadgk(t2, 0, a)
    return result / a, errorA
end

# Function to calculate current at low bias voltage using integration
function calculate_current(m_eff, E_kin, PhiB, Voltage, a)
    function t3(xE)
        E, x = xE
        k1 = sqrt(2.0 * m_eff * E_kin / (SMJ_Constants.hbar^2))
        kappa_argument = 2.0 * m_eff * ((PhiB + E + Voltage * SMJ_Constants.eV * x / a) - E_kin) / (SMJ_Constants.hbar^2)
        kappa = sqrt(max(0.0, kappa_argument))
        k2_argument = 2.0 * m_eff * (E_kin - E - Voltage * SMJ_Constants.eV) / (SMJ_Constants.hbar^2)
        k2 = sqrt(max(0.0, k2_argument))
        return 4.0 * k1 * kappa^2 * k2 / ((k1^2 + kappa^2) * (k2^2 + kappa^2) * sinh(kappa*a)^2 + (k1*kappa + k2*kappa)^2)
    end
    result, errorB = hcubature(t3, [0, 0], [-Voltage*SMJ_Constants.eV, a])
    current_lowbias = (2*SMJ_Constants.e/SMJ_Constants.h) * result*SMJ_Constants.I_SI/a
    return current_lowbias, errorB
end

# Function to print results to output files
function print_results(m_eff, a, E_Fermi, E_Transport, PhiB,
                       Transmission_zerobias,
                       Transmission_lowbias,
                       errorA,
                       Current_lowbias,
                       errorB,
                       Voltage,
                       working_path)
    output_file_path = joinpath(working_path, "Conductance_OUTPUT.txt")
    open(output_file_path, "w") do io
        @printf(io, "E_Fermi = %.8e AU\n", E_Fermi)
        @printf(io, "E_Transport = %.8e AU\n", E_Transport)
        @printf(io, "Barrier height = %.8e eV\n", PhiB/SMJ_Constants.eV)
        @printf(io, "Barrier width = %.8e Angstrom\n", a/SMJ_Constants.Angstrom)
        @printf(io, "Effective mass = %.8e m_e\n", m_eff/SMJ_Constants.m_e)
        @printf(io, "Bias voltage = %.8f V\n\n", Voltage)
        @printf(io, "The tunneling coefficient at zero bias voltage is: %.8e G0\n", Transmission_zerobias)
        @printf(io, "The tunneling coefficient at low bias voltage is: %.8e G0\n", Transmission_lowbias)
        @printf(io, "The error for transmission at low bias is: %.8e\n", errorA)
        @printf(io, "The current at low bias voltage is: %.8e A\n", Current_lowbias)
        @printf(io, "The error in the current at low bias is: %.8e\n", errorB)
    end

    println("Files written to: $working_path")
end

# Export functions so they can be accessed directly
export calculate_iv_curve

# Function to calculate I-V curve and write results to file
function calculate_iv_curve(
    m_eff::Float64,
    E_kin::Float64,
    PhiB::Float64,
    a::Float64,
    working_path::String;
    Voltage_increment::Float64=0.001,
    range_start::Int=-500,
    range_end::Int=501
)
    # Path for output file
    iv_file_path = joinpath(working_path, "I-V-OUTPUT.txt")
    
    # Open file for writing
    open(iv_file_path, "w") do io
        for Voltage_counter in range_start:range_end # Loop over voltage range
            Voltage = Voltage_increment * -abs(Voltage_counter)
            
            # Define the t3 function for integration
            function t3(xE)
                E, x = xE
                k1 = sqrt(2.0 * m_eff * E_kin / (SMJ_Constants.hbar^2))
                kappa = sqrt(2.0 * m_eff * ((PhiB + E + Voltage * SMJ_Constants.eV * x / a) - E_kin) / (SMJ_Constants.hbar^2))
                k2 = sqrt(2.0 * m_eff * (E_kin - E - Voltage * SMJ_Constants.eV) / (SMJ_Constants.hbar^2))
                return 4.0 * k1 * kappa^2 * k2 / ((k1^2 + kappa^2) * (k2^2 + kappa^2) * sinh(kappa * a)^2 + (k1 * kappa + k2 * kappa)^2)
            end
            
            # Perform integration using hcubature
            result3, errorC = hcubature(t3, [0, 0], [-Voltage * SMJ_Constants.eV, a])
            
            # Calculate current at low bias voltage
            Current_lowbias = (2 * SMJ_Constants.e / SMJ_Constants.h) * result3 * SMJ_Constants.I_SI / a
            
            # Print voltage with 3 decimal places and current with 8 decimal places
            println(io, "$(@sprintf("%.3f", Voltage_increment*Voltage_counter)) $(@sprintf("%.8f", sign(Voltage_counter)*Current_lowbias*1e9))")
        end
    end
    
    println("I-V curve written to: $iv_file_path")
end

end # module TransmissionFunctions