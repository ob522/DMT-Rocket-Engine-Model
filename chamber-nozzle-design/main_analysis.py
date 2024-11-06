# Main script to run the analysis
from data_loader import assign_constants
from fluid_functions import chamber_nozzle_sizing, conical_radius_function, conical_area_function, mach_no_as_function_of_distance, pressure_as_function_of_distance,temperature_as_function_of_distance,density_as_function_of_distance,velocity_as_function_of_distance
from heat_transfer_functions import hg_as_function_of_distance, constant_hl_value, SS_HT_analysis_along_chamber, SS_HT_itterative_analysis_along_chamber
from plotting import plot_chamber_nozzle_geometry,plot_multiple_sets
import constants as c  # Import constants from the constants file

#import matplotlib.pyplot as plt


def main():
    (Pe, P0, F, CR, char_L, T0, rho0, Cp, gamma, R, alpha) = assign_constants(c)

    (A0, At, Ae, Lc) = chamber_nozzle_sizing(Pe, P0, F, CR, char_L, T0, gamma, R)

    (x_values1, r_values) = conical_radius_function(A0, At, Ae, Lc, alpha)[:2]
    #(x_values, r_values) = conical_area_function(A0, At, Ae, Lc, alpha)

    (x_values2, M_values) = mach_no_as_function_of_distance(A0, At, Ae, Lc, alpha, gamma)
    (P_values) = pressure_as_function_of_distance(A0, At, Ae, Lc, alpha, gamma,P0)[1]
    (rho_values) = density_as_function_of_distance(A0, At, Ae, Lc, alpha, gamma,rho0)[1]
    (T_values) = temperature_as_function_of_distance(A0, At, Ae, Lc, alpha, gamma,T0)[1]
    (V_values) = velocity_as_function_of_distance(A0, At, Ae, Lc, alpha, gamma,R,T0)[1]

    A_values,throat_index = conical_area_function(A0, At, Ae, Lc, alpha)[1:]


    Tw = 2500 ################# this is the firts estimate of wall temp. ######## will be itterated !!!!!!!!!
    (hg_values) = hg_as_function_of_distance(P0,A0, At, Ae, Lc, alpha, gamma,R,Cp,T0,Tw)[1]
    print(hg_values[throat_index])
    print("........................")

    """
    Tw = 2500 ################# this is the firts estimate of wall temp. ######## will be itterated !!!!!!!!!
    (hg_values) = hg_as_function_of_distance(P0,A0, At, Ae, Lc, alpha, gamma,R,Cp,T0,Tw)[1]

    hl = constant_hl_value(3000, 1036*10**(-6), 0.17, 0.06, 9*10**(-6), 0.003)  ### test valuesS
    tw = 0.003
    Kw = 401
    q_values, Tl_values, Twg_values, Twl_values = SS_HT_analysis_along_chamber(x_values2, hg_values, T_values, tw, Kw, hl)
    """

    """"
    tw = 0.003
    Kw = 401
    hl = constant_hl_value(3000, 1036*10**(-6), 0.17, 0.06, 9*10**(-6), 0.003)  ### test values
    Twg_start_val = 2000
    q_values, Tl_values, Twg_values, Twl_values = SS_HT_itterative_analysis_along_chamber( A_values,M_values,T_values,tw,Kw,hl,   P0,At,gamma,R,Cp,T0,Twg_start_val)
    """




    print (x_values1)
    print (r_values)

    #Plots:

    #plot_chamber_nozzle_geometry(x_values1, r_values)

    plots = (
        #(x_values2, M_values, 'X Values', 'Mach Number'),
        #(x_values2, P_values, 'X Values', 'Pressure (P)'),
        #(x_values2, rho_values, 'X Values', 'Density (ρ)'),
        (x_values2, T_values, 'X Values', 'Temperature (T)'),
        #(x_values2, V_values, 'X Values', 'Velocity (V)'),
        (x_values2, hg_values, 'X Values', 'hg'),
        #(x_values2, q_values, 'X Values', 'q'),
        #(x_values2, Tl_values, 'X Values', 'Tl'),
        #(x_values2, Twg_values, 'X Values', 'Twg'),
        #(x_values2, Twl_values, 'X Values', 'Twl')

    )

    plot_multiple_sets(plots)


if __name__ == "__main__":
    main()
