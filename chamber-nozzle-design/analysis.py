# Main script to run the analysis
from data_loader import assign_constants
from nozzle_functions import chamber_nozzle_sizing, conical_radius_function, conical_area_function,parabolic_radius_function,parabolic_area_function, mach_no_as_function_of_distance, pressure_as_function_of_distance,temperature_as_function_of_distance,density_as_function_of_distance,velocity_as_function_of_distance
from regen_channel_functions import channel_width_as_function_of_distance, channel_area_as_function_of_distance, Re_no_as_function_of_distance, get_rho_and_velocity_values
from heat_transfer_functions import hg_as_function_of_distance, constant_hl_value,calculate_hg, Sequential_Heat_Transfer_Analysis, Sequential_Heat_Transfer_Analysis_Extended, Sequential_Heat_Transfer_Analysis_Extended2, Sequential_Heat_Transfer_Analysis_Extended3, adjust_for_film_cooling
from plotting import plot_chamber_nozzle_geometry,plot_multiple_sets, chamber_nozzle_geometry, plot_on_same_axes, plot_temperature_data, plot_Ts_vs_x, plot_Ts_vs_x2, plot_and_save_multiple_datasets, plot_multiple_datasets_varied_x, plot_multiple_datasets_varied_x_2, chamber_nozzle_geometry_2
from data_handler import write_arrays_to_file, filter_and_extract_numbers, split_data_columns, split_tests_by_time, find_closest_index
from stress_analysis import stress_analysis_calculations
import constants as c  # Import constants from the constants file
from external_heat_transfer import calc_external_wall_temp, lumped_capacitance, get_channel_CHT_coef, transient_conduction_1D

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os


def design_GW_analysis():
    (Pe, P0, F, CR, char_L, T0, rho0, Cp, gamma, R, alpha) = assign_constants(c)

    (A0, At, Ae, Lc, m_dot) = chamber_nozzle_sizing(Pe, P0, F, CR, char_L, T0, gamma, R)

    (x_values, r_values, throat_index) = conical_radius_function(A0, At, Ae, Lc, alpha)
    A_values = conical_area_function(A0, At, Ae, Lc, alpha)[1]

    M_values = mach_no_as_function_of_distance(A0, At, Ae, Lc, alpha, gamma)[1]
    P_values = pressure_as_function_of_distance(A0, At, Ae, Lc, alpha, gamma,P0)[1]
    rho_values = density_as_function_of_distance(A0, At, Ae, Lc, alpha, gamma,rho0)[1]
    T_values = temperature_as_function_of_distance(A0, At, Ae, Lc, alpha, gamma,T0)[1]
    V_values = velocity_as_function_of_distance(A0, At, Ae, Lc, alpha, gamma,R,T0)[1]
    #hg_values = hg_as_function_of_distance(P0,A0, At, Ae, Lc, alpha, gamma,R,Cp,T0,1000)[1]

    #Tw = 2500 ################# this is the firts estimate of wall temp. ######## will be itterated !!!!!!!!!
    #(hg_values) = hg_as_function_of_distance(P0,A0, At, Ae, Lc, alpha, gamma,R,Cp,T0,Tw)[1]
    #print(hg_values[throat_index])
    #print("........................")

    """
    Tw = 2500 ################# this is the firts estimate of wall temp. ######## will be itterated !!!!!!!!!
    (hg_values) = hg_as_function_of_distance(P0,A0, At, Ae, Lc, alpha, gamma,R,Cp,T0,Tw)[1]

    hl = constant_hl_value(3000, 1036*10**(-6), 0.17, 0.06, 9*10**(-6), 0.003)  ### test valuesS
    tw = 0.003
    Kw = 401
    q_values, Tl_values, Twg_values, Twl_values = SS_HT_analysis_along_chamber(x_values2, hg_values, T_values, tw, Kw, hl)
    """

    
    #tw = 0.003
    #Kw = 401
    hl = constant_hl_value(3000, 1036*10**(-6), 0.17, 0.06, 9*10**(-6), 0.003)  ### test values
    #Twg_start_val = 2000
    #q_values, Tl_values, Twg_values, Twl_values = SS_HT_itterative_analysis_along_chamber( A_values,M_values,T_values,tw,Kw,hl,   P0,At,gamma,R,Cp,T0,Twg_start_val)
    




    #print (x_values1)
    #print (r_values)

    #Plots:

    #plot_chamber_nozzle_geometry(x_values, r_values)

    chamber_nozzle_geometry(x_values, r_values, [(1, "Rc"),(throat_index, "Rt"),(-1, "Re")], x_label='Axial Position (x)', r_label='Radius (r)', title='Rocket Engine Thrust Chamber Geometry')
    
    plots = (
        
        #(x_values, A_values, 'X Values', 'Area'),
        #(x_values, M_values, 'X Values', 'Mach Number'),
        #(x_values, P_values, 'X Values', 'Pressure (P)'),
        #(x_values, rho_values, 'X Values', 'Density (ρ)'),
        (x_values, T_values, 'X Values', 'Temperature (T)'),
        #(x_values, V_values, 'X Values', 'Velocity (V)'),

        #(x_values, hg_values, 'X Values', 'hg')
        #(x_values2, q_values, 'X Values', 'q'),
        #(x_values2, Tl_values, 'X Values', 'Tl'),
        #(x_values2, Twg_values, 'X Values', 'Twg'),
        #(x_values2, Twl_values, 'X Values', 'Twl')

    )
    

    plot_multiple_sets(plots)
    
    print("At:", A_values[throat_index])
    print("D0",r_values[1]*2)
    print("Dt",r_values[throat_index]*2)
    print("De",r_values[-1]*2)
    print("Mt:", M_values[throat_index])
    print("Pt:", P_values[throat_index])
    print("rhot:", rho_values[throat_index])
    print("Tt:", T_values[throat_index])
    print("Vt:", V_values[throat_index])
    
    Tw = 1000
    At = A_values[throat_index]
    Mt = M_values[throat_index]

    hgt = calculate_hg(P0, At, gamma,R,Cp,T0,Tw, At, Mt)
    print("hgt:", hgt)
    Tt= T_values[throat_index]
    twt= 0.003
    Kwt= 100
    #(Cp, mu, k, m_dot, area, D)
    hlt= constant_hl_value(3000, 1036*10**(-6), 0.17, 0.14, 4*10**(-6)*50, 0.002)
    (qt, Tlt, Twgt, Twlt) = SS_HT_analysis_equation_solver(hgt,Tt,twt,Kwt,hlt)
    print("T wall *:",Twgt)
    print("T coolant *:",Tlt)
    print(" hl coolant *:", hl)
    print("Isp=", V_values[-1]/9.81)









def Heat_Transfer_Analysis_1():
    (Pe, P0, F, CR, char_L, T0, rho0, Cp, gamma, R, alpha) = assign_constants(c)

    (A0, At, Ae, Lc, m_dot) = chamber_nozzle_sizing(Pe, P0, F, CR, char_L, T0, gamma, R)

    (x_values, r_values, throat_index, chamber_end_index) = conical_radius_function(A0, At, Ae, Lc, alpha)
    A_values = conical_area_function(A0, At, Ae, Lc, alpha)[1]

    M_values = mach_no_as_function_of_distance(A0, At, Ae, Lc, alpha, gamma)[1]
    P_values = pressure_as_function_of_distance(A0, At, Ae, Lc, alpha, gamma,P0)[1]
    rho_values = density_as_function_of_distance(A0, At, Ae, Lc, alpha, gamma,rho0)[1]
    T_values = temperature_as_function_of_distance(A0, At, Ae, Lc, alpha, gamma,T0)[1]
    V_values = velocity_as_function_of_distance(A0, At, Ae, Lc, alpha, gamma,R,T0)[1]
    hg_values = hg_as_function_of_distance(P0,A0, At, Ae, Lc, alpha, gamma,R,Cp,T0,1000)[1]

    plots1 = (
        
        (x_values, T_values, 'X Values', 'Temperature (K)'),
        (x_values, hg_values, 'X Values', 'hg')
        #(x_values2, q_values, 'X Values', 'q'),
        #(x_values2, Tl_values, 'X Values', 'Tl'),
        #(x_values2, Twg_values, 'X Values', 'Twg'),
        #(x_values2, Twl_values, 'X Values', 'Twl')

    )
    plot_multiple_sets(plots1)

    Nc = 27        # No. channels
    hc = 0.002     # channel height (2mm)

    Wc_values = channel_width_as_function_of_distance(Nc, r_values)
    Ac_values = channel_area_as_function_of_distance(hc, Wc_values)

    tw = 0.001   # inner wall thickness (1mm)
    Kw = 25     # wall conductivity of steel

    Tl_in = 290   # room temp 20 degrees c


    (Twg_values, Twl_values, Tl_values, q_dot_values) = Sequential_Heat_Transfer_Analysis(x_values, r_values, M_values, T_values, hg_values, gamma, Wc_values, Ac_values, tw, Kw, Tl_in)

    Coolant_Re_values = Re_no_as_function_of_distance(Wc_values, Ac_values,Tl_values)
    plots = (
        
        (x_values, T_values, 'X Values', 'gas Temperature (K)'),
        (x_values, Twg_values, 'X Values', 'gas side wall remp (K)'),

        (x_values, Twl_values, 'X Values', 'coolant side wall temp'),
        (x_values, Tl_values, 'X Values', 'coolant temp'),
        (x_values, q_dot_values, 'X Values', 'heat flux'),
        (x_values, Coolant_Re_values, 'X Values', 'Coolant Re'),
        #(x_values2, Twl_values, 'X Values', 'Twl')

    )

    plot_multiple_sets(plots)


def Heat_Transfer_Analysis_2():
    (Pe, P0, F, CR, char_L, T0, rho0, Cp, gamma, R, alpha) = assign_constants(c)

    (A0, At, Ae, Lc,m_dot) = chamber_nozzle_sizing(Pe, P0, F, CR, char_L, T0, gamma, R)

    print("m_dot: ",m_dot)

    Lch = 47.32 *10**(-3)
    Lt = 120 *10**(-3)
    L = 161.26 *10**(-3)
    b = 30*np.pi/180
    R2 = 90.41 *10**(-3)
    R1 = 12.98 *10**(-3)
    Rn = 3.31 *10**(-3)
    Tn = 18*np.pi/180
    Te = 8*np.pi/180

    (x_values, r_values, throat_index, chamber_end_index) = parabolic_radius_function(A0, At, Ae, Lch, Lt, L, b, R2, R1, Rn, Tn, Te)
    A_values = parabolic_area_function(A0, At, Ae, Lch, Lt, L, b, R2, R1, Rn, Tn, Te)[1]

    M_values = mach_no_as_function_of_distance(A0, At, Ae, Lch, Lt, L, b, R2, R1, Rn, Tn, Te, gamma)[1]
    P_values = pressure_as_function_of_distance(A0, At, Ae, Lch, Lt, L, b, R2, R1, Rn, Tn, Te, gamma,P0)[1]
    rho_values = density_as_function_of_distance(A0, At, Ae, Lch, Lt, L, b, R2, R1, Rn, Tn, Te, gamma,rho0)[1]
    T_values = temperature_as_function_of_distance(A0, At, Ae, Lch, Lt, L, b, R2, R1, Rn, Tn, Te, gamma,T0)[1]
    V_values = velocity_as_function_of_distance(A0, At, Ae, Lch, Lt, L, b, R2, R1, Rn, Tn, Te, gamma,R,T0)[1]
    hg_values = hg_as_function_of_distance(P0,A0, At, Ae, Lch, Lt, L, b, R2, R1, Rn, Tn, Te, gamma,R,Cp,T0,1000)[1]

    plots1 = (
        
        (x_values, T_values, 'X Values', 'Temperature (K)'),
        (x_values, hg_values, 'X Values', 'hg')
        #(x_values2, q_values, 'X Values', 'q'),
        #(x_values2, Tl_values, 'X Values', 'Tl'),
        #(x_values2, Twg_values, 'X Values', 'Twg'),
        #(x_values2, Twl_values, 'X Values', 'Twl')

    )
    plot_multiple_sets(plots1)

    Nc = 27        # No. channels
    hc = 0.002     # channel height (2mm)

    Wc_values = channel_width_as_function_of_distance(Nc, r_values)
    Ac_values = channel_area_as_function_of_distance(hc, Wc_values)

    tw = 0.001   # inner wall thickness (1mm)
    Kw = 20     # wall conductivity of steel

    Tl_in = 290   # room temp 20 degrees c
    Pl_in = 50*101325
    m_dot = 0.16/27


    (Twg_values, Twl_values, Tl_values, Pl_values, q_dot_values , q_dot_g_conv_values, q_dot_g_rad_values) = Sequential_Heat_Transfer_Analysis_Extended(x_values, r_values, M_values, T_values, hg_values, gamma, Wc_values, Ac_values, tw, Kw, Tl_in, Pl_in, m_dot)

    Coolant_Re_values = Re_no_as_function_of_distance(Wc_values, Ac_values, Tl_values, Pl_values, m_dot)

    ( rho_l_values, Ul_values) = get_rho_and_velocity_values(m_dot, Ac_values,Tl_values, Pl_values)

    write_arrays_to_file("test.txt", x_values*10**3, r_values*10**3, hg_values*10**(-3), q_dot_g_conv_values*10**(-3), q_dot_g_rad_values*10**(-3), q_dot_values*10**(-3),
                     Twg_values, Twg_values, Twl_values, Tl_values, Pl_values*10**(-6), Ul_values, rho_l_values)

    plots = (
        
        (x_values, T_values, 'X Values', 'gas Temperature (K)'),
        (x_values, Twg_values, 'X Values', 'gas side wall remp (K)'),

        (x_values, Twl_values, 'X Values', 'coolant side wall temp'),
        (x_values, Tl_values, 'X Values', 'coolant temp'),
        (x_values, Pl_values, 'X Values', 'coolant P'),
        (x_values, q_dot_g_conv_values, 'X Values', 'conv heat flux'),
        (x_values, q_dot_g_rad_values, 'X Values', 'rad heat flux'),
        (x_values, q_dot_values, 'X Values', 'heat flux'),
        (x_values, Coolant_Re_values, 'X Values', 'Coolant Re'),
        (x_values, Ul_values, 'X Values', 'Coolant vel'),
        (x_values, rho_l_values, 'X Values', 'Coolant rho'),
        #(x_values2, Twl_values, 'X Values', 'Twl')

    )

    plot_multiple_sets(plots)





def geometry_analysis():
    (Pe, P0, F, CR, char_L, T0, rho0, Cp, gamma, R, alpha) = assign_constants(c)

    (A0, At, Ae, Lc,m_dot) = chamber_nozzle_sizing(Pe, P0, F, CR, char_L, T0, gamma, R)

    Lch = 47.32 *10**(-3)
    Lt = 120 *10**(-3)
    L = 161.26 *10**(-3)
    b = 30*np.pi/180
    R2 = 90.41 *10**(-3)
    R1 = 12.98 *10**(-3)
    Rn = 3.31 *10**(-3)
    Tn = 18*np.pi/180
    Te = 8*np.pi/180
  
    (x_values, r_values, throat_index, chamber_end_index) = parabolic_radius_function(A0, At, Ae, Lch, Lt, L, b, R2, R1, Rn, Tn, Te)
    A_values = parabolic_area_function(A0, At, Ae, Lch, Lt, L, b, R2, R1, Rn, Tn, Te)[1]

    #M_values = mach_no_as_function_of_distance(A0, At, Ae, Lc, alpha, gamma)[1]
    #P_values = pressure_as_function_of_distance(A0, At, Ae, Lc, alpha, gamma,P0)[1]
    #rho_values = density_as_function_of_distance(A0, At, Ae, Lc, alpha, gamma,rho0)[1]
    #T_values = temperature_as_function_of_distance(A0, At, Ae, Lc, alpha, gamma,T0)[1]
    #V_values = velocity_as_function_of_distance(A0, At, Ae, Lc, alpha, gamma,R,T0)[1]
    #hg_values = hg_as_function_of_distance(P0,A0, At, Ae, Lc, alpha, gamma,R,Cp,T0,1000)[1]


    chamber_nozzle_geometry(x_values, r_values, [(1, "Rc"),(throat_index, "Rt"),(-1, "Re")], x_label='Axial Position (x)', r_label='Radius (r)', title='Rocket Engine Thrust Chamber Geometry')

    print("Lc: ", Lc)
    print("Lt: ", x_values[throat_index])
    print("L: ", x_values[-1])


def Stress_Analysis():
    (Pe, P0, F, CR, char_L, T0, rho0, Cp, gamma, R, alpha) = assign_constants(c)

    (A0, At, Ae, Lc,m_dot) = chamber_nozzle_sizing(Pe, P0, F, CR, char_L, T0, gamma, R)

    Lch = 47.32 *10**(-3)
    Lt = 120 *10**(-3)
    L = 161.26 *10**(-3)
    b = 30*np.pi/180
    R2 = 90.41 *10**(-3)
    R1 = 12.98 *10**(-3)
    Rn = 3.31 *10**(-3)
    Tn = 18*np.pi/180
    Te = 8*np.pi/180

    (x_values, r_values, throat_index, chamber_end_index) = parabolic_radius_function(A0, At, Ae, Lch, Lt, L, b, R2, R1, Rn, Tn, Te)
    A_values = parabolic_area_function(A0, At, Ae, Lch, Lt, L, b, R2, R1, Rn, Tn, Te)[1]

    M_values = mach_no_as_function_of_distance(A0, At, Ae, Lch, Lt, L, b, R2, R1, Rn, Tn, Te, gamma)[1]
    P_values = pressure_as_function_of_distance(A0, At, Ae, Lch, Lt, L, b, R2, R1, Rn, Tn, Te, gamma,P0)[1]
    rho_values = density_as_function_of_distance(A0, At, Ae, Lch, Lt, L, b, R2, R1, Rn, Tn, Te, gamma,rho0)[1]
    T_values = temperature_as_function_of_distance(A0, At, Ae, Lch, Lt, L, b, R2, R1, Rn, Tn, Te, gamma,T0)[1]
    V_values = velocity_as_function_of_distance(A0, At, Ae, Lch, Lt, L, b, R2, R1, Rn, Tn, Te, gamma,R,T0)[1]
    hg_values = hg_as_function_of_distance(P0,A0, At, Ae, Lch, Lt, L, b, R2, R1, Rn, Tn, Te, gamma,R,Cp,T0,1000)[1]

    plots1 = (
        
        (x_values, T_values, 'X Values', 'Temperature (K)'),
        (x_values, hg_values, 'X Values', 'hg')
        #(x_values2, q_values, 'X Values', 'q'),
        #(x_values2, Tl_values, 'X Values', 'Tl'),
        #(x_values2, Twg_values, 'X Values', 'Twg'),
        #(x_values2, Twl_values, 'X Values', 'Twl')

    )
    plot_multiple_sets(plots1)

    Nc = 27        # No. channels
    hc = 0.002     # channel height (2mm)

    Wc_values = channel_width_as_function_of_distance(Nc, r_values)
    Ac_values = channel_area_as_function_of_distance(hc, Wc_values)

    tw = 0.001   # inner wall thickness (1mm)
    Kw = 20     # wall conductivity of steel
    E = 190*10**(9)
    v = 0.29
    t_outer = 0.003
    Alpha = 16*10**(-6)


    Tl_in = 290   # room temp 20 degrees c
    Pl_in = 50*101325
    m_dot = 0.16/27



    (Twg_values, Twl_values, Tl_values, Pl_values, q_dot_values , q_dot_g_conv_values, q_dot_g_rad_values) = Sequential_Heat_Transfer_Analysis_Extended(x_values, r_values, M_values, T_values, hg_values, gamma, Wc_values, Ac_values, tw, Kw, Tl_in, Pl_in, m_dot)

    (sigma_hoop_values, sigma_r_values, sigma_axial_max_values, sigma_VM_values) = stress_analysis_calculations(Twg_values, Twl_values, q_dot_values, Pl_values, P_values, Wc_values, tw, E,v, Alpha,Kw, Pl_in, Ac_values[-1], P_values[0], A0, At,r_values[throat_index],hc, t_outer)

    data_sets = (
        (sigma_hoop_values*10**(-6), 'sigma_hoop_values'),
        (sigma_r_values*10**(-6), 'sigma_r_values'),
        (sigma_axial_max_values*10**(-6), 'sigma_axial_max_values'),
        (sigma_VM_values*10**(-6), 'sigma_VM_values')
    )

    plot_on_same_axes(x_values, data_sets)







def Manufacture_GW_analysis():
    (Pe, P0, F, CR, char_L, T0, rho0, Cp, gamma, R, alpha) = assign_constants(c)

    (A0, At, Ae, Lc,m_dot) = chamber_nozzle_sizing(Pe, P0, F, CR, char_L, T0, gamma, R)

    print("m_dot: ",m_dot)

    Lch = 47.32 *10**(-3)
    Lt = 120 *10**(-3)
    L = 161.26 *10**(-3)
    b = 30*np.pi/180
    R2 = 90.41 *10**(-3)
    R1 = 12.98 *10**(-3)
    Rn = 3.31 *10**(-3)
    Tn = 18*np.pi/180
    Te = 8*np.pi/180

    (x_values, r_values, throat_index, chamber_end_index) = parabolic_radius_function(A0, At, Ae, Lch, Lt, L, b, R2, R1, Rn, Tn, Te)
    A_values = parabolic_area_function(A0, At, Ae, Lch, Lt, L, b, R2, R1, Rn, Tn, Te)[1]

    M_values = mach_no_as_function_of_distance(A0, At, Ae, Lch, Lt, L, b, R2, R1, Rn, Tn, Te, gamma)[1]
    P_values = pressure_as_function_of_distance(A0, At, Ae, Lch, Lt, L, b, R2, R1, Rn, Tn, Te, gamma,P0)[1]
    rho_values = density_as_function_of_distance(A0, At, Ae, Lch, Lt, L, b, R2, R1, Rn, Tn, Te, gamma,rho0)[1]
    T_values = temperature_as_function_of_distance(A0, At, Ae, Lch, Lt, L, b, R2, R1, Rn, Tn, Te, gamma,T0)[1]
    V_values = velocity_as_function_of_distance(A0, At, Ae, Lch, Lt, L, b, R2, R1, Rn, Tn, Te, gamma,R,T0)[1]
    hg_values = hg_as_function_of_distance(P0,A0, At, Ae, Lch, Lt, L, b, R2, R1, Rn, Tn, Te, gamma,R,Cp,T0,1000)[1]

    Nc = 27        # No. channels
    hc = 0.002     # channel height (2mm)

    Wc_values = channel_width_as_function_of_distance(Nc, r_values)
    Ac_values = channel_area_as_function_of_distance(hc, Wc_values)

    tw = 0.00075   # inner wall thickness (1mm)
    Kw = 20     # wall conductivity of steel

    Tl_in = 290   # room temp 20 degrees c
    Pl_in = 50*101325
    m_dot = 0.16/27


    (Twg_values, Twl_values, Tl_values, Pl_values, q_dot_values , q_dot_g_conv_values, q_dot_g_rad_values) = Sequential_Heat_Transfer_Analysis_Extended(x_values, r_values, M_values, T_values, hg_values, gamma, Wc_values, Ac_values, tw, Kw, Tl_in, Pl_in, m_dot)

    Coolant_Re_values = Re_no_as_function_of_distance(Wc_values, Ac_values, Tl_values, Pl_values, m_dot)

    ( rho_l_values, Ul_values) = get_rho_and_velocity_values(m_dot, Ac_values,Tl_values, Pl_values)

    T_max_working = 1323

    #write_arrays_to_file("tw_0.1.txt", x_values*10**3, r_values*10**3, hg_values*10**(-3), q_dot_g_conv_values*10**(-3), q_dot_g_rad_values*10**(-3), q_dot_values*10**(-3),
    #                 Twg_values, Twg_values, Twl_values, Tl_values, Pl_values*10**(-6), Ul_values, rho_l_values)

    ### plots ###

    chamber_nozzle_geometry(x_values*10**3, r_values*10**3, [(1, "Rc"),(throat_index, "Rt"),(-1, "Re")], x_label='Axial Position (x)', r_label='Radius (r)', title='Rocket Engine Thrust Chamber Geometry')
    print("Twg_max :", Twg_values[throat_index])
    #plot_temperature_data(T_values, Twg_values, Twl_values, Tl_values, x_values*10**3, 
    #                      T_max_working, labels=["Gas temp, Tg.", "Gas side wall temp, Twg", "Coolant side wall temp, Twl", "Coolant temp, Tl"], x_label="Axial Distance (mm)", y_label="Temperature (K)", title="Temperature Profiles Along Nozzle", file_name="temperature_plot_test.png")

def process_RPA_output():

    input_path_no_film = "results/RPA-output-no-film.txt"
    output_path_no_film = "results/RPA-results-no-film.txt"
    input_path_with_film = "results/RPA-output-with-film.txt"
    output_path_with_film = "results/RPA-results-with-film.txt"

    ########### only run one time:
    """
    filter_and_extract_numbers(input_path_no_film, output_path_no_film)
    filter_and_extract_numbers(input_path_with_film, output_path_with_film)
    """
    ###########

    x_vals_no, r_vals_no, q_vals_no, Twg_vals_no, Twc_vals_no, Tc_vals_no, Pc_vals_no = split_data_columns(output_path_no_film)

    data_sets_no = (
        (Twg_vals_no, 'Twg'),
        (Twc_vals_no, 'Twc'),
        (Tc_vals_no, 'Tc')
    )

    q_data_sets_no = (
        (q_vals_no, 'q'),
        #(Twc_vals_no, 'Twc'),
    )

    x_vals_with, r_vals_with, q_vals_with, Twg_vals_with, Twc_vals_with, Tc_vals_with, Pc_vals_with = split_data_columns(output_path_with_film)
    q_data_sets_with = (
        (q_vals_with, 'q'),
        #(Twc_vals_no, 'Twc'),
    )

    Ta = 293.15   ######## input #########
    Pa = 101325 

    k_w = 15
    t_w = 0.003

    Ts_values = calc_external_wall_temp(q_vals_with,
                        Twc_vals_with,
                        Tc_vals_with,
                        r_vals_with,
                        Ta, Pa,
                        k_w,    # wall thermal conductivity
                        t_w)

    data_sets_with = (
        (Twg_vals_with, 'Twg'),
        (Twc_vals_with, 'Twc'),
        (Tc_vals_with, 'Tc'),
        (Ts_values, 'Ts')
    )

    #plot_on_same_axes(x_vals_no, data_sets_no)
    #plot_on_same_axes(x_vals_no, q_data_sets_no)
    plot_on_same_axes(x_vals_with, data_sets_with)
    plot_on_same_axes(x_vals_with, q_data_sets_with)

    x1 = 0.02 *10**3
    x2 = 0.12*10**3

    #i1,i2 = plot_Ts_vs_x(x_vals_with, Ts_values, x1, x2)
    i1,i2 = plot_Ts_vs_x2(x_vals_with, Ts_values, x1, x2, "Ts-plot")

def Final_Analysis():
    (Pe, P0, F, CR, char_L, T0, rho0, Cp, gamma, R, alpha) = assign_constants(c)

    (A0, At, Ae, Lc,m_dot) = chamber_nozzle_sizing(Pe, P0, F, CR, char_L, T0, gamma, R)

    # parabolic nozzle geometries
    Lch = 47.32 *10**(-3)
    Lt = 120 *10**(-3)
    L = 161.26 *10**(-3)
    b = 30*np.pi/180
    R2 = 90.41 *10**(-3)
    R1 = 12.98 *10**(-3)
    Rn = 3.31 *10**(-3)
    Tn = 18*np.pi/180
    Te = 8*np.pi/180

    (x_values, r_values, throat_index, chamber_end_index) = parabolic_radius_function(A0, At, Ae, Lch, Lt, L, b, R2, R1, Rn, Tn, Te)
    A_values = parabolic_area_function(A0, At, Ae, Lch, Lt, L, b, R2, R1, Rn, Tn, Te)[1]

    M_values = mach_no_as_function_of_distance(A0, At, Ae, Lch, Lt, L, b, R2, R1, Rn, Tn, Te, gamma)[1]
    P_values = pressure_as_function_of_distance(A0, At, Ae, Lch, Lt, L, b, R2, R1, Rn, Tn, Te, gamma,P0)[1]
    rho_values = density_as_function_of_distance(A0, At, Ae, Lch, Lt, L, b, R2, R1, Rn, Tn, Te, gamma,rho0)[1]
    T_values = temperature_as_function_of_distance(A0, At, Ae, Lch, Lt, L, b, R2, R1, Rn, Tn, Te, gamma,T0)[1]
    V_values = velocity_as_function_of_distance(A0, At, Ae, Lch, Lt, L, b, R2, R1, Rn, Tn, Te, gamma,R,T0)[1]
    hg_values = hg_as_function_of_distance(P0,A0, At, Ae, Lch, Lt, L, b, R2, R1, Rn, Tn, Te, gamma,R,Cp,T0,1000)[1]

    Nc = 27        # No. channels
    hc = 0.002     # channel height (2mm)

    Wc_values = channel_width_as_function_of_distance(Nc, r_values)
    Ac_values = channel_area_as_function_of_distance(hc, Wc_values)

    tw = 0.0003   # inner wall thickness (0.3mm)
    Kw = 20     # wall conductivity of steel

    Tl_in = 290   # room temp 20 degrees c
    Pl_in = 40*101325
    m_dot_l = 0.16/27

    # (Twg_values, Twl_values, Tl_values, Pl_values, q_dot_values , q_dot_g_conv_values, q_dot_g_rad_values) = Sequential_Heat_Transfer_Analysis_Extended(x_values, r_values, M_values, T_values, hg_values, gamma, Wc_values, Ac_values, tw, Kw, Tl_in, Pl_in, m_dot)

    # Temp_data_sets = (
    #     (Twg_values, 'Twg'),
    #     (Twl_values, 'Twl'),
    #     (Tl_values, 'Tl')
    # )

    # plot_on_same_axes(x_values, Temp_data_sets)

    # data_sets = (
    #     (M_values, 'rho'),
    # )

    # plot_on_same_axes(x_values, data_sets)

    (Twg_values, Twl_values, Tl_values, Pl_values, q_dot_values , q_dot_g_conv_values, q_dot_g_rad_values) = Sequential_Heat_Transfer_Analysis_Extended2(x_values, r_values, M_values, T_values, hg_values, gamma, Wc_values, Ac_values, tw, Kw, Tl_in, Pl_in, m_dot_l, A_values, chamber_end_index, throat_index, R, P0, Cp, m_dot)
    #(Twg_values, Twl_values, Tl_values, Pl_values, q_dot_values , q_dot_g_conv_values, q_dot_g_rad_values) = Sequential_Heat_Transfer_Analysis_Extended3(x_values, r_values, M_values,V_values, T_values, P_values,rho_values, hg_values, gamma, Wc_values, Ac_values, tw, Kw, Tl_in, Pl_in, m_dot_l, A_values, chamber_end_index, throat_index, R, P0, Cp, m_dot)  # this anaysis is wrong

    Temp_data_sets = (
        (Twg_values, 'Twg'),
        (Twl_values, 'Twl'),
        (Tl_values, 'Tl')
    )

    q_data_sets = (
        (q_dot_values, 'q'),

    )

    plot_on_same_axes(x_values, Temp_data_sets)
    plot_on_same_axes(x_values, q_data_sets)

    #print( hg_values)
    #print((A0, At, Ae, Lc,m_dot))

    Ta = 293.15   ######## input #########
    Pa = 101325 

    k_w = 15
    t_w = 0.003

    Ts_values = calc_external_wall_temp(q_dot_values,
                        Twl_values,
                        Tl_values,
                        r_values,
                        Ta, Pa,
                        k_w,    # wall thermal conductivity
                        t_w)
    
    x1 = 0.02 *10**3
    x2 = 0.12 *10**3

    #i1,i2 = plot_Ts_vs_x(x_vals_with, Ts_values, x1, x2)
    i1,i2 = plot_Ts_vs_x2(x_values*10**3, Ts_values, x1, x2, "Ts-plot-model")

    print(M_values)

"""
def Hot_Fire_data():
    # Define inputs
    input_file = "data/Hot fire data 20250524_DMT_Full_Dataset.xlsx"
    output_folder = "processed_data"
    time_windows = [(3500, 5000), (7500,9000 )]  # Example start/end times


    # Run the function once
    split_tests_by_time(input_file, output_folder, time_windows)

    df = pd.read_excel("processed_data/test1_raw.xlsx")
    print("Available columns:", df.columns.tolist())
    time = df['Time(s)']
    p_chamber = df['P_Chamber(bar)']
    t_throat = df['T_Throat_Wall(C)']
    t_chamber = df['T_Chamber_Wall(C)']

    Temp_data_sets = (
        (t_chamber, 'Tc'),

    )

    plot_on_same_axes(time, Temp_data_sets)
    
"""
def find_column(columns, key_substring):
    for col in columns:
        if key_substring.lower() in col.lower():
            return col
    raise ValueError(f"Column with keyword '{key_substring}' not found.")

def Hot_Fire_data():
    # --- Helper function scoped inside ---
    def find_column(columns, key_substring):
        for col in columns:
            if key_substring.lower() in col.lower():
                return col
        raise ValueError(f"Column with keyword '{key_substring}' not found.")

    # Define inputs
    input_file = "data/Hot fire data 20250524_DMT_Full_Dataset.xlsx"
    output_folder = "data"  # Updated to save in the same folder as the input

    time_windows = [(3726., 3755), (7748.9, 7820)]  # Start/end times for each test

    # Run the function once to split raw test data
    #split_tests_by_time(input_file, output_folder, time_windows)

    # Load Test 2 data (after split)
    df = pd.read_excel(os.path.join(output_folder, "test2_raw.xlsx"))
    df.columns = df.columns.str.strip()  # Clean column names

    print("Available columns:", df.columns.tolist())  # For inspection

    # Use robust substring matching to find the needed columns
    time = df[find_column(df.columns, "time")]
    p_chamber = df[find_column(df.columns, "p_chamber")]
    p_fuel_tank = df[find_column(df.columns, "p_fuel_tank")]
    p_ox_tank = df[find_column(df.columns, "p_ox_tank")]
    t_throat = df[find_column(df.columns, "t_throat")]
    t_chamber = df[find_column(df.columns, "t_chamber")]

    # Package data for plotting
    Temp_data_sets = [
        (t_chamber, 'tc'),
        (t_throat, 'tt'),
        #(p_fuel_tank, 'pft'),
        #(p_ox_tank, 'pot'),
        # Optionally add more series here
    ]

    # Call your existing plotting function
    plot_on_same_axes(time, Temp_data_sets)

def Report_Analysis():
    ##############################
    ##### 1 # modal analysis #####
    ##############################
    (Pe, P0, F, CR, char_L, T0, rho0, Cp, gamma, R, alpha) = assign_constants(c)

    (A0, At, Ae, Lc, m_dot) = chamber_nozzle_sizing(Pe, P0, F, CR, char_L, T0, gamma, R)

    # parabolic nozzle geometries
    Lch = 47.32 *10**(-3)
    Lt = 120 *10**(-3)
    L = 161.26 *10**(-3)
    b = 30*np.pi/180
    R2 = 90.41 *10**(-3)
    R1 = 12.98 *10**(-3)
    Rn = 3.31 *10**(-3)
    Tn = 18*np.pi/180
    Te = 8*np.pi/180

    (x_values, r_values, throat_index, chamber_end_index) = parabolic_radius_function(A0, At, Ae, Lch, Lt, L, b, R2, R1, Rn, Tn, Te)
    A_values = parabolic_area_function(A0, At, Ae, Lch, Lt, L, b, R2, R1, Rn, Tn, Te)[1]

    M_values = mach_no_as_function_of_distance(A0, At, Ae, Lch, Lt, L, b, R2, R1, Rn, Tn, Te, gamma)[1]
    P_values = pressure_as_function_of_distance(A0, At, Ae, Lch, Lt, L, b, R2, R1, Rn, Tn, Te, gamma,P0)[1]
    rho_values = density_as_function_of_distance(A0, At, Ae, Lch, Lt, L, b, R2, R1, Rn, Tn, Te, gamma,rho0)[1]
    T_values = temperature_as_function_of_distance(A0, At, Ae, Lch, Lt, L, b, R2, R1, Rn, Tn, Te, gamma,T0)[1]
    V_values = velocity_as_function_of_distance(A0, At, Ae, Lch, Lt, L, b, R2, R1, Rn, Tn, Te, gamma,R,T0)[1]
    hg_values = hg_as_function_of_distance(P0,A0, At, Ae, Lch, Lt, L, b, R2, R1, Rn, Tn, Te, gamma,R,Cp,T0,1000)[1]

    print(M_values)

    Nc = 27        # No. channels
    hc = 0.002     # channel height (2mm)

    Wc_values = channel_width_as_function_of_distance(Nc, r_values)
    Ac_values = channel_area_as_function_of_distance(hc, Wc_values)

    tw = 0.0003   # inner wall thickness (0.3mm)
    Kw = 20     # wall conductivity of steel

    Tl_in = 290   # room temp 20 degrees c
    Pl_in = 40*101325
    m_dot_l = 0.16/27

    (Twg_values, Twl_values, Tl_values, Pl_values, q_dot_values , q_dot_g_conv_values, q_dot_g_rad_values) = Sequential_Heat_Transfer_Analysis_Extended2(x_values, r_values, M_values, T_values, hg_values, gamma, Wc_values, Ac_values, tw, Kw, Tl_in, Pl_in, m_dot_l, A_values, chamber_end_index, throat_index, R, P0, Cp, m_dot)
    #(Twg_values, Twl_values, Tl_values, Pl_values, q_dot_values , q_dot_g_conv_values, q_dot_g_rad_values) = Sequential_Heat_Transfer_Analysis_Extended3(x_values, r_values, M_values,V_values, T_values, P_values,rho_values, hg_values, gamma, Wc_values, Ac_values, tw, Kw, Tl_in, Pl_in, m_dot_l, A_values, chamber_end_index, throat_index, R, P0, Cp, m_dot)  # this anaysis is wrong

    ############################
    ##### 2 # RPA analysis #####
    ############################
    output_path_no_film = "results/RPA-results-no-film.txt"
    output_path_with_film = "results/RPA-results-with-film.txt"
    x_vals_no, r_vals_no, q_vals_no, Twg_vals_no, Twc_vals_no, Tc_vals_no, Pc_vals_no = split_data_columns(output_path_no_film)
    x_vals_with, r_vals_with, q_vals_with, Twg_vals_with, Twc_vals_with, Tc_vals_with, Pc_vals_with = split_data_columns(output_path_with_film)
    scale = 0.12/0.09448
    x_vals_no = np.array(x_vals_no)*scale/1000 # correct offset
    x_vals_with = np.array(x_vals_with)*scale/1000

    Twg_vals_no_int = np.interp(x_values, x_vals_no, Twg_vals_no)
    Twc_vals_no_int = np.interp(x_values, x_vals_no, Twc_vals_no)
    Tc_vals_no_int = np.interp(x_values, x_vals_no, Tc_vals_no)
    Twg_vals_with_int = np.interp(x_values, x_vals_with, Twg_vals_with)
    Twc_vals_with_int = np.interp(x_values, x_vals_with, Twc_vals_with)
    Tc_vals_with_int = np.interp(x_values, x_vals_with, Tc_vals_with)
    #Twg_values_film, Twl_values_film, Tl_values_film = Twg_values - (Twg_vals_no_int - Twg_vals_with_int), Twl_values - (Twc_vals_no_int - Twc_vals_with_int), Tl_values - (Tc_vals_no_int - Tc_vals_with_int)
    # DTwg = Twg_vals_no_int[throat_index] - Twg_vals_with_int[throat_index]
    # DTwl = Twc_vals_no_int[throat_index] - Twc_vals_with_int[throat_index]
    # DTl = Tc_vals_no_int[throat_index] - Tc_vals_with_int[throat_index]
    dt1 = np.max(Twg_vals_no_int) - np.max(Twg_vals_with_int)
    dt2 = np.max(Tc_vals_no_int) - np.max(Tc_vals_with_int)
    Twg_values_film = Twg_values - dt1
    Twl_values_film = Twl_values - dt1
    Tl_values_film  = Tl_values  - dt2

    print(M_values)
    
    ###########################################
    ##### 3 # temp plot for modal and RPA ##### plots DONE and saved
    ###########################################
    # T_max_working = 1323
    # plot_temperature_data(T_values, Twg_values_film, Twl_values_film, Tl_values_film, x_values*10**3, 
    #     T_max_working, labels=["Gas temp, Tg.", "Gas side wall temp, Twg", "Coolant side wall temp, Twl", "Coolant temp, Tl"], x_label="Axial Distance (mm)", y_label="Temperature (K)", title="Modelled Temperature Profiles Along Nozzle", file_name="modal_temperature_plots.png")
    # plot_temperature_data(T_values, Twg_vals_with_int, Twc_vals_with_int, Tc_vals_with_int, x_values*10**3, 
    #     T_max_working, labels=["Gas temp, Tg.", "Gas side wall temp, Twg", "Coolant side wall temp, Twl", "Coolant temp, Tl"], x_label="Axial Distance (mm)", y_label="Temperature (K)", title="RPA Temperature Profiles Along Nozzle", file_name="RPA_temperature_plots.png")

    #########################################
    ##### 4 # Ts plot for modal and RPA #####
    #########################################
    Ta = 293.15   ######## input #########
    Pa = 101325   # atm pressure

    k_w = 15  # at this temp range
    t_w = 0.003  # 3mm

    # q_vals_with_int = np.interp(x_values, x_vals_with, q_vals_with)
    # Ts_values_Modal = calc_external_wall_temp(np.array(q_vals_with_int)*1000,
    #                     Twl_values_film,
    #                     Tl_values_film,
    #                     r_values ,
    #                     Ta, Pa,
    #                     k_w,    # wall thermal conductivity
    #                     t_w)
    
    Ts_values_RPA = calc_external_wall_temp(np.array(q_vals_with)*1000,
                    Twc_vals_with,
                    Tc_vals_with,
                    r_vals_with,
                    Ta, Pa,
                    k_w,    # wall thermal conductivity
                    t_w)
    
    x1 = 0.02 *10**3
    x2 = 0.12 *10**3

    #i1,i2 = plot_Ts_vs_x(x_vals_with, Ts_values, x1, x2)
    #i1,i2 = plot_Ts_vs_x2(x_values*10**3, Ts_values_Modal, x1, x2, "Ts-plot-model") #
    i1,i2 = plot_Ts_vs_x2(x_vals_with*10**3, Ts_values_RPA, x1, x2, "Ts-plot-RPA") # done and saved



    ###################################################
    ##### 4 # Transient Ts plot for RPA #####
    ###################################################
    
    # at chamber
    h_c = get_channel_CHT_coef(q_vals_with[i1]*1000, Twc_vals_with[i1], Tc_vals_with[i1])/1.23 #temp correction
    A_c = 1
    V_c = t_w
    rho_steel = 7900
    c_p_steel = 500
    T_i = 293
    Tc_c = Tc_vals_with[i1]

    dTdt_c = h_c*A_c*(Tc_c-T_i)/(rho_steel*V_c*c_p_steel)

    times_c, temps_c = lumped_capacitance(h_c, A_c, rho_steel, c_p_steel, V_c, T_i, Tc_c, t_step=0.1) ## need a better model
    #times_c, temps_c = transient_conduction_1D(k_w, rho_steel, c_p_steel, t_w, T_i, Tc_c, t_step=0.1)

    # at throat
    h_t = get_channel_CHT_coef(q_vals_with[i2]*1000, Twc_vals_with[i2], Tc_vals_with[i2])/1.23 #temp correction
    A_t = A_c
    V_t = V_c
    Tc_t = Tc_vals_with[i2]

    dTdt_t = h_t*A_t*(Tc_t-T_i)/(rho_steel*V_t*c_p_steel)
    
    times_t, temps_t = lumped_capacitance(h_t, A_t, rho_steel, c_p_steel, V_t, T_i, Tc_t, t_step=0.1)
    #times_t, temps_t = transient_conduction_1D(k_w, rho_steel, c_p_steel, t_w, T_i, Tc_t, t_step=0.1)

    # Calculate how many nans to add
    padding = len(temps_c) - len(temps_t)

    # Append np.nan values
    if padding > 0:
        temps_t = np.append(temps_t, [np.nan] * padding)
    
    # Temp_data_sets = [
    #     (temps_c, 'tc'),
    #     (temps_t, 'tt'),
    #     #(p_fuel_tank, 'pft'),
    #     #(p_ox_tank, 'pot'),
    #     # Optionally add more series here
    # ]

    # # Call your existing plotting function
    # plot_on_same_axes(times_c, Temp_data_sets)

    #x = [1, 2, 3, 4, 5]
    #dataset1 = ("Dataset 1", [2, 4, 6, 8, 10])
    #dataset2 = ("Dataset 2", [1, 3, 5, 7, 9])
    datasets=[
        ("Ts at chamber", temps_c), 
        ("Ts at throat", temps_t)
    ]

    plot_and_save_multiple_datasets(
        times_c,
        datasets,
        x_label="Time (s)",
        y_label="Surface Temp, Ts (K)",
        title="Transient Heat Modal",
        filename="full_transient_plot.png"
    )

    ###################################################
    ##### 5 # Hot fire 10 second plots #####
    ###################################################

    # measured data
    # Load Test 1 data (after split)
    df = pd.read_excel(os.path.join('data', "test1_raw.xlsx"))
    df.columns = df.columns.str.strip()  # Clean column names

    # Use robust substring matching to find the needed columns
    time1 = df[find_column(df.columns, "time")]
    #p_chamber = df[find_column(df.columns, "p_chamber")]
    #p_fuel_tank = df[find_column(df.columns, "p_fuel_tank")]
    #p_ox_tank = df[find_column(df.columns, "p_ox_tank")]
    t_throat1 = df[find_column(df.columns, "t_throat")]
    t_chamber1 = df[find_column(df.columns, "t_chamber")]

    # start_index1 = find_closest_index(time1, 3726.25+3)
    # end_index1 = find_closest_index(time1, 3733.75+3)

    start_index1 = find_closest_index(time1, 3726.25)
    end_index1 = find_closest_index(time1, 3755)

    # modal
    #Ts_modal_index = find_closest_index(times_c, 3733.75-3726.25)
    Ts_modal_index = find_closest_index(times_c, time1.iloc[-1]-3726.25)

    dataset1a = {
        'label': 'Measured at Chamber',
        'x_values': time1[start_index1:end_index1]- time1[start_index1],
        'y_values': t_chamber1[start_index1:end_index1]+273.15,
        'color': 'blue',
        'linestyle': '-'
    }
    dataset1b = {
        'label': 'Predicted at Chamber',
        'x_values': times_c[:Ts_modal_index],
        'y_values': temps_c[:Ts_modal_index],
        'color': 'blue',
        'linestyle': '--'
    }

    dataset2a = {
        'label': 'Measured at Throat',
        'x_values': time1[start_index1:end_index1]- time1[start_index1],
        'y_values': t_throat1[start_index1:end_index1]+273.15,
        'color': 'red',
        'linestyle': '-'
    }
    dataset2b = {
        'label': 'Predicted at Throat',
        'x_values': times_c[:Ts_modal_index],
        'y_values': temps_t[:Ts_modal_index],
        'color': 'red',
        'linestyle': '--'
    }

    events = [
        {
            'time': 7.5,
            'color': 'gray',
            'linestyle': '--',
            'label': 'Combustion Ends'
        }
    ]

    plot_multiple_datasets_varied_x_2(
        datasets=[dataset1a, dataset1b, dataset2a, dataset2b],
        x_label='Time from Ignition (s)',
        y_label='Surface Temperature (K)',
        title='Transient Surface Temperature (Hot Fire 1)',
        filename='10s_plot_test1.png',
        events = events
    )

    # Load Test 2 data (after split)
    df = pd.read_excel(os.path.join('data', "test2_raw.xlsx"))
    df.columns = df.columns.str.strip()  # Clean column names

    # Use robust substring matching to find the needed columns
    time2 = df[find_column(df.columns, "time")]
    #p_chamber = df[find_column(df.columns, "p_chamber")]
    #p_fuel_tank = df[find_column(df.columns, "p_fuel_tank")]
    #p_ox_tank = df[find_column(df.columns, "p_ox_tank")]
    t_throat2 = df[find_column(df.columns, "t_throat")]
    t_chamber2 = df[find_column(df.columns, "t_chamber")]

    # start_index2 = find_closest_index(time2, 7749 + 3)
    # end_index2 = find_closest_index(time2, 7758.5 + 3)

    start_index2 = find_closest_index(time2, 7749)
    end_index2 = find_closest_index(time2, 7820)

    # modal
    #Ts_modal_index = find_closest_index(times_c, 7758.5-7749)
    Ts_modal_index = find_closest_index(times_c, 70)

    dataset1a = {
        'label': 'Measured at Chamber',
        'x_values': time2[start_index2:end_index2]- time2[start_index2],
        'y_values': t_chamber2[start_index2:end_index2]+273.15,
        'color': 'blue',
        'linestyle': '-'
    }
    dataset1b = {
        'label': 'Predicted at Chamber',
        'x_values': times_c[:Ts_modal_index],
        'y_values': temps_c[:Ts_modal_index],
        'color': 'blue',
        'linestyle': '--'
    }

    dataset2a = {
        'label': 'Measured at Throat',
        'x_values': time2[start_index2:end_index2]- time2[start_index2],
        'y_values': t_throat2[start_index2:end_index2]+273.15,
        'color': 'red',
        'linestyle': '-'
    }
    dataset2b = {
        'label': 'Predicted at Throat',
        'x_values': times_c[:Ts_modal_index],
        'y_values': temps_t[:Ts_modal_index],
        'color': 'red',
        'linestyle': '--'
    }


    events = [
        {
            'time': 9.5,
            'color': 'gray',
            'linestyle': '--',
            'label': 'Combustion Ends'
        }
    ]

    plot_multiple_datasets_varied_x_2(
        datasets=[dataset1a, dataset1b, dataset2a, dataset2b],
        x_label='Time from Ignition (s)',
        y_label='Surface Temperature (K)',
        title='Transient Surface Temperature (Hot Fire 2)',
        filename='10s_plot_test2.png',
        events=events
    )

    # print('### 1 ###: throat temp %: ', (np.max(temps_t) - np.max(t_throat1[start_index1:end_index1]+273.15))*100/np.max(temps_t))
    # print('### 1 ###: chamber temp %: ', (abs(np.max(temps_c) - np.max(t_chamber1[start_index1:end_index1]+273.15)))*100/np.max(temps_c))
    # print('### 2 ###: throat temp %: ', (abs(np.max(temps_t) - np.max(t_throat2)+273.15))*100/np.max(temps_t))
    # print('### 2 ###: chamber temp %: ', (abs(np.max(temps_c) - np.max(t_chamber2)+273.15))*100/np.max(temps_c))

    print('grad at c:',dTdt_c)
    print('grad at t:',dTdt_t)


def Presentation_Plots():
    ##############################
    ##### 1 # modal analysis #####
    ##############################
    (Pe, P0, F, CR, char_L, T0, rho0, Cp, gamma, R, alpha) = assign_constants(c)

    (A0, At, Ae, Lc, m_dot) = chamber_nozzle_sizing(Pe, P0, F, CR, char_L, T0, gamma, R)

    # parabolic nozzle geometries
    Lch = 47.32 *10**(-3)
    Lt = 120 *10**(-3)
    L = 161.26 *10**(-3)
    b = 30*np.pi/180
    R2 = 90.41 *10**(-3)
    R1 = 12.98 *10**(-3)
    Rn = 3.31 *10**(-3)
    Tn = 18*np.pi/180
    Te = 8*np.pi/180

    (x_values, r_values, throat_index, chamber_end_index) = parabolic_radius_function(A0, At, Ae, Lch, Lt, L, b, R2, R1, Rn, Tn, Te)
    A_values = parabolic_area_function(A0, At, Ae, Lch, Lt, L, b, R2, R1, Rn, Tn, Te)[1]

    M_values = mach_no_as_function_of_distance(A0, At, Ae, Lch, Lt, L, b, R2, R1, Rn, Tn, Te, gamma)[1]
    P_values = pressure_as_function_of_distance(A0, At, Ae, Lch, Lt, L, b, R2, R1, Rn, Tn, Te, gamma,P0)[1]
    rho_values = density_as_function_of_distance(A0, At, Ae, Lch, Lt, L, b, R2, R1, Rn, Tn, Te, gamma,rho0)[1]
    T_values = temperature_as_function_of_distance(A0, At, Ae, Lch, Lt, L, b, R2, R1, Rn, Tn, Te, gamma,T0)[1]
    V_values = velocity_as_function_of_distance(A0, At, Ae, Lch, Lt, L, b, R2, R1, Rn, Tn, Te, gamma,R,T0)[1]
    hg_values = hg_as_function_of_distance(P0,A0, At, Ae, Lch, Lt, L, b, R2, R1, Rn, Tn, Te, gamma,R,Cp,T0,1000)[1]

    Nc = 27        # No. channels
    hc = 0.002     # channel height (2mm)

    Wc_values = channel_width_as_function_of_distance(Nc, r_values)
    Ac_values = channel_area_as_function_of_distance(hc, Wc_values)

    tw = 0.0003   # inner wall thickness (0.3mm)
    Kw = 20     # wall conductivity of steel

    Tl_in = 290   # room temp 20 degrees c
    Pl_in = 40*101325
    m_dot_l = 0.16/27

    (Twg_values, Twl_values, Tl_values, Pl_values, q_dot_values , q_dot_g_conv_values, q_dot_g_rad_values) = Sequential_Heat_Transfer_Analysis_Extended2(x_values, r_values, M_values, T_values, hg_values, gamma, Wc_values, Ac_values, tw, Kw, Tl_in, Pl_in, m_dot_l, A_values, chamber_end_index, throat_index, R, P0, Cp, m_dot)
    #(Twg_values, Twl_values, Tl_values, Pl_values, q_dot_values , q_dot_g_conv_values, q_dot_g_rad_values) = Sequential_Heat_Transfer_Analysis_Extended3(x_values, r_values, M_values,V_values, T_values, P_values,rho_values, hg_values, gamma, Wc_values, Ac_values, tw, Kw, Tl_in, Pl_in, m_dot_l, A_values, chamber_end_index, throat_index, R, P0, Cp, m_dot)  # this anaysis is wrong

    #chamber_nozzle_geometry(x_values*10**3, r_values*10**3, [(1, "Rc"),(throat_index, "Rt"),(-1, "Re")], x_label='Axial Position (x)', r_label='Radius (r)', title='Rocket Engine Thrust Chamber Geometry')
    chamber_nozzle_geometry_2(x_values*10**3, r_values*10**3, [(1, "Rc"),(throat_index, "Rt"),(-1, "Re")], x_label='Axial Position (mm)', r_label='Radius (mm)', title='Rocket Engine Thrust Chamber Geometry', output_filename='chamber_geometry.png')