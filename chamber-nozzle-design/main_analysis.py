# Main script to run the analysis
from analysis import design_GW_analysis, geometry_analysis, Heat_Transfer_Analysis_1, Heat_Transfer_Analysis_2, Stress_Analysis, Manufacture_GW_analysis, process_RPA_output, Final_Analysis, Hot_Fire_data, Report_Analysis, Presentation_Plots

from data_loader import assign_constants
from nozzle_functions import chamber_nozzle_sizing, conical_radius_function, conical_area_function, mach_no_as_function_of_distance, pressure_as_function_of_distance,temperature_as_function_of_distance,density_as_function_of_distance,velocity_as_function_of_distance
from heat_transfer_functions import hg_as_function_of_distance, constant_hl_value, calculate_hg
from plotting import plot_chamber_nozzle_geometry,plot_multiple_sets, chamber_nozzle_geometry
import constants as c  # Import constants from the constants file

#import matplotlib.pyplot as plt


def main():

    #design_GW_analysis()
    #geometry_analysis()
    #Heat_Transfer_Analysis_1()
    #Heat_Transfer_Analysis_2()
    #Stress_Analysis()
    #Manufacture_GW_analysis()
    #process_RPA_output()
    #Final_Analysis()
    #Hot_Fire_data()
    #Report_Analysis()
    Presentation_Plots()

    



if __name__ == "__main__":
    main()
