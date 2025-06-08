#plotting
import matplotlib.pyplot as plt
import os
import numpy as np
import math

def plot_chamber_nozzle_geometry(x_values,r_values):
    plt.plot(x_values, r_values, 'ko', markersize=2)
    plt.plot(x_values, -r_values, 'ko', markersize=2)

    plt.grid(True)
    plt.axis('equal')
    plt.show()


def chamber_nozzle_geometry(x_values, r_values, important_indices=None, x_label='Axial Position (x)', r_label='Radius (r)', title='Rocket Engine Thrust Chamber Geometry'):
    """
    Plots the chamber geometry and nozzle shape.

    Parameters:
        x_values (list or array): x-coordinates of the geometry.
        r_values (list or array): Corresponding radii values.
        important_indices (list of tuples, optional): Indices of important radii with format [(index, label), ...].
        x_label (str): Label for the x-axis.
        r_label (str): Label for the r-axis.
        title (str): Title of the plot.
    """
    # Plot the chamber geometry as smaller points
    plt.plot(x_values, r_values, 'ko', markersize=1, label="Positive r")
    plt.plot(x_values, -r_values, 'ko', markersize=1, label="Negative r")
    
    # Mark important radii if provided
    if important_indices:
        for index, label in important_indices:
            # Normalize negative indices
            if index < 0:
                index = len(x_values) + index
            
            # Check for valid indices
            if 0 <= index < len(x_values):
                x = x_values[index]
                r = r_values[index]
                # Plot vertical line only up to the positive radius
                plt.vlines(x=x, ymin=0, ymax=r, color='r', linestyle='--', linewidth=1)
                # Annotate the radius value
                plt.text(x+0.005, -5, f'{label} ={r:.4f}', color='r', fontsize=8, ha='center')
            else:
                print(f"Warning: Index {index} is out of range.")
    
    # Labels and title
    plt.xlabel(x_label)
    plt.ylabel(r_label)
    plt.title(title)

    # Additional plot settings
    plt.grid(True)
    plt.axis('equal')
    plt.legend()
    plt.show()

import os
import matplotlib.pyplot as plt

def chamber_nozzle_geometry_2(x_values, r_values, important_indices=None, x_label='Axial Position (x)', r_label='Radius (r)', title='Rocket Engine Thrust Chamber Geometry', output_filename=None):
    """
    Plots the chamber geometry and nozzle shape.

    Parameters:
        x_values (list or array): x-coordinates of the geometry.
        r_values (list or array): Corresponding radii values.
        important_indices (list of tuples, optional): Indices of important radii with format [(index, label), ...].
        x_label (str): Label for the x-axis.
        r_label (str): Label for the r-axis.
        title (str): Title of the plot.
        output_filename (str, optional): Filename to save the plot in the 'results' folder. If None, the plot will just be displayed.
    """
    plt.figure(figsize=(10, 6))
    # Plot the chamber geometry as smaller points
    plt.plot(x_values, r_values, 'ko', markersize=2, label="Positive r")
    plt.plot(x_values, -r_values, 'ko', markersize=2, label="Negative r")
    
    # Mark important radii if provided
    if important_indices:
        for index, label in important_indices:
            # Normalize negative indices
            if index < 0:
                index = len(x_values) + index
            
            # Check for valid indices
            if 0 <= index < len(x_values):
                x = x_values[index]
                r = r_values[index]
                # Plot vertical line only up to the positive radius
                plt.vlines(x=x, ymin=0, ymax=r, color='r', linestyle='--', linewidth=1)
                # Annotate the radius value
                #plt.text(x+0.005, -5, f'{label} = {r:.4f}', color='r', fontsize=8, ha='center')
                plt.text(x+0.005, -5, f'{label} = {r:.1f}', color='r', fontsize=8, ha='center') # 1dp

            else:
                print(f"Warning: Index {index} is out of range.")
    
    # Labels and title
    plt.xlabel(x_label)
    plt.ylabel(r_label)
    plt.title(title)

    # Additional plot settings
    plt.grid(True)
    plt.axis('equal')
    #plt.legend()
    plt.tight_layout()  # Ensures everything fits nicely

    # Save or show the plot
    if output_filename:
        # Ensure the 'results' directory exists
        os.makedirs('results', exist_ok=True)
        save_path = os.path.join('results', output_filename)
        #plt.savefig(save_path)
        plt.savefig(save_path, dpi=600)

        print(f"Plot saved to {save_path}")
        plt.close()  # Close the plot to avoid overlapping when running multiple times
    else:
        plt.show()




def plot_multiple_sets(data_sets):
    """
    This function takes in a tuple of data sets where each data set is a tuple of:
    (x_values, y_values, x_label, y_label).
    example:

    data_sets = (
        (x_values1, M_values, 'X Values', 'Mach Number'),
        (x_values1, P_values, 'X Values', 'Pressure (P)'),
        (x_values1, rho_values, 'X Values', 'Density (ρ)')
    )

    It generates a plot for each set of values using different colors.
    
    Parameters:
    data_sets: tuple of tuples, each containing (x_values, y_values, x_label, y_label)
    """
    # List of marker colors for different plots
    colors = [
        'r',  # Red
        'b',  # Blue
        'g',  # Green
        'm',  # Magenta
        'c',  # Cyan
        'y',  # Yellow
        'k',  # Black
        'orange',  # Orange
        'purple',  # Purple
        'brown',  # Brown
        'pink',  # Pink
        'lime',  # Lime
        'teal',  # Teal
        'navy',  # Navy
        'violet',  # Violet
        'darkred',  # Dark Red
        'gold',  # Gold
        'darkgreen',  # Dark Green
        'gray',  # Gray
        'olive',  # Olive
        'lightblue',  # Light Blue
        'coral',  # Coral
        'turquoise',  # Turquoise
        'salmon',  # Salmon
    ]

    # Loop through each data set and generate plots
    for i, data_set in enumerate(data_sets):
        # Unpack the data set tuple
        x_values, y_values, x_label, y_label = data_set
        
        # Choose a color for the plot (loop over the colors list if there are more sets than colors)
        color = colors[i % len(colors)]
        
        # Create a new figure for each plot
        plt.figure()
        plt.plot(x_values, y_values, marker='o', color=color, markersize=2, label=f'{y_label} vs {x_label}')

        
        # Add labels and grid
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.grid(True)
        plt.legend()
        plt.title(f'{y_label} vs {x_label}')
        
        # Show the plot
        plt.show()

import matplotlib.pyplot as plt

def plot_on_same_axes(x_values, data_sets):
    """
    This function takes in a single set of x values and a tuple of data sets where each data set is a tuple of:
    (y_values, label).
    
    All data sets are plotted on the same set of axes with unique colors and markers.

    Example:

    x_values = [1, 2, 3, 4]
    data_sets = (
        (M_values, 'Mach Number'),
        (P_values, 'Pressure (P)'),
        (rho_values, 'Density (ρ)')
    )

    Parameters:
    x_values: list or array of x-axis values
    data_sets: tuple of tuples, each containing (y_values, label)
    """
    # List of marker colors for different plots
    colors = [
        'r',  # Red
        'b',  # Blue
        'g',  # Green
        'm',  # Magenta
        'c',  # Cyan
        'y',  # Yellow
        'k',  # Black
        'orange',  # Orange
        'purple',  # Purple
        'brown',  # Brown
        'pink',  # Pink
        'lime',  # Lime
        'teal',  # Teal
        'navy',  # Navy
        'violet',  # Violet
        'darkred',  # Dark Red
        'gold',  # Gold
        'darkgreen',  # Dark Green
        'gray',  # Gray
        'olive',  # Olive
        'lightblue',  # Light Blue
        'coral',  # Coral
        'turquoise',  # Turquoise
        'salmon',  # Salmon
    ]

    # Create a new figure for the combined plot
    plt.figure()

    # Loop through each data set and plot on the same axes
    for i, data_set in enumerate(data_sets):
        # Unpack the data set tuple
        y_values, label = data_set

        # Choose a color for the plot (loop over the colors list if there are more sets than colors)
        color = colors[i % len(colors)]

        # Plot the data set on the same axes
        plt.plot(x_values, y_values, marker='o', color=color, markersize=4, label=label)

    # Add labels, legend, grid, and title
    plt.xlabel('X Values')
    plt.ylabel('Y Values')
    plt.grid(True)
    plt.legend()
    plt.title('Combined Plot of All Data Sets')

    # Show the combined plot
    plt.show()

# Example usage
# x_values = [1, 2, 3, 4]
# data_sets = (
#     (M_values, 'Mach Number'),
#     (P_values, 'Pressure (P)'),
#     (rho_values, 'Density (ρ)')
# )
# plot_on_same_axes(x_values, data_sets)



def plot_temperature_data(T_gas_values, Twg_values, Twl_values, Tl_values, x_values, 
                          T_max_working, labels, x_label, y_label, title, file_name):
    """
    Plots four temperature datasets against distance and saves the plot.

    Args:
        T_gas_values (list): Temperature values for the gas.
        Twg_values (list): Temperature values for wall gas.
        Twl_values (list): Temperature values for wall liquid.
        Tl_values (list): Temperature values for the liquid.
        x_values (list): Distance values.
        T_max_working (float): Maximum working temperature.
        labels (list): List of labels for the four datasets.
        x_label (str): Label for the x-axis.
        y_label (str): Label for the y-axis.
        title (str): Title of the plot.
        file_name (str): Name of the file to save the plot as (in 'results' folder).
    """
    # Ensure the 'results' directory exists
    results_dir = "results"
    os.makedirs(results_dir, exist_ok=True)

    # Full path for the file
    file_path = os.path.join(results_dir, file_name)

    # If the file already exists, rename it to include '.backup'
    if os.path.exists(file_path):
        backup_path = file_path + "backup.png"
        os.rename(file_path, backup_path)

    # Plot the temperature values
    plt.figure(figsize=(10, 6))
    plt.plot(x_values, T_gas_values, label=labels[0], linestyle='-', marker='o', color='r', markersize=1)
    plt.plot(x_values, Twg_values, label=labels[1], linestyle='-', marker='o', color='orange', markersize=1)
    plt.plot(x_values, Twl_values, label=labels[2], linestyle='-', marker='o', color='purple', markersize=1)
    plt.plot(x_values, Tl_values, label=labels[3], linestyle='-', marker='o', color='b', markersize=1)

    # Plot the horizontal line for T_max_working
    plt.axhline(y=T_max_working, color='r', linestyle='--', label=f'316L max working temp ({T_max_working} K)')

    # Add grid lines with more increments
    plt.grid(which='both', linestyle='--', linewidth=0.5)

    # Find min and max values for the axes
    x_min, x_max = min(x_values), max(x_values)
    y_min = min(min(T_gas_values), min(Twg_values), min(Twl_values), min(Tl_values), T_max_working)
    y_max = max(max(T_gas_values), max(Twg_values), max(Twl_values), max(Tl_values), T_max_working)

    # Adjust min and max values to nearest factors
    x_min = math.floor(x_min / 10) * 10
    x_max = math.ceil(x_max / 10) * 10
    y_min = math.floor(y_min / 200) * 200
    y_max = math.ceil(y_max / 200) * 200

    # Generate ticks that are factors of 10 for x and 200 for y
    x_ticks = list(range(x_min, x_max + 10, 10))  # Factors of 10
    y_ticks = list(range(y_min, y_max + 200, 200))  # Factors of 200
    plt.xticks(x_ticks)
    plt.yticks(y_ticks)

    # Add labels, title, and legend
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title)
    plt.legend()
    plt.grid(True)

    # Save the plot
    plt.tight_layout()
    #plt.savefig(file_path)
    plt.savefig(file_path, dpi=600)  # or dpi=600 for very high resolution

    plt.close()

# Example usage:
# plot_temperature_data(
#     T_gas_values=[300, 320, 340, 360],
#     Twg_values=[310, 330, 350, 370],
#     Twl_values=[305, 325, 345, 365],
#     Tl_values=[295, 315, 335, 355],
#     x_values=[0, 1, 2, 3],
#     T_max_working=350,
#     labels=["T_gas", "Twg", "Twl", "Tl"],
#     x_label="Distance (m)",
#     y_label="Temperature (K)",
#     title="Temperature Profiles vs Distance",
#     file_name="temperature_plot.png"
# )


def plot_Ts_vs_x(x_values, Ts_values, x1, x2):
    """
    Plot Ts_values vs. x_values, highlight the points closest to x1 and x2, 
    and label their Ts.

    Parameters
    ----------
    x_values : array-like of shape (n,)
        The x-coordinates.
    Ts_values : array-like of shape (n,)
        The temperature values at each x.
    x1, x2 : float
        Two x-positions to locate on the curve.

    Returns
    -------
    idx1, idx2 : int
        Indices of the closest points to x1 and x2, respectively.
    """
    # ensure numpy arrays
    x  = np.asarray(x_values, dtype=float)
    Ts = np.asarray(Ts_values, dtype=float)

    # find the indices of the closest matches
    idx1 = np.argmin(np.abs(x - x1))
    idx2 = np.argmin(np.abs(x - x2))

    x1_match, Ts1 = x[idx1], Ts[idx1]
    x2_match, Ts2 = x[idx2], Ts[idx2]

    # make the plot
    plt.figure()
    plt.plot(x, Ts, '-o', label='Ts vs x', lw=1, ms=4)
    plt.plot(x1_match, Ts1, 'ro', label=f'Closest to {x1}')
    plt.plot(x2_match, Ts2, 'go', label=f'Closest to {x2}')

    # annotate the two points
    plt.annotate(f'{Ts1:.3g}',
                 xy=(x1_match, Ts1),
                 xytext=(5, 5),
                 textcoords='offset points')
    plt.annotate(f'{Ts2:.3g}',
                 xy=(x2_match, Ts2),
                 xytext=(5, -10),
                 textcoords='offset points')

    plt.xlabel('x')
    plt.ylabel('Ts')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    return idx1, idx2

def plot_Ts_vs_x2(x_values, Ts_values, x1, x2, filename):
    """
    Plot Ts_values vs. x_values, highlight the points closest to x1 and x2,
    label their Ts values, and save the figure to results/filename.

    Parameters
    ----------
    x_values : array-like of shape (n,)
        The x-coordinates.
    Ts_values : array-like of shape (n,)
        The temperature values at each x.
    x1, x2 : float
        Two x-positions to locate on the curve.
    filename : str
        Filename (with extension) to save in the 'results' folder.

    Returns
    -------
    idx1, idx2 : int
        Indices of the closest points to x1 and x2, respectively.
    """
    # ensure numpy arrays
    x  = np.asarray(x_values, dtype=float)
    Ts = np.asarray(Ts_values, dtype=float)

    # find the indices of the closest matches
    idx1 = np.argmin(np.abs(x - x1))
    idx2 = np.argmin(np.abs(x - x2))

    x1_match, Ts1 = x[idx1], Ts[idx1]
    x2_match, Ts2 = x[idx2], Ts[idx2]

    # make the plot
    plt.figure()
    plt.plot(x, Ts, '-o', label='Ts vs x', lw=1, ms=4)
    plt.plot(x1_match, Ts1, 'ro', label=f'Point 1: {x1}mm')
    plt.plot(x2_match, Ts2, 'go', label=f'Point 2: {x2}mm')

    # annotate the two points
    plt.annotate(f'{Ts1:.3g}',
                 xy=(x1_match, Ts1),
                 xytext=(5, 5),
                 textcoords='offset points')
    plt.annotate(f'{Ts2:.3g}',
                 xy=(x2_match, Ts2),
                 xytext=(5, -10),
                 textcoords='offset points')

    plt.xlabel('Axial Distance (mm)')
    plt.ylabel('Surface Temp, Ts (K)')
    plt.title('Steady State Surface Temp.')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()

    # ensure 'results' directory exists and save the plot
    os.makedirs('results', exist_ok=True)
    filepath = os.path.join('results', filename)
    plt.savefig(filepath, dpi=600, bbox_inches='tight')

    # show the plot
    plt.show()

    return idx1, idx2

import matplotlib.pyplot as plt
import os

def plot_and_save_multiple_datasets(x_values, datasets, x_label, y_label, title, filename):
    """
    Plots multiple datasets with the same x-values.

    Parameters:
    - x_values: list or array of x values (shared by all datasets)
    - datasets: list of tuples (label, y_values)
        Example: [("Dataset 1", y_values1), ("Dataset 2", y_values2)]
    - x_label: label for the x-axis
    - y_label: label for the y-axis
    - title: plot title
    - filename: name of the file (saved in 'data/' folder)

    Saves the plot as a PNG file in the 'data/' folder.
    """
    plt.figure(figsize=(10, 6))
    
    # Cycle through different colors automatically
    for label, y_values in datasets:
        plt.plot(x_values, y_values, label=label)
    
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title)
    plt.legend()
    plt.grid(True)
    
    # Ensure 'data' directory exists
    os.makedirs('data', exist_ok=True)
    
    # Save the plot
    output_path = os.path.join('results', filename)
    plt.savefig(output_path, dpi=600, format='png')
    plt.close()
    print(f"Plot saved to {output_path}")

# import matplotlib.pyplot as plt
# import os

def plot_multiple_datasets_varied_x(datasets, x_label, y_label, title, filename):
    """
    Plots multiple datasets on the same axes, even if their x-values differ in range and resolution.

    Parameters:
    - datasets: list of dictionaries, each with keys:
        - 'label': string, label for the dataset
        - 'x_values': list or array of x-values
        - 'y_values': list or array of y-values
        - 'color': string, color for the line
        - 'linestyle': string, matplotlib line style (e.g. '-', '--', ':', '-.')
    - x_label: label for the x-axis
    - y_label: label for the y-axis
    - title: plot title
    - filename: name of the file (saved in 'data/' folder)

    Saves the plot as a PNG file in the 'data/' folder.
    """
    plt.figure(figsize=(10, 6))
    
    for dataset in datasets:
        plt.plot(
            dataset['x_values'], 
            dataset['y_values'], 
            label=dataset['label'], 
            color=dataset['color'], 
            linestyle=dataset['linestyle']
        )
    
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title)
    plt.legend()
    plt.grid(True)
    
    # Ensure 'data' directory exists
    os.makedirs('data', exist_ok=True)
    
    # Save the plot
    output_path = os.path.join('data', filename)
    plt.savefig(output_path,dpi=600, format='png')
    plt.close()
    print(f"Plot saved to {output_path}")


def plot_multiple_datasets_varied_x_2(
    datasets, 
    x_label, 
    y_label, 
    title, 
    filename,
    events=None
):
    """
    Plots multiple datasets on the same axes, even if their x-values differ in range and resolution.
    Optionally adds vertical lines for one or more single-point events.

    Parameters:
    - datasets: list of dictionaries, each with keys:
        - 'label': string, label for the dataset
        - 'x_values': list or array of x-values
        - 'y_values': list or array of y-values
        - 'color': string, color for the line
        - 'linestyle': string, matplotlib line style (e.g. '-', '--', ':', '-.')
    - x_label: label for the x-axis
    - y_label: label for the y-axis
    - title: plot title
    - filename: name of the file (saved in 'data/' folder)
    - events: optional list of event dictionaries, each with keys:
        - 'time': float, x-value at which the event occurs
        - 'color': string, color for the vertical line(s)
        - 'linestyle': string, matplotlib line style (e.g. '--', ':', '-.')
        - 'label': string, label for the legend (optional)
    """
    plt.figure(figsize=(10, 6))
    
    for dataset in datasets:
        plt.plot(
            dataset['x_values'], 
            dataset['y_values'], 
            label=dataset['label'], 
            color=dataset['color'], 
            linestyle=dataset['linestyle']
        )

    # Add vertical lines for events if specified
    if events is not None:
        for event in events:
            time = event.get('time', None)
            color = event.get('color', 'gray')
            linestyle = event.get('linestyle', '--')
            label = event.get('label', None)
            
            if time is not None:
                plt.axvline(
                    x=time, 
                    color=color, 
                    linestyle=linestyle, 
                    label=label
                )

    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title)
    plt.legend()
    plt.grid(True)

    #plt.figure(figsize=(10, 6))  # Taller plot

    # Ensure 'data' directory exists
    os.makedirs('data', exist_ok=True)
    
    # Save the plot
    output_path = os.path.join('data', filename)
    plt.savefig(output_path,dpi=600, format='png')
    plt.close()
    print(f"Plot saved to {output_path}")
