#plotting
import matplotlib.pyplot as plt

def plot_chamber_nozzle_geometry(x_values,r_values):
    plt.plot(x_values, r_values, 'ko', markersize=2)
    plt.plot(x_values, -r_values, 'ko', markersize=2)

    plt.grid(True)
    plt.axis('equal')
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
        plt.plot(x_values, y_values, color + 'o', markersize=2, label=f'{y_label} vs {x_label}')
        
        # Add labels and grid
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.grid(True)
        plt.legend()
        plt.title(f'{y_label} vs {x_label}')
        
        # Show the plot
        plt.show()
