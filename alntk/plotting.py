import matplotlib as mpl
import matplotlib.pyplot as plt

def default_plot_style():
    """
    Call function at the beginning of the file in order to set the default style.
    """
    mpl.rcParams['lines.linewidth'] = 1
    mpl.rcParams['font.family'] = 'Arial'
    mpl.rcParams['font.size'] = 12
    
    mpl.rcParams['axes.spines.right'] = False
    mpl.rcParams['axes.spines.top'] = False
    mpl.rcParams['axes.linewidth'] = 0.4
    mpl.rcParams['axes.titlepad'] = 10.
    mpl.rcParams['axes.labelpad'] = 8.
    
    mpl.rcParams['patch.force_edgecolor'] = True
    mpl.rcParams['patch.edgecolor'] = 'black'
    mpl.rcParams['patch.linewidth'] = 0.4
    mpl.rcParams['patch.facecolor'] = 'white'
    
    color_cycle = ['377eb8', 'e41a1c', '4daf4a', '984ea3', 'ff7f00', 'ffff33', 'a65628', 'f781bf', '999999']
    mpl.rcParams["axes.prop_cycle"] = mpl.cycler("color", color_cycle)
    return color_cycle