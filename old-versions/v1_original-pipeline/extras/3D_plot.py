#
# The function to make a 3D plot of the data. Has been taken out of main
# pipeline as matplotlib (numpy) is currently not compatible.
#

#~ import matplotlib.pyplot as plt
#~ import matplotlib.ticker as ticker
#~ from mpl_toolkits.mplot3d import Axes3D
#~ import numpy as np
#~ import matplotlib

def draw_3D_plot(self, out_file, top, out_id):
    """Draws a 3D barchart comparing V,J usage for the top clusters taking
       into account the different Vs used by each cluster.
    """
    # Collect plot data
    plot_data = {}
    v_members = set([])
    j_members = set([])
    for record in self.record_list[:top]:
        clus_prop = self.get_proportion(record)
        # Get J gene name
        j_seg = record.get_top_segment('J')
        j_members.add(j_seg)
        # Get V gene name
        v_seg = record.get_top_segment('V')
        v_members.add(v_seg)
        # Add data to dict
        if not v_seg in plot_data:
            plot_data[v_seg] = {}
        try:
            plot_data[v_seg][j_seg] += clus_prop
        except KeyError:
            plot_data[v_seg][j_seg] = clus_prop
    # Make sorted lists of V's and J's
    v_list = sorted(list(v_members))
    j_list = sorted(list(j_members))
    # Create the axes
    colours_dict = {1: 'blue', 2: 'green', 3: 'red', 4: 'cyan',
                    5: 'magenta', 6: 'yellow', 7: 'orange', 8: 'purple'}
    pattern = re.compile('[A-Z](\d{1,2})')
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    # Select first J segment (Z coordinate)
    for z in range(len(j_list)):
        x_coords = range(len(v_list))
        # Build list of y-values (bar heights)
        y_coords = []
        for x in x_coords:
            y_coords.append(plot_data[v_list[x]].get(j_list[z], 0))
        # Build a list of colours
        colour_list = []
        for x in x_coords:
            colour_list.append(v_colour(v_list[x], pattern, colours_dict))
        # Plot V bars for current J segment
        ax.bar(x_coords, np.array(y_coords), zs=z, zdir='y', alpha=0.8,
               color=colour_list, width=1, linewidth=0.2)
    # Change font family
    font = {'family' : 'serif'}
    matplotlib.rc('font', **font)
    # Set title and change/format labels on plot
    ax.set_title('V and J usage for top {1} clusters\n{0}'.format(
            out_id, len(self.record_list[:top]), fontsize=15))
    ax.w_xaxis.set_major_locator(ticker.MultipleLocator(1))
    # Note: have been trying to get tickers correctly places,
    # ha=cen is correct but va=top is not
    ax.w_xaxis.set_ticklabels(v_list, rotation=70, va='top', ha='center') 
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(6)
    ax.w_yaxis.set_major_locator(ticker.MultipleLocator(1))
    ax.w_yaxis.set_ticklabels(j_list, rotation=-70, va='center', ha='left')
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(6)
    for tick in ax.zaxis.get_major_ticks():
        tick.label.set_fontsize(6)
    ax.set_xlabel('V segments', fontsize=10)
    ax.set_ylabel('J segments', fontsize=10)
    ax.set_zlabel('Proportion of reads', fontsize=10)
    # Save plot as PDF
    fig.savefig(out_file)
    
    return 0 
