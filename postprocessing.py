import matplotlib.pyplot as plt
import numpy as np
import math
import matplotlib.colors as mcolors
from matplotlib.patches import Circle
from matplotlib.gridspec import GridSpec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from capacity import Capacity

class Plot:
    def __init__(self):
        # Set figure size
        self.fig = plt.figure(figsize=(12, 12), constrained_layout=False)

        # GridSpec for flexibility
        # 2 rows and 2 columns
        # - 1st cell: 3D plot
        # - 2nd cell: MM plot
        # - 3rd cell: PM plot for thetas 0°, 90°, 180°, 270°
        # - 4th cell: PM plot for multiple loads
        gs = GridSpec(2, 2, figure=self.fig, 
                    height_ratios=[1, 1],      # top row and bottom row has same height
                    width_ratios=[1, 1],       # left and right column has same width
                    hspace=0.18,               # vertikal space between row (height space)
                    top=0.92,                  # prevent suptitle covering the figure
                    bottom=0.04,
                    wspace=0.18)               # horizontal space between column (width space)

        # 1st cell: 3D Plot
        self.ax_3d = self.fig.add_subplot(gs[0, 0], projection='3d')
        self.ax_3d.set_box_aspect([1, 1, 1])
        self.ax_3d.set_autoscale_on(True)

        # Inset to show the concrete section at the top left of the 3D plot
        self.ax_section = inset_axes(self.ax_3d,
                                width="26%",
                                height="26%",
                                loc='upper left',
                                borderpad=2)

        # 2nd cell: MM plot
        self.ax_mm = self.fig.add_subplot(gs[0,1]) # MM
        self.ax_mm.set_title('MM Interaction Diagram')
        self.ax_mm.set_xlabel('Mnx (KN.m)')
        self.ax_mm.set_ylabel('Mny (KN.m)', labelpad=-5)
        self.ax_mm.axhline(y=0, color='k', linestyle='-')
        self.ax_mm.axvline(x=0, color='k', linestyle='-')

        # 3rd cell: PM plot for thetas 0°, 90°, 180°, 270°
        self.ax_pm_thetas = self.fig.add_subplot(gs[1,0])
        self.ax_pm_thetas.set_title('PM Interaction Diagram\n(0°, 90°, 180°, 270°)')
        self.ax_pm_thetas.set_xlabel('Mn (KN.m)')
        self.ax_pm_thetas.set_ylabel('Pn (KN)', labelpad=-5)
        self.ax_pm_thetas.axhline(y=0, color='k', linestyle='-')

        # 4th cell: PM plot for multiple loads
        self.ax_pm_loads = self.fig.add_subplot(gs[1,1])
        self.ax_pm_loads.set_title('PM Interaction Diagram and Ultimate Loads')
        self.ax_pm_loads.set_xlabel('Mn (KN.m)')
        self.ax_pm_loads.set_ylabel('Pn (KN)', labelpad=-5)
        self.ax_pm_loads.axhline(y=0, color='k', linestyle='-')

        # Set aspect
        self.ax_3d.set_aspect('equal', adjustable='box', anchor='C')
        self.ax_mm.set_aspect('equal', adjustable='box', anchor='C')
        self.ax_pm_thetas.set_aspect('auto', adjustable='box', anchor='C')
        self.ax_pm_loads.set_aspect('auto', adjustable='box', anchor='C')

        # Set title at the top
        self.fig.suptitle('Concrete Column Interaction Diagrams', fontsize=16, y=0.96)

    def draw_section(self, solution):
        """
        Draw the concrete section
        """
        ax = self.ax_section
        ax.plot(np.append(solution.concrete_boundary['x_norm'], solution.concrete_boundary['x_norm'][0]),
                    np.append(solution.concrete_boundary['y_norm'], solution.concrete_boundary['y_norm'][0]), 'k-')
        for xi, yi, d in zip(solution.rebar['x_norm'], solution.rebar['y_norm'], solution.rebar['diameter']):
            radius = d / 2.0                     # karena Circle butuh radius
            circle = Circle((xi, yi), radius=radius, 
                            facecolor='maroon', 
                            edgecolor='black', 
                            linewidth=0.1, 
                            alpha=0.8)
            ax.add_patch(circle)
        ax.set_aspect('equal')

    def plot3d(self, rc: Capacity):
        """
        Plot 3D
        """
        # nominalCapacity = nc.capacity
        reducedCapacity = rc.capacity
        ax = self.ax_3d

        for i, _ in enumerate(reducedCapacity):
            # Skip the first value.
            # ax.fill_between is working by connecting the 2 groups of points.
            # Hence, we start plotting from the second value.
            if i == 0:
                continue

            # First group (collection of Mnx and Mny at Pn1)
            Pn_1st = reducedCapacity[i-1]['Pn']
            Mnx_1st = reducedCapacity[i-1]['Mnx']
            Mny_1st = reducedCapacity[i-1]['Mny']
            # Close the loop, so that there is no gap on the plot surface
            Pn_1st = [*Pn_1st, Pn_1st[0]]
            Mnx_1st = [*Mnx_1st, Mnx_1st[0]]
            Mny_1st = [*Mny_1st, Mny_1st[0]]

            # Second group (collection of Mnx and Mny at Pn2)
            Pn_2nd = reducedCapacity[i]['Pn']
            Mnx_2nd = reducedCapacity[i]['Mnx']
            Mny_2nd = reducedCapacity[i]['Mny']
            # Close the loop, so that there is no gap on the plot surface
            Pn_2nd = [*Pn_2nd, Pn_2nd[0]]
            Mnx_2nd = [*Mnx_2nd, Mnx_2nd[0]]
            Mny_2nd = [*Mny_2nd, Mny_2nd[0]]

            # Plot
            color = 'thistle'
            edgecolor ='crimson'
            alpha = 0.3
            linewidth = 0.3
            if rc.reduced:
                color = 'lightblue'
                edgecolor ='b'
                alpha = 0.5
                linewidth = 0.3
            ax.fill_between(Mnx_1st, Mny_1st, Pn_1st, Mnx_2nd, Mny_2nd, Pn_2nd,
                            alpha=alpha, color=color, edgecolor=edgecolor, linewidth=linewidth)

        # Add labels and title
        ax.set_xlabel('Mnx (KN.m)')
        ax.set_ylabel('Mny (KN.m)')
        ax.set_zlabel('Pn (KN)')
        ax.set_title('3D Interaction Diagram')
        ax.set_aspect('equalxy')

    def plot_PM_from_thetas(self, Mns: list, Pns: list, thetas: float, reduced=False):
        ax = self.ax_pm_thetas
        legends = []
        x_max = 0
        y_max = 0
        y_min = 0

        # For each load and capacity, plot the interaction diagram
        colors = list(mcolors.TABLEAU_COLORS.values())
        linestyle = '--'
        opacity = 0.3
        if reduced:
            linestyle = '-'
            opacity = 1
        for i, (Mn, Pn, theta) in enumerate(zip(Mns, Pns, thetas)):
            # Plot nominal capacity
            ax.plot(Mns[i], Pns[i], c=colors[i], linestyle=linestyle, alpha=opacity)

            # Collect legend
            legends.append(f'{theta:.1f}°')

            # Find the max value of Mn and Pn
            x_max = max(x_max, max(Mn))
            y_max = max(y_max, max(Pn))
            y_min = min(y_min, min(Pn))

        leg = ax.legend(legends, loc='upper right')
        leg.set_draggable(True)
        ax.set_xlim(0, x_max*1.1)
        ax.set_ylim(y_min, y_max)

    def plot_PM_from_loads(self, Mns: list, Pns: list, loads: dict, thetas: float, reduced=False):
        ax = self.ax_pm_loads
        legends = []
        x_max = 0
        y_max = 0
        y_min = 0

        # For each load and capacity, plot the interaction diagram
        colors = list(mcolors.TABLEAU_COLORS.values())
        linestyle = '--'
        opacity = 0.3
        if reduced:
            linestyle = '-'
            opacity = 1
        for i, (Pu, Mux, Muy, safe) in enumerate(zip(loads["Pu"], loads["Mux"], loads["Muy"], loads["safe"])):
            # Calculate vector moment
            Mu = math.sqrt(Mux**2 + Muy**2)

            # Plot nominal capacity
            ax.plot(Mns[i], Pns[i], c=colors[i], linestyle=linestyle, alpha=opacity)

            # Plot ultimate load
            marker = 'o' if safe else 'x'
            ax.plot([Mu], [Pu], c=colors[i], marker=marker)

            # Collect legend
            legends.append(f'{thetas[i]:.1f}°_Mu:{Mu:.1f}_Pu:{Pu:.1f}')
    
            # Find the max value of Mn and Pn
            x_max = max(x_max, max(Mns[i]), Mu)
            y_max = max(y_max, max(Pns[i]), Pu)
            y_min = min(y_min, min(Pns[i]), Pu)

        leg = ax.legend(legends, loc='upper right')
        leg.set_draggable(True)
        ax.set_xlim(0, x_max*1.1)
        ax.set_ylim(y_min, y_max)

    def plot_MM(self, c: Capacity):
        results = c.capacity
        ax = self.ax_mm
        legends = []

        # Find the max value of Mn.
        # We will include this index in the plot.
        i_max = 0
        Mn_max = 0
        for i, result in enumerate(results):
            Mn = max(result['Mnx'])
            if Mn > Mn_max:
                Mn_max = Mn
                i_max = i

        # For each load and capacity, plot the interaction diagram
        colors = list(mcolors.TABLEAU_COLORS.values())
        for i, result in enumerate(results):
            # Skip if i is 0 (full compression) and i is the last index (full tension)
            # since these are only points
            if i == 0 or i == len(results)-1:
                continue
            # Only show these indices to avoid crowded plot
            # 1 & 2: tension area
            # others: selected indices in between compression
            selected_index = {1, 2, 3, 5, 8, 11}
            selected_index.update({i_max})
            if i not in list(selected_index):
                continue

            # Plot nominal capacity
            Pn = result['Pn'][i]
            Mnx = [*result['Mnx'], result['Mnx'][0]]
            Mny = [*result['Mny'], result['Mny'][0]]
            linestyle = '--' if Pn < 0 else '-'

            i_col = i % len(colors)
            ax.plot(Mnx, Mny, c=colors[i_col], linestyle=linestyle)

            # Collect legend
            legends.append(f'Pn:{result["Pn"][i]:.1f}')

        leg = ax.legend(legends, loc='upper right')
        leg.set_draggable(True)

    def show(self):
        self.ax_section.grid(linewidth=0.5)
        self.ax_3d.grid(linewidth=0.5)
        self.ax_pm_thetas.grid(linewidth=0.5)
        self.ax_pm_loads.grid(linewidth=0.5)
        self.ax_mm.grid(linewidth=0.5)
        plt.grid(linewidth=0.5)
        plt.show()
