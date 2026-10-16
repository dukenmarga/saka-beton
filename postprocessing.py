import matplotlib.pyplot as plt
import numpy as np
import math
import os
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

        # Drawings remembered for the report figures
        self._section_solution = None
        self._pm_loads_calls = []

        # Set aspect
        self.ax_3d.set_aspect('equal', adjustable='box', anchor='C')
        self.ax_mm.set_aspect('equal', adjustable='box', anchor='C')
        self.ax_pm_thetas.set_aspect('auto', adjustable='box', anchor='C')
        self.ax_pm_loads.set_aspect('auto', adjustable='box', anchor='C')

        # Set title at the top
        self.fig.suptitle('Concrete Column Interaction Diagrams', fontsize=16, y=0.96)

    def draw_section(self, solution, ax=None):
        """
        Draw the concrete section
        """
        if ax is None:
            ax = self.ax_section
            self._section_solution = solution
        ax.plot(np.append(solution._original_concrete_boundary['x_norm'], solution._original_concrete_boundary['x_norm'][0]),
                    np.append(solution._original_concrete_boundary['y_norm'], solution._original_concrete_boundary['y_norm'][0]), 'k-')
        for i in range(len(solution.openings)):
            ax.plot(np.append(solution.openings[i]['x_norm'], solution.openings[i]['x_norm'][0]),
                        np.append(solution.openings[i]['y_norm'], solution.openings[i]['y_norm'][0]), 'k-')
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
        x_min = 0   # tracks the most negative Mn (left side of diagram)
        x_max = 0   # tracks the most positive Mn (right side of diagram)
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
            prefix = 'ϕ' if reduced else ''
            ax.plot(Mn, Pn, c=colors[i], linestyle=linestyle, alpha=opacity, label=f'{prefix}{theta:.1f}°')

            # Find the max and min values of Mn and Pn
            x_min = min(x_min, min(Mn))
            x_max = max(x_max, max(Mn))
            y_max = max(y_max, max(Pn))
            y_min = min(y_min, min(Pn))

        leg = ax.legend(loc='upper right')
        leg.set_draggable(True)
        ax.set_xlim(x_min*1.1, x_max*1.1)  # 10% margin on both sides
        ax.set_ylim(y_min*1.05, y_max*1.05)

    def plot_PM_from_loads(self, Mns: list, Pns: list, loads: dict, thetas: float, reduced=False, ax=None):
        if ax is None:
            ax = self.ax_pm_loads
            self._pm_loads_calls.append((Mns, Pns, loads, thetas, reduced))
        legends = []
        x_min = 0   # tracks the most negative Mn
        x_max = 0   # tracks the most positive Mn
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
            # Calculate vector moment (always positive — it is the resultant magnitude)
            Mu = math.sqrt(Mux**2 + Muy**2)

            # Plot nominal capacity
            ax.plot(Mns[i], Pns[i], c=colors[i], linestyle=linestyle, alpha=opacity)

            # Plot ultimate load (always on the positive Mn side since Mu is a magnitude)
            marker = 'o' if safe else 'x'
            ax.plot([Mu], [Pu], c=colors[i], marker=marker)

            # Collect legend
            legends.append(f'{thetas[i]:.1f}°_Mu:{Mu:.1f}_Pu:{Pu:.1f}')

            # Find the max and min values of Mn and Pn
            x_min = min(x_min, min(Mns[i]))
            x_max = max(x_max, max(Mns[i]), Mu)
            y_max = max(y_max, max(Pns[i]), Pu)
            y_min = min(y_min, min(Pns[i]), Pu)

        leg = ax.legend(legends, loc='upper right')
        leg.set_draggable(True)
        ax.set_xlim(x_min*1.1, x_max*1.1)  # 10% margin on both sides
        ax.set_ylim(y_min*1.05, y_max*1.05)

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
            Pn = result['Pn'][0]
            Mnx = [*result['Mnx'], result['Mnx'][0]]
            Mny = [*result['Mny'], result['Mny'][0]]
            linestyle = '--' if Pn < 0 else '-'

            i_col = i % len(colors)
            ax.plot(Mnx, Mny, c=colors[i_col], linestyle=linestyle)

            # Collect legend
            legends.append(f'Pn:{result["Pn"][0]:.1f}')

        leg = ax.legend(legends, loc='upper right')
        leg.set_draggable(True)

    def build_triangulated_mesh_3d(self, c: Capacity):
        """
        Render the 3D surface as a proper triangulated mesh using the structured
        k×alpha grid. Adjacent alphas in the same ring and adjacent rings form
        quads → 2 triangles each. Avoids Delaunay topology errors from projecting
        the closed 3D shell onto the Mnx-Mny plane.
        """
        capacity = c.capacity
        n_levels = len(capacity)
        n_alpha  = len(capacity[0]['Mnx'])

        verts_list = []
        for k in range(n_levels):
            result   = capacity[k]
            pn_level = result['Pn'][0]
            for i in range(n_alpha):
                verts_list.append([result['Mnx'][i], result['Mny'][i], pn_level])
        verts = np.array(verts_list, dtype=np.float64)

        faces_list = []
        for k in range(n_levels - 1):
            for i in range(n_alpha):
                j   = (i + 1) % n_alpha
                v00 = k * n_alpha + i
                v01 = k * n_alpha + j
                v10 = (k + 1) * n_alpha + i
                v11 = (k + 1) * n_alpha + j
                faces_list.append([v00, v01, v11])
                faces_list.append([v00, v11, v10])
        faces = np.array(faces_list, dtype=np.int32)

        ax = self.ax_3d
        ax.plot_trisurf(verts[:, 0], verts[:, 1], verts[:, 2],
                        triangles=faces, cmap='viridis', alpha=0.7, linewidth=0)
        ax.set_xlabel('Mnx (KN.m)')
        ax.set_ylabel('Mny (KN.m)')
        ax.set_zlabel('Pn (KN)')
        ax.set_title('3D Interaction Surface')
        ax.set_aspect('equalxy')

    def save_for_report(self, directory: str = "tmp", prefix: str = "report"):
        """
        Draw the section and the PM diagram on figures of their own and save them
        as SVG for the PDF report. This replays the drawings already made for the
        screen figure, so nothing is recalculated.
        """
        os.makedirs(directory, exist_ok=True)
        images = {}

        fig, ax = plt.subplots(figsize=(4, 4))
        self.draw_section(self._section_solution, ax=ax)
        ax.set_xlabel('x (mm)')
        ax.set_ylabel('y (mm)')
        images["section"] = self._save_svg(fig, directory, f"{prefix}_section")

        # No title on the figure itself: the report puts a heading above it.
        fig, ax = plt.subplots(figsize=(7, 5))
        ax.set_xlabel('Mn (KN.m)')
        ax.set_ylabel('Pn (KN)')
        ax.axhline(y=0, color='k', linestyle='-')
        for Mns, Pns, loads, thetas, reduced in self._pm_loads_calls:
            self.plot_PM_from_loads(Mns, Pns, loads, thetas, reduced=reduced, ax=ax)
        images["diagram"] = self._save_svg(fig, directory, f"{prefix}_diagram")

        return images

    def _save_svg(self, fig, directory: str, name: str):
        path = os.path.join(directory, f"{name}.svg")
        # Matplotlib writes an RDF <metadata> block that fpdf2 cannot read and
        # warns about; passing empty metadata leaves it out.
        fig.savefig(path, format="svg", bbox_inches="tight",
                    metadata={"Creator": None, "Date": None, "Format": None, "Type": None})
        plt.close(fig)
        return path

    def show(self):
        self.ax_section.grid(linewidth=0.5)
        self.ax_3d.grid(linewidth=0.5)
        self.ax_pm_thetas.grid(linewidth=0.5)
        self.ax_pm_loads.grid(linewidth=0.5)
        self.ax_mm.grid(linewidth=0.5)
        plt.grid(linewidth=0.5)
        plt.show()
