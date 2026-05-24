import numpy as np
import pandas as pd
import datetime as dt
import multiprocessing
import shapely
import math
import pygmsh
import copy
from scipy.optimize import brentq
from shapely import LineString, Polygon, Point
from shapely.affinity import rotate
from typing import Literal

from sni_2847 import SNI_2847_2019
import geometry
from postprocessing import Plot
from material import Material
from capacity import Capacity

import faulthandler
faulthandler.enable()

class ConcreteInteractionDiagram:
    """
    Calculate concrete interaction diagram (biaxial analysis)
    """
    def __init__(
            self, concrete_boundary: pd.DataFrame,
            rebar: pd.DataFrame,
            fc_ = 30, # MPa
            fy = 400, # MPa
            fsu = None, # MPa
            epsilon_0 = 0.002,
            epsilon_u = 0.0038,
            gamma = 0.15,
            Es = 200000, # MPa 
            eps_sh = 0.015, # 0.005 - 0.025, depend on material
            eps_su = 0.05,
            transverse: Literal["ties", "spiral"] = "ties",
            loads: dict = {
                "Pu": [0],
                "Mux": [0],
                "Muy": [0],
            },
        ):
        # Define concrete boundary and rebar coordinates
        self._concrete_boundary = concrete_boundary
        self._rebar = rebar

        # Define material properties
        self.material = Material(fc_, fy, fsu, epsilon_0, epsilon_u, gamma, Es, eps_sh, eps_su)
        self.sni = SNI_2847_2019(transverse)

        # Define ultimate loads
        self.loads = loads

    def preanalysis(self):
        """
        Preanalysis to find the geometric centroid of the concrete boundary and transform the coordinates relative to the centroid
        """

        # Calculate geometric centroid
        # When defining the coordinate of concrete boundary above, we randomly choose origin O to be (0, 0),
        # so that we can define the concrete boundary relative to the origin.
        # But, for stability & convergence reasons, we need to redefine the "correct" origin of the concrete boundary.
        # We will be using the geometric centroid of the concrete boundary as the new origin.
        # Other method is by using plastic centroid (not used here).
        XO, YO = self.calculate_centroids(self.concrete_boundary[['x', 'y']].values)

        # Transform concrete and rebar coordinates relative to
        # the new origin
        self.transform_to_new_origin(XO, YO)

        # Close loop
        self.close_boundary_loop()

        self.transform_to_O()

    @property
    def concrete_boundary(self):
        return self._concrete_boundary

    @property
    def rebar(self):
        return self._rebar

    def _brent_method(self, Pn_target: float, alpha: float = 0,
                      h_min: float = -100, h_max: float = 100,
                      tol: float = 1e-1, max_iter: float = 100):
        """
        Find combination of Mnx and Mny, such that Pn is close to Pn_target
            by utilizing Brent method.
        """
        d, r = brentq(
            f=self.find_Pn,
            a=h_min,
            b=h_max,
            args=(alpha, Pn_target),
            full_output=True,
            disp=False,
            xtol=tol,
            rtol=tol,
            maxiter=max_iter
        )
        h = float(d)
        return self.find_Pn_Mn(h, alpha)

    def calculate_centroids(self, points: np.ndarray):
        """
        Calculate gross section's geometric centroid.
            We use pygmsh to mesh the concrete boundary
            in aims to calculate the centroid. No need fine meshing.
        """
        if self.first_and_last_node_same(points):
            points = points[:-1, :]

        with pygmsh.geo.Geometry() as geom:
            poly = geom.add_polygon(
                points,
                mesh_size=5000,
            )
            geom.add_physical(poly, label="surface")
            mesh = geom.generate_mesh()

        vertices = mesh.points
        triangles = None
        for cell in mesh.cells:
            if cell.type == "triangle":
                triangles = cell.data
                break

        # centroid will be the centroid of each meshed triangle
        centroid = np.zeros((len(triangles), 2))
        # areas will be the area of each meshed triangle
        areas = np.zeros(len(triangles))

        # For each triangle, calculate the centroid and area
        for i in range(len(triangles)):
            # Calculate area
            x1 = vertices[triangles[i][0]][0]
            x2 = vertices[triangles[i][1]][0]
            x3 = vertices[triangles[i][2]][0]
            y1 = vertices[triangles[i][0]][1]
            y2 = vertices[triangles[i][1]][1]
            y3 = vertices[triangles[i][2]][1]
            area = 0.5 * np.abs(x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2))

            x_sum = x1 + x2 + x3
            y_sum = y1 + y2 + y3

            # Calculate centroid
            x_cen = x_sum / 3
            y_cen = y_sum / 3

            # Assign values
            centroid[i] = [x_cen, y_cen]
            areas[i] = area

        # Calculate total area
        total_area = np.sum(areas)

        # Calculate global centroid
        X_cen = np.sum(centroid[:,0] * areas) / total_area
        Y_cen = np.sum(centroid[:,1] * areas) / total_area

        return X_cen, Y_cen


    def calculate_depth(self, na_line: LineString):
        """
        Calculate depth (height) of compression zone
        from the neutral axis to the most compressed point
        (i.e. the top (extreme) of the concrete boundary).
        point_on_na: x and y coordinate of a point on neutral axis
        alpha: angle of the neutral axis that passes through point_on_na (clockwise)

        point_on_na is X
        alpha is a

             NA
              \
               \
              a \
            -----X----
                  \
                   \
        """
        # Localising
        concrete_boundary = self.concrete_boundary
        XO = self.XO
        YO = self.YO

        # Find midline along neutral axis and set it as local origin
        mid_na = na_line.interpolate(0.5, True)

        # Extract x,y coordinates
        x_na, y_na = mid_na.x, mid_na.y

        # Map concrete boundary from global to local coordinate
        # (relative to neutral axis).
        # Local origin O' is set to the midline of the neutral axis
        concrete_boundary_local = geometry.map_to_local(concrete_boundary[['x', 'y']].values, XO, YO, x_na, y_na, na_line)

        # Get Xe, Ye: the farthest point on the compression side.
        # tens_node: farthest node on tension side (not always true, but it is the least value)
        # comp_node: farthest node on compression side
        # Initial assumption is the first is in tension, the second is in compression
        tens_node, comp_node = geometry.find_farthest_node(concrete_boundary_local)
        # xec, yec: x and y coordinate of the farthest compression point (local coordinate)
        xec, yec = comp_node[0], comp_node[1]
        # xet, yet: x and y coordinate of the farthest tension point (local coordinate)
        xet, yet = tens_node[0], tens_node[1]

        # yec must be positive, since compression height is positive above the neutral axis.
        # If yec negative, then our initial assumption is wrong. We just need to swap
        # compression and tension.
        if yec < 0:
            xec, yec, xet, yet = xet, yet, xec, yec

        # Make a point on the farthest compression point and tension (local)
        xeye_c = np.array([[xec, yec]])
        xeye_t = np.array([[xet, yet]])

        # Map the points from local to global
        XeYe_c = geometry.map_to_global(xeye_c, XO, YO, x_na, y_na, na_line)
        XeYe_t = geometry.map_to_global(xeye_t, XO, YO, x_na, y_na, na_line)

        # Xe, Ye is the farthest compression point in global coordinate
        Xec, Yec = XeYe_c[0][0], XeYe_c[0][1]
        Xet, Yet = XeYe_t[0][0], XeYe_t[0][1]

        return xec, yec, Xec, Yec, xet, yet,Xet, Yet

    def calculate_strain(self, rebar_y: np.ndarray, y_comp: float):
        """
        Calculate strain for each rebar.
            rebar_y: y coordinate of rebar (local coordinate system)
            y_comp: maximum depth of concrete compression
            epsilon_u: ultimate strain of concrete (usually 0.0038)
        """
        if y_comp == 0:
            raise ValueError("concrete compression depth cannot be zero")

        return rebar_y / y_comp * self.material.epsilon_u

    def close_boundary_loop(self):
        """
        Check if the boundary is closed.
        If not, close it
        """
        if not self.first_and_last_node_same(self.concrete_boundary[['x', 'y']].values):
            first_point = self.concrete_boundary.loc[[0]]
            self._concrete_boundary = pd.concat([self.concrete_boundary, first_point], axis=0, ignore_index=True)

    def construct_neutral_axis(self, point_on_na: tuple[float], alpha: float) -> LineString:
        # Localising
        concrete_boundary = self.concrete_boundary
        # Extract x,y
        x_na, y_na = point_on_na
        xmin = concrete_boundary['x'].min()
        xmax = concrete_boundary['x'].max()

        # Create a horizontal line that passes through point_on_na
        # We extend the line to cover the whole concrete boundary on the x-axis
        # (using xmin and xmax)
        line = LineString([(xmin, y_na), (xmax, y_na)])
        # Rotate the line about the point_on_na by alpha
        # We negate alpha, since rotate operation rotates counter clockwise,
        # but we want to rotate clockwise (our calculation is based on clockwise rotation)
        rotated_line = rotate(line, alpha, origin=(x_na, y_na))
        x, y = rotated_line.xy
        # delta_x and delta_y are used to extend the line.
        # The intent is to get the line that covers the whole concrete boundary on both axis
        delta_x = 10*(x[1] - x[0])
        delta_y = 10*(y[1] - y[0])
        # Create a Neutral axis that rotates to alpha that passes through the point_on_na
        na_line = LineString([(x[0] - delta_x, y[0] - delta_y), (x[1] + delta_x, y[1] + delta_y)])
        return na_line

    def find_lowest_left_node(self):
        """
        Find the lowest left node
        """
        min_x = min(self.concrete_boundary['x'])
        min_y = min(self.concrete_boundary['y'])
        
        return min_x, min_y

    def _find_Mnx_Mny_at_Pn(self, Pn_target, alpha, h_min, h_max, method="secant", tol=1e-1, max_iter=100):
        """
        Find Mnx and Mny such that Pn is close to Pn_target
        """
        if method == "brent":
            tol = min(1e-4, tol)
            return self._brent_method(Pn_target, alpha, h_min, h_max, tol, max_iter)
        if method == "secant":
            return self._secant_method(Pn_target, alpha, h_min, h_max, tol, max_iter)
        raise ValueError(f"Unknown method: {method}")

    def find_Mx_My_at_P(self, P: float, alpha: float):
        """
        TODO: Need to rename this function. It will be used to
        calculate the phi Mn and phi Pn as well.
        Find combination of Mx and My that gives Pn close to P and alpha
        """
        # Construct neutral axis that passes through point (0, 0)
        na = self.construct_neutral_axis(point_on_na=(0, 0), alpha=alpha)

        # Get maximum depth from neutral axis to the farthest point of compression and tension
        _, yec, _, _, _, yet, _, _ = self.calculate_depth(na)

        # h_min and h_max are limit where we will iterate to find Mnx and Mny.
        # h is measure from top of compression zone.
        # h = 5 mm means the distance from top of compression zone to the neutral axis.
        # Avoid h at the farthest point of compression and tension
        # to avoid numerical error (no compression or tension area is formed).
        # Hence we offset 1mm from top and bottom   
        h_total = yec + (-yet)
        h_min =  1
        h_max = h_total - 1

        Pn, Mnx, Mny, eps_t = self._find_Mnx_Mny_at_Pn(P, alpha, h_min, h_max)

        # Calculate theta, the angle between Mnx and Mny
        theta = math.atan2(Mny, Mnx) * 180 / math.pi

        # Normalize theta
        theta = self.normalize_theta(theta)

        # eps_t = [0] * len(alphas)
        eps_y = self.material.fy / self.material.Es
        sni = self.sni
        phi = sni.compute_phi(reduction_type="moment_axial", eps_t=eps_t, eps_y=eps_y)

        return Pn, Mnx, Mny, theta, phi

    def find_Pn(self, h: float = 10, alpha: float = 0, Pn_target: float = 0):
        Pn, _, _, _ = self.find_Pn_Mn(h, alpha)
        return Pn - Pn_target

    def find_Pn_Mn(self, h: float = 10,alpha: float = 0):
        """
        Calculate Pn, Mnx, Mny from alpha and h
            alpha is in degrees
            h is height from origin O to neutral axis
            c is depth of compression from neutral axis to the most compressed point
            h + c = c_max
        """
        XO = self.XO
        YO = self.YO

        # Construct benchmark line from a (0,0) that rotates to alpha
        helper_line = self.construct_neutral_axis(point_on_na=(XO, YO), alpha=alpha)

        # Calculate:
        # xec: x coordinate of the most compressed point (local axis)
        # yec: y coordinate of the most compressed point (local axis)
        # yec = compression depth (same as c)
        # Xec: x coordinate of the most compressed point (global axis)
        # Yec: y coordinate of the most compressed point (global axis)
        xec, yec, Xec, Yec, _, _, _, _ = self.calculate_depth(na_line=helper_line)

        # This is a point in neutral axis that is h units away from the most compressed point
        # It still uses local coordinates.
        na_point = np.array([[xec, yec-h]])

        # Map from local to global.
        # The NA point uses global coordinates.
        na_point = geometry.map_to_global(na_point, XO, YO, XO, YO, helper_line)

        # Construct neutral axis (away h units from xe, ye)
        na_line = self.construct_neutral_axis(point_on_na=(na_point[0][0], na_point[0][1]), alpha=alpha)

        # Isolate concrete boundary and rebar locally to this function only
        # - Readability
        concrete_boundary = copy.copy(self.concrete_boundary)
        rebar = copy.copy(self.rebar)

        # Make a polygon from the concrete boundary
        concrete_polygon = Polygon(concrete_boundary[['x', 'y']].values)

        # Calculate coordinate O' (Xa, Ya) in global axis
        mid_na_line = na_line.interpolate(0.5, True)
        Xa = mid_na_line.x
        Ya = mid_na_line.y

        # Remap concrete boundary and rebar coordinates in the rotated coordinate system (Neutral Axis) using new O'.
        rebar_local = geometry.map_to_local(rebar[['x', 'y']].values, XO, YO, Xa, Ya, na_line)

        rebar_local_pd = pd.DataFrame(rebar_local, columns=['x', 'y'])

        # Get strain (epsilon_s) for each point on the rebar, linearly distributed from Xe, Ye.
        # We use h, since it the location of the neutral axis
        rebar_strain_local = self.calculate_strain(rebar_local_pd[['y']].values, h)
        # Concat strain (relative to Neutral Axis line)
        rebar['strain'] = rebar_strain_local

        # Calculate area of each rebar
        rebar_area_As = geometry.calculate_rebar_area(rebar[['diameter']])
        # Concat rebar area
        rebar['As'] = rebar_area_As

        # Get compression part from the concrete if all the section is crossed by
        # the neutral axis
        # Xe, Ye: indicator of the most compressed point, so it knows which part
        # is in compression
        compression_parts = geometry.get_compression_part(concrete_polygon, na_line, (Xec, Yec))

        # TODO: Whitney coefficient is determined by the fc'
        beta = 0.85

        Pcn = 0
        Mcnx = 0
        Mcny = 0

        # Calculate Pn, Mnx, Mny for each compression part 
        for compression_parts in compression_parts:
            Pcn_part, Mcnx_part, Mcny_part = self.Pcn_Mcnx_Mcny(compression_parts, XO, YO)
            Pcn +=  beta * Pcn_part
            Mcnx +=  beta * Mcnx_part
            Mcny +=  beta * Mcny_part

        ##########
        # STRESS #
        ##########

        # Calculate rebar stress based on strain and area
        rebar_fs = self.material.steel_rebar_stress_strain_ufunc(rebar[['strain']])
        rebar['fs'] = rebar_fs

        # Calculate concrete stress at reinforcement rebar (to correct hole made by rebar)
        rebar_fci = self.material.concrete_stress_strain_ufunc(rebar[['strain']])
        rebar['fci'] = rebar_fci

        #########
        # FORCE #
        #########

        # Rebar force: Fs = As * fs
        # As * fs (mm2 * MPa / 1000 -> KN)
        rebar_As_fs = rebar[['fs']].values * rebar[['As']].values / 1000
        rebar['As_fs'] = rebar_As_fs

        # Concrete force at rebar position (to correct hole made by rebar)
        # Fci = As * fci
        # As * fci (mm2 * MPa / 1000 -> KN)
        rebar_As_fci = rebar[['fci']].values * rebar[['As']].values / 1000
        rebar['As_fci'] = rebar_As_fci

        #############
        # SUMMATION #
        #############

        # Total nominal concrete capacity
        # As * fci -> KN
        sum_As_fci = np.sum(rebar[['As_fci']])

        # Total nominal concrete moment capacity about global X-axis
        # As * fci * Ysi -> KN * mm / 1e6 = KN.m
        sum_As_fci_Ysi = np.sum(rebar[['As_fci']].values * rebar[['y']].values / 1e3)

        # Total nominal concrete moment capacity about global Y-axis
        # As * fci * Xsi -> KN * mm / 1e6 = KN.m
        sum_As_fci_Xsi = np.sum(rebar[['As_fci']].values * rebar[['x']].values / 1e3)

        # Total nominal rebar capacity
        # As * fs -> KN
        sum_As_fs = np.sum(rebar[['As_fs']])

        # Total nominal rebar moment capacity about global X-axis 
        # As * fs * Ysi -> KN * mm / 1e6 = KN.m
        sum_As_fs_Ysi = np.sum(rebar[['As_fs']].values * rebar[['y']].values / 1e3)

        # Total nominal rebar moment capacity about global Y-axis
        # As * fs * Xsi -> KN * mm / 1e6 = KN.m
        sum_As_fs_Xsi = np.sum(rebar[['As_fs']].values * rebar[['x']].values / 1e3)

        # Total nominal capacity
        # Pn: Nominal axial capacity
        Pn = Pcn + sum_As_fs - sum_As_fci
        # Mnx: Nominal moment capacity about global X-axis
        MnX = Mcnx + sum_As_fs_Ysi - sum_As_fci_Ysi
        # Mny: Nominal moment capacity about global Y-axis
        MnY = Mcny + sum_As_fs_Xsi - sum_As_fci_Xsi

        # Get strain of the most tensile rebar
        eps_t = -min(rebar[['strain']].values)[0]

        # Return nominal capacity based on selected alpha and h
        return Pn, MnX, MnY, eps_t

    def first_and_last_node_same(self, points: np.ndarray = None):
        """
        Check if the first and last node are the same
        """
        first_point = points[0]
        last_point = points[-1]
        return (first_point == last_point).all()

    def _secant_method(self, Pn_target: float, alpha: float = 0,
                      h_min: float = -100, h_max: float = 100,
                      tol: float = 1e-1, max_iter: float = 100):
        """
        Find combination of Mnx and Mny, such that Pn is close to Pn_target
            by utilizing secant method.
        """
        iter_count = 0
        h1 = h_min
        h2 = h_max
        h3 = 0
        Pn1, _, _, _ = self.find_Pn_Mn(h1, alpha)

        max_iter = 300
        while iter_count < max_iter:
            iter_count += 1
            error = False

            try:
                Pn2, Mnx2, Mny2, eps_t = self.find_Pn_Mn(h2, alpha)
            except ValueError:
                error = True

            # Avoid division by zero
            if abs(Pn2 - Pn1) < 1e-6:
                error = True

            # If h1 and h2 are out of range (getting bigger for each iteration),
            # it indicates that the method is not converging.
            # The possibility is initial h1 and h2 are wrong.
            # Hence, we need to swap h1 and h2 (based on h_min and h_max).
            # We choose 5_000 as limit, since physically it is very unlikely
            # to have a concrete dimension greater than 5_000 mm
            elif abs(h1) > 5_000 or abs(h2) > 5_000:
                error = True

            # If there is an error, we are attempting to swap h1 and h2
            if error:
                # Reset and swap values (opposite of initial value before looping)
                h1 = h_max
                h2 = h_min
                Pn1, _, _, _ = self.find_Pn_Mn(h1, alpha)
                continue

            # Update the new value by using secant formula:
            # h_new = h2 - [f(h2) - target] * (h2 - h1) / (f(h2) - f(h1))
            h3 = h2 - (Pn2 - Pn_target) * (h2 - h1) / (Pn2 - Pn1)

            # Return the result if the iteration has reached the convergence
            if abs(Pn2 - Pn_target) < tol:
                return Pn2, Mnx2, Mny2, eps_t

            # Update the values for next iteration if there is no error
            # and still not reaching the convergence
            h1, Pn1 = h2, Pn2
            h2 = h3

        return 0, 0, 0, 0

    def interpolate_PM(self, theta_target: float, results: Capacity, element: list = ["Pn", "Mnx", "Mny"]):
        Mn = []
        Pn = []
        Pn_string = element[0]
        Mnx_string = element[1]
        Mny_string = element[2]

        for k, result in enumerate(results.capacity):
            for i, r in enumerate(result["theta"]):
                j = (i + 1) % len(result["theta"])

                bigger = result["theta"][i]
                smaller = result["theta"][j]

                in_between = smaller <= theta_target < bigger
                if in_between:
                    # Create a virtual line between origin O to theta.
                    # We choose x = 5_000, as hard limit as physically it is
                    # rare to have concrete section larger than 5_000 mm
                    if 0 < theta_target < 90 or 180 < theta_target < 270:
                        x_inf = 5_000
                    else:
                        x_inf = -5_000
                    y_inf = math.tan(theta_target / 180 * math.pi) * x_inf

                    # Special cases if theta is 0, 90, 180, 270.
                    # At this degree, tan is undefined, hence we need to
                    # handle this situation and set a large value.
                    if theta_target == 0:
                        x_inf = 5_000
                        y_inf = 0
                    elif theta_target == 180:
                        x_inf = -5_000
                        y_inf = 0
                    elif theta_target == 90:
                        x_inf = 0
                        y_inf = 5_000
                    elif theta_target == 270:
                        x_inf = 0
                        y_inf = -5_000
                    if theta_target in [0, 90, 180, 270]:
                        result[Mnx_string][i] = 0 if abs(result[Mnx_string][i]) < 0.000001 else result[Mnx_string][i]
                        result[Mny_string][i] = 0 if abs(result[Mny_string][i]) < 0.000001 else result[Mny_string][i]
                        result[Mnx_string][j] = 0 if abs(result[Mnx_string][j]) < 0.000001 else result[Mnx_string][j]
                        result[Mny_string][j] = 0 if abs(result[Mny_string][j]) < 0.000001 else result[Mny_string][j]

                    # Build a line from failure surface
                    x = result[Mnx_string][i], result[Mnx_string][j]
                    y = result[Mny_string][i], result[Mny_string][j]
                    failure_line = LineString([(x[0], y[0]), (x[1], y[1])])

                    # Build a line that form theta target
                    theta_line = LineString([(0, 0), (x_inf, y_inf)])

                    # Get the intersection.
                    # This intersection is the Mx and My at theta target
                    intersection = failure_line.intersection(theta_line)

                    # extract Mx and My
                    Mx, My = intersection.x, intersection.y

                    # Calculate vector moment from Mx and My
                    M = math.sqrt(Mx**2 + My**2)
                    Mn.append(M)
                    Pn.append(result[Pn_string][i])
                    break
            else:
                Mn.append(0)
                Pn.append(result[Pn_string][i])
        return Mn, Pn

    def interpolate_PM_at_loads(self, loads: dict, results: list, reduced=False):
        # Mns and Pns is collection of Mn and Pn.
        # Example value:
        # Mns = [
        #   [10, 20, 30, ...], -> Mn at 0°
        #   [20, 30, 40, ...], -> Mn at 30°
        #   [30, 40, 50, ...], -> Mn at 60°
        # ]
        # theta = [0, 30, 60]
        Mns = []
        Pns = []
        thetas = []

        # Assume all loads are unsafe
        loads.update({"safe": [False] * len(loads["Pu"])})

        for i, _ in enumerate(loads["Pu"]):
            Pu = loads["Pu"][i]
            Mux = loads["Mux"][i]
            Muy = loads["Muy"][i]

            theta = math.atan2(Muy, Mux) * 180 / math.pi

            # Normalize theta
            theta = self.normalize_theta(theta)
            
            Mu = math.sqrt(Mux**2 + Muy**2)

            # Mn & Pn is a list of points that forms the failure surface at theta
            Mn, Pn = self.interpolate_PM(theta, results)
            Mns.append(Mn)
            Pns.append(Pn)
            thetas.append(theta)

            # If the capacity is reduced (has been multiplied by phi),
            # then we check if the capacity is still safe for the load
            if results.reduced:
                polygon = Polygon(np.column_stack((Mns[-1], Pns[-1])))
                point = Point(Mu, Pu)
                loads["safe"][i] = polygon.contains(point)

        return Mns, Pns, thetas, loads

    def interpolate_PM_at_thetas(self, thetas: list, results: list):
        # Mns and Pns is collection of Mn and Pn.
        # Example value:
        # Mns = [
        #   [10, 20, 30, ...], -> Mn at 0°
        #   [20, 30, 40, ...], -> Mn at 30°
        #   [30, 40, 50, ...], -> Mn at 60°
        # ]
        # theta = [0, 30, 60]
        Mns = []
        Pns = []
        for i, theta in enumerate(thetas):
            # Mn & Pn is a list of points that forms the failure surface at theta
            Mn, Pn = self.interpolate_PM(theta, results)
            Mns.append(Mn)
            Pns.append(Pn)
        return Mns, Pns

    def normalize_theta(self, theta):
        # Handle very small angles around 0 degrees
        if abs(theta) < 0.0000001:
            theta = 0
        
        # Handle negative angles.
        # atan2 returns angle in [-pi, pi], but we want angle in [0, 360]
        # For 0°, we convert it to 360°, so that the final theta (on parent function) is descending from 360° to >0°
        if theta <= 0:
            theta += 360
        return theta

    def Pc_max(self):
        """
        Nominal compressive strength
        """
        fc_ = self.material.fc_
        fy = self.material.fy
        concrete_area = geometry.calculate_concrete_area(self.concrete_boundary[['x', 'y']].values)
        rebar_area = np.sum(geometry.calculate_rebar_area(self.rebar['diameter'].values))
        # TODO: Whitney coefficient is determined from the fc'
        beta = 0.85
        return (beta * 0.85 * fc_ * (concrete_area-rebar_area) + fy * rebar_area) / 1e3

    def Pcn_Mcnx_Mcny(self, concrete_compression_boundary_xy: np.ndarray, XO: float, YO: float):
        """
        Calculate concrete nominal compression axial capacity, Pn, Mnx, Mny based on
        Whitney stress block.
        
            concrete_boundary_xy: coordinates of concrete boundary (global axis)
            XO, YO: origin to calculate moment length
        """
        Pn: float = 0
        Mnx: float = 0
        Mny: float = 0

        # Find x and y centroid of points in global axis
        x_cen_compress, y_centroid_compress = self.calculate_centroids(concrete_compression_boundary_xy)

        # Calculate compression area
        compression_area = geometry.calculate_concrete_area(concrete_compression_boundary_xy)

        # Calculate moment length (distance from centroid to neutral axis)
        x_depth = x_cen_compress - XO
        y_depth = y_centroid_compress - YO

        # Calculate Pn, Mnx, Mny
        # 0.85 is reduction due to size effect based on Hognestad research
        Pn = 0.85 * compression_area * self.material.fc_ / 1e3
        Mnx = 0.85 * compression_area * self.material.fc_ * y_depth / 1e6
        Mny = 0.85 * compression_area * self.material.fc_ * x_depth / 1e6

        return Pn, Mnx, Mny

    def Pt_max(self):
        """
        Calculate maximum concrete tensile strength.
        The value will be in negative.
        """
        rebar_area = np.sum(geometry.calculate_rebar_area(self.rebar['diameter'].values))
        Pn_max = self.material.fy * rebar_area / 1e3
        return -Pn_max

    def solve_axial_from_thetas(self, thetas: list, results: list):
        # Mns and Pns is collection of Mn and Pn.
        # Example value:
        # Mns = [
        #   [10, 20, 30, ...], -> Mn at 0°
        #   [20, 30, 40, ...], -> Mn at 30°
        #   [30, 40, 50, ...], -> Mn at 60°
        # ]
        # theta = [0, 30, 60]
        Mns = []
        Pns = []
        for i, theta in enumerate(thetas):
            # Mn & Pn is a list of points that forms the failure surface at theta
            Mn, Pn = self.interpolate_PM(theta, results)
            Mns.append(Mn)
            Pns.append(Pn)
        return Mns, Pns

    def solve_biaxial_analysis(self):
        """
        Solve biaxial analysis from full tension to full compression
        """
        alphas = [alpha for alpha in range(0, 360, 10)]
        sni = self.sni

        # Supress scientific notation
        np.set_printoptions(suppress=True)
        Pn_max = self.Pc_max()
        Pn_min = self.Pt_max()

        # Add result for full tension.
        # At this level, Mnx = Mny = 0, Pn = Pn_min
        zeros = [0] * len(alphas)
        Pnt = [Pn_min] * len(alphas)

        state = ['tension'] * len(alphas)
        phis = sni.compute_phi_vector(reduction_type=state)

        # Nominal capacity and reduced capacity
        nc = Capacity(False, [])
        rc = Capacity(True, [])

        nc.capacity = [{
            'Pn': Pnt,
            'Mnx': zeros,
            'Mny': zeros,
            'theta': zeros,
        }]
        rc.capacity = [{
            'Pn': phis * Pnt,
            'Mnx': zeros,
            'Mny': zeros,
            'theta': zeros,
        }]

        # Add result between full tension and full compression.

        # Set Pn for extension zone (Pn < 0)
        Pns_target = [P for P in range(int(Pn_min*0.6), 0, abs(int(Pn_min/3)))]
        # Set Pn for compression zone (Pn >= 0)
        Pns_target = Pns_target + [P for P in range(0, int(Pn_max*0.85), int(Pn_max/15))]

        for _, Pn_target in enumerate(Pns_target):
            Pn0 = []
            Mnx0 = []
            Mny0 = []
            theta0 = []
            phi0 = []
            # eps_t0 = []

            with multiprocessing.Pool(processes=16) as pool:
                # Map the function to the inputs and get results in order
                for Pn, Mnx, Mny, theta, phi in pool.starmap(self.find_Mx_My_at_P, [(Pn_target, alpha) for alpha in alphas]):
                    Pn0.append(Pn)
                    Mnx0.append(Mnx)
                    Mny0.append(Mny)
                    theta0.append(theta)
                    phi0.append(phi)
            nc.capacity.append({
                "Pn": Pn0,
                "Mnx": Mnx0,
                "Mny": Mny0,
                "theta": theta0,
            })

            phix = sni.compute_phi(reduction_type="compression", eps_t=0, eps_y=self.material.fy / self.material.Es)
            phiPn_max = sni.compute_phiPn_max(phix * Pn_max)
            rc.capacity.append({
                "Pn": [sni.compute_phiPn(phi, Pn, phiPn_max) for phi, Pn in zip(phi0, Pn0)],
                "Mnx": [phi * Mnx for phi, Mnx in zip(phi0,  Mnx0)],
                "Mny": [phi * Mny for phi, Mny in zip(phi0, Mny0)],
                "theta": theta0,
            })

        # Add result for full compression
        # At this level, Mnx = Mny = 0, Pn = Pn_max
        Pnc = [Pn_max] * len(alphas)

        eps_t = [0] * len(alphas)
        eps_y = [self.material.fy / self.material.Es] * len(alphas)
        state = ["compression"] * len(alphas)
        phis = sni.compute_phi_vector(reduction_type=state, eps_t=eps_t, eps_y=eps_y)

        nc.capacity.append({
            "Pn": Pnc,
            "Mnx": zeros,
            "Mny": zeros,
            "theta": zeros,
        })
        phix = sni.compute_phi(reduction_type="compression", eps_t=0, eps_y=self.material.fy / self.material.Es)
        phiPn_max = sni.compute_phiPn_max(phix * Pn_max)
        rc.capacity.append({
            "Pn": [sni.compute_phiPn(phi, Pn, phiPn_max) for phi, Pn in zip(phis, Pnc)],
            "Mnx": zeros,
            "Mny": zeros,
            "theta": zeros,
        })

        return nc, rc

    def transform_to_new_origin(self, XO: float = 0, YO: float = 0):
        """
        Transform (translate) any coordinate relative to the origin XO, YO.

        Example:
            Coordinate (10, 20) if transformed to new origin (5, 5) will be (5, 15)
        """
        # Translate concrete boundary relative to the origin
        self._concrete_boundary['x'] = self._concrete_boundary['x'] - XO
        self._concrete_boundary['y'] = self._concrete_boundary['y'] - YO

        # Translate concrete rebar relative to the origin
        self._rebar['x'] = self._rebar['x'] - XO
        self._rebar['y'] = self._rebar['y'] - YO

        # Reset XO, YO to 0 after coordinate transformation
        self.XO = 0
        self.YO = 0
        self.Origin = shapely.geometry.Point(XO, YO)

    def transform_to_O(self):
        x, y = self.find_lowest_left_node()
        self._concrete_boundary['x_norm'] = self._concrete_boundary['x'] - x
        self._concrete_boundary['y_norm'] = self._concrete_boundary['y'] - y

        # Translate concrete rebar relative to the origin
        self._rebar['x_norm'] = self._rebar['x'] - x
        self._rebar['y_norm'] = self._rebar['y'] - y

def run(solution: ConcreteInteractionDiagram):
    now = dt.datetime.now()

    solution.preanalysis()

    # Solve Biaxial
    nominalCapacity, reducedCapacity = solution.solve_biaxial_analysis()
    print(f"Biaxial analysis time:\t\t\t{dt.datetime.now() - now}")

    # Interpolate PM graphs at multiple thetas
    now = dt.datetime.now()
    target_thetas = [0, 90, 180, 270]
    Mns1, Pns1 = solution.interpolate_PM_at_thetas(target_thetas, nominalCapacity)
    phiMns1, phiPns1 = solution.interpolate_PM_at_thetas(target_thetas, reducedCapacity)

    # Interpolate PM graphs at multiple loads
    Mns2, Pns2, theta2, loads = solution.interpolate_PM_at_loads(solution.loads, nominalCapacity)
    phiMns2, phiPns2, phitheta2, philoads = solution.interpolate_PM_at_loads(loads, reducedCapacity)

    # Plot
    plot = Plot()
    plot.draw_section(solution)
    plot.plot3d(reducedCapacity)
    plot.plot3d(nominalCapacity)

    # Plot PM.
    # We plot reduced capacity first, then nominal capacity.
    # Larger values of nominal capacity will zoom out the final plot.
    plot.plot_PM_from_thetas(phiMns1, phiPns1, target_thetas, reduced=True)
    plot.plot_PM_from_thetas(Mns1, Pns1, target_thetas, reduced=False)

    # Plot PM
    plot.plot_PM_from_loads(phiMns2, phiPns2, philoads, phitheta2, reduced=True)
    plot.plot_PM_from_loads(Mns2, Pns2, loads, theta2, reduced=False)

    # Plot MM
    plot.plot_MM(nominalCapacity)
    print(f"Interpolation and plotting time:\t{dt.datetime.now() - now}")

    plot.show()
