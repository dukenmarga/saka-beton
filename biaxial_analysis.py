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
from functools import cached_property

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

    def check_clear_concrete_cover(self):
        # Group concrete boundary into edges
        concrete_edges: list[LineString] = []
        for i, p1 in enumerate(self.concrete_boundary[['x', 'y']].values):
            j = (i + 1) % len(self._concrete_boundary)
            p2 = self.concrete_boundary.iloc[j].values
            edge = LineString([(p1[0], p1[1]), (p2[0], p2[1])])
            concrete_edges.append(edge)

        # Check if any rebar is too close to concrete boundary
        for rebar in self._rebar.itertuples():
            p = Point(rebar.x, rebar.y)
            for edge in concrete_edges:
                # Calculate clear concrete cover from concrete boundary and outer side of rebar
                clear_cover = shapely.distance(edge, p)
                if (clear_cover - rebar.diameter / 2) < 20:
                    import sys
                    print("""WARNING: One or more rebar is too close to concrete boundary (clear cover less than 20mm).
                          The solution might not converge.""",
                        file=sys.stderr, flush=True)
                    return

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
        XO, YO = self.compute_centroid(self.concrete_boundary[['x', 'y']].values)

        # Transform concrete and rebar coordinates relative to
        # the new origin
        self.transform_to_new_origin(XO, YO)

        # Close loop
        self.close_boundary_loop()

        # Transform concrete and rebar coordinates relative to
        # the origin
        self.transform_to_O()

        # Calculate beta1
        self.beta1 = self.sni.compute_beta1(fc=self.material.fc_)

    @property
    def concrete_boundary(self):
        return self._concrete_boundary

    @property
    def rebar(self):
        return self._rebar

    def _brent_method(self, Pn_target: float, alpha: float = 0,
                      c_min: float = -100, c_max: float = 100,
                      tol: float = 1e-1, max_iter: float = 100):
        """
        Find combination of Mnx and Mny, such that Pn is close to Pn_target
            by utilizing Brent method.
        """
        self.Pn_target = Pn_target
        d, r = brentq(
            f=self.compute_Pn_brent,
            a=c_min,
            b=c_max,
            args=(alpha, Pn_target),
            full_output=True,
            disp=False,
            xtol=tol,
            rtol=tol,
            maxiter=max_iter
        )
        c = float(d)
        return self.compute_Pn_Mn(c, alpha)

    def compute_centroid(self, points: np.ndarray):
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

    def calculate_strain(self, rebar_y: np.ndarray, c: float):
        """
        Calculate strain for each rebar.
            rebar_y: y coordinate of rebar (local coordinate system)
            c: depth of compressed concrete from neutral axis
            epsilon_u: ultimate strain of concrete (usually 0.003)
        """
        # c zero means full compression
        if c == 0:
            return -self.material.epsilon_u

        return rebar_y / c * self.material.epsilon_u

    def calculate_strain_full_compression(self):
        """
        Calculate strain in full compression.
        """
        return self.material.epsilon_u

    def close_boundary_loop(self):
        """
        Check if the boundary is closed.
        If not, close it
        """
        if not self.first_and_last_node_same(self.concrete_boundary[['x', 'y']].values):
            first_point = self.concrete_boundary.loc[[0]]
            self._concrete_boundary = pd.concat([self.concrete_boundary, first_point], axis=0, ignore_index=True)

    def construct_line(self, point_on_line: tuple[float], alpha: float) -> LineString:
        """
        Construct a line that passes through point_on_line and is rotated by alpha.
        It can be used as virtual axis or neutral axis.
        """
        # Localising
        concrete_boundary = self.concrete_boundary
        # Extract x,y
        x_na, y_na = point_on_line
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

        # Build a line that rotates to alpha that passes through the point_on_na
        line = LineString([(x[0] - delta_x, y[0] - delta_y), (x[1] + delta_x, y[1] + delta_y)])
        return line

    def find_lowest_left_node(self):
        """
        Find the lowest left node
        """
        min_x = min(self.concrete_boundary['x'])
        min_y = min(self.concrete_boundary['y'])
        
        return min_x, min_y

    def _find_Mnx_Mny_at_Pn(self, Pn_target, alpha, c_min, c_max, method="secant", tol=1e-1, max_iter=100):
        """
        Find Mnx and Mny such that Pn is close to Pn_target
        """
        if method == "brent":
            tol = min(1e-4, tol)
            return self._brent_method(Pn_target, alpha, c_min, c_max, tol, max_iter)
        if method == "secant":
            return self._secant_method(Pn_target, alpha, c_min, c_max, tol, max_iter)
        raise ValueError(f"Unknown method: {method}")

    def find_Mx_My_at_P(self, P: float, alpha: float):
        """
        TODO: Need to rename this function. It will be used to
        calculate the phi Mn and phi Pn as well.
        Find combination of Mx and My that gives Pn close to P and alpha
        """
        # Construct virtual line that passes through point (0, 0)
        va = self.construct_line(point_on_line=(0, 0), alpha=alpha)

        # Get maximum depth from virtual axis to the farthest point of compression and tension
        _, yec, _, _, _, yet, _, _ = self.calculate_depth(va)
        h_max = yec + (-yet)
        self.h_max = h_max

        # c_min and c_max are limit where we will iterate to find Mnx and Mny.
        # c is measure from top of compression zone.
        # c = 5 mm means the distance from top of compression zone to the neutral axis.
        c_min = 0
        # c_max is the maximum possibility of neutral axis distance from the extreme of compression fibre.
        # Theoritically c_max is unlimited if we are at full compression.
        # But, there is a point where the concrete is assumed to be fully compressed.
        # The starting point is where c_max is equal to h_max/beta1. Beyond this point, we assume
        # the concrete is fully compressed.
        # So c_max theoritically can be h_max/beta1.
        # We add 1% to make sure we are not at the boundary.
        c_max = h_max / self.beta1 * 1.01

        Pn, Mnx, Mny, eps_t = self._find_Mnx_Mny_at_Pn(P, alpha, c_min, c_max)

        # Calculate theta, the angle between Mnx and Mny
        theta = math.atan2(Mny, Mnx) * 180 / math.pi

        # Normalize theta
        theta = self.normalize_theta(theta)

        # eps_t = [0] * len(alphas)
        eps_y = self.material.fy / self.material.Es
        sni = self.sni
        phi = sni.compute_phi(reduction_type="moment_axial", eps_t=eps_t, eps_y=eps_y)

        return Pn, Mnx, Mny, eps_t, theta, phi

    def compute_Pn_brent(self, c: float = 10, alpha: float = 0, Pn_target: float = 0):
        Pn = self.compute_Pn(c, alpha)
        return Pn - Pn_target

    def compute_Pn(self, c: float = 10, alpha: float = 0):
        """
        Calculate Pn from alpha and h
            alpha is in degrees
            c is distance from the most compressed point to the neutral axis
        """
        XO = self.XO
        YO = self.YO

        # Construct compression region.
        # Compression region is same as Whitney stress block.
        # In rectangle, this area is equal to Area = c * beta1 * bw
        compression_parts, na_line = self.construct_concrete_compression_region(c, alpha)

        # Calculate Pn for each compressed part
        Pcn = 0
        for compression_parts in compression_parts:
            Pcn_part = self.Pcn(compression_parts)
            Pcn += Pcn_part

        # Calculate coordinate O' (Xa, Ya) in global axis
        mid_na_line = na_line.interpolate(0.5, True)
        Xa = mid_na_line.x
        Ya = mid_na_line.y

        # Remap rebar coordinates in the rotated coordinate system (Neutral Axis) using new O'.
        rebar = copy.copy(self.rebar)
        rebar_local = geometry.map_to_local(rebar[['x', 'y']].values, XO, YO, Xa, Ya, na_line)
        rebar_local_pd = pd.DataFrame(rebar_local, columns=['x', 'y'])

        # Get strain (epsilon_s) for each point on the rebar, linearly distributed from Xe, Ye.
        # We use c, since it is the location of the neutral axis.
        # For full compression, we return the ultimate concrete strain
        if c >= self.h_max/self.beta1:
            rebar_strain_local = self.calculate_strain_full_compression()
        else:
            rebar_strain_local = self.calculate_strain(rebar_local_pd[['y']].values, c)
        # Concat strain (relative to Neutral Axis line)
        rebar['strain'] = rebar_strain_local

        # Calculate area of each rebar
        rebar_area_As = geometry.calculate_rebar_area(rebar[['diameter']])
        # Concat rebar area
        rebar['As'] = rebar_area_As

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
        # Fci = 0.85*fci * As
        # fci * As (MPa * mm2 / 1000 -> KN)
        rebar_As_fci = 0.85*rebar[['fci']].values * rebar[['As']].values / 1000
        rebar['As_fci'] = rebar_As_fci

        #############
        # SUMMATION #
        #############

        # Total nominal concrete capacity
        # As * fci -> KN
        sum_As_fci = np.sum(rebar[['As_fci']])

        # Total nominal rebar capacity
        # As * fs -> KN
        sum_As_fs = np.sum(rebar[['As_fs']])

        # Total nominal capacity
        # Pn: Nominal axial capacity
        # (to be consistent with formula from SNI/ACI subsection 22.4.2.2)
        # We introduce 0.001 KN to avoid error when using Brent method,
        # where f(b) might be very close to 0 but in negative side.
        # So, to prevent both f(a) and f(b) have same negative sign,
        # adding 0.001 KN is necessary to move f(b) to positive side.
        Pn = (Pcn - sum_As_fci) + sum_As_fs + 0.001

        return Pn

    def compute_Pn_Mn(self, c: float = 10, alpha: float = 0):
        """
        Calculate Pn, Mnx, Mny from alpha and h
            alpha is in degrees
            c is distance from the most compressed point to the neutral axis
        """
        # Initialize
        XO = self.XO
        YO = self.YO
        Pcn = 0
        Mcnx = 0
        Mcny = 0

        # Construct compression region.
        # Compression region is same as Whitney stress block.
        # In rectangle, this area is equal to Area = c * beta1 * bw
        compressed_parts, na_line = self.construct_concrete_compression_region(c, alpha)

        # Calculate the centroid of Whitney's compression stress block.
        ecc = self.compute_centroid_from_polygons(compressed_parts)

        for compressed_part in compressed_parts:
            Pcn_part = self.Pcn(compressed_part)
            Pcn += Pcn_part

            # Extract moment length
            x_depth = ecc[0]
            y_depth = ecc[1]

            # Moment capacity is always calculated using beta1
            Mcnx += Pcn_part * y_depth / 1e3
            Mcny += Pcn_part * x_depth / 1e3

        # Calculate coordinate O' (Xa, Ya) in global axis
        mid_na_line = na_line.interpolate(0.5, True)
        Xa = mid_na_line.x
        Ya = mid_na_line.y

        # Remap rebar coordinates in the rotated coordinate system (Neutral Axis) using new O'.
        rebar = copy.copy(self.rebar)
        rebar_local = geometry.map_to_local(rebar[['x', 'y']].values, XO, YO, Xa, Ya, na_line)
        rebar_local_pd = pd.DataFrame(rebar_local, columns=['x', 'y'])

        # Get strain (epsilon_s) for each point on the rebar, linearly distributed from Xe, Ye.
        # We use c, since it is the location of the neutral axis.
        # For full compression, we return the ultimate concrete strain
        if c >= self.h_max/self.beta1:
            rebar_strain_local = self.calculate_strain_full_compression()
        else:
            rebar_strain_local = self.calculate_strain(rebar_local_pd[['y']].values, c)
        # Concat strain (relative to Neutral Axis line)
        rebar['strain'] = rebar_strain_local

        # Calculate area of each rebar
        rebar_area_As = geometry.calculate_rebar_area(rebar[['diameter']])
        # Concat rebar area
        rebar['As'] = rebar_area_As

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
        # Fci = 0.85*fci * As
        # fci * As (MPa * mm2 / 1000 -> KN)
        rebar_As_fci = 0.85*rebar[['fci']].values * rebar[['As']].values / 1000
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
        # (to be consistent with formula from SNI/ACI subsection 22.4.2.2)
        # 0.85 is fc' reduction due to size effect based on Hognestad research
        Pn = (Pcn - sum_As_fci) + sum_As_fs
        # Mnx: Nominal moment capacity about global X-axis
        MnX = (Mcnx - sum_As_fci_Ysi) + sum_As_fs_Ysi
        # Mny: Nominal moment capacity about global Y-axis
        MnY = (Mcny - sum_As_fci_Xsi) + sum_As_fs_Xsi

        # Get strain of the most tensile rebar.
        # Absolute value of strain.
        eps_t = abs(min(rebar[['strain']].values)[0])

        # Return nominal capacity based on selected alpha and h,
        # and strain of the most tensile rebar
        return Pn, MnX, MnY, eps_t

    def compute_centroid_from_polygons(self, polygons: np.ndarray):
        """
        Compute centroid of polygons.

        The formula:
        x_cen: sum(x_i * A_i) / sum(A_i)
        y_cen: sum(y_i * A_i) / sum(A_i)
        """
        sum_x = 0
        sum_y = 0
        sum_area = 0
        for compressed_part in polygons:
            x_cen_part, y_cen_part = self.compute_centroid(compressed_part)
            area = geometry.calculate_concrete_area(compressed_part)
            sum_x += x_cen_part * area
            sum_y += y_cen_part * area
            sum_area += area

        # Avoid division by 0
        if sum_area == 0:
            return 0, 0

        x_cen = sum_x / sum_area
        y_cen = sum_y / sum_area

        centroids = [x_cen, y_cen]
        return centroids

    def construct_concrete_compression_region(self, c: float, alpha: float = 0):
        """
        Construct concrete compression region.

        There are 2 cases:
        1. c < h_max: the neutral axis is inside the compression region
        2. c >= h_max: the neutral axis is outside the compression region
        But, for the second case, we limit the compression region to h_max.
        Hence, compressed region is not always aligned with the neutral axis.
        We will distinguish these two lines as: compression line and neutral axis line.
        - c < h_max: neutral axis is aligned with compression line
        - c >= h_max: neutral axis is not aligned with compression line, compression line
          can't go lower than h_max, but neutral axis can go lower than h_max
        """
        XO = self.XO
        YO = self.YO

        # Construct benchmark line from a (0,0) that rotates to alpha
        helper_line = self.construct_line(point_on_line=(XO, YO), alpha=alpha)

        # Calculate:
        # xec: x coordinate of the most compressed point (local axis)
        # yec: y coordinate of the most compressed point (local axis)
        # yec = compression depth (same as c)
        # Xec: x coordinate of the most compressed point (global axis)
        # Yec: y coordinate of the most compressed point (global axis)
        xec, yec, Xec, Yec, _, yet, _, _ = self.calculate_depth(na_line=helper_line)

        # 1. NEUTRAL AXIS LINE

        na_point = np.array([[xec, yec-c]])

        # Map from local to global.
        # The NA point uses global coordinates.
        na_point = geometry.map_to_global(na_point, XO, YO, XO, YO, helper_line)

        # Construct neutral axis (measured from from xe, ye)
        na_line = self.construct_line(point_on_line=(na_point[0][0], na_point[0][1]), alpha=alpha)

        # 2. COMPRESSION LINE and COMPRESSION REGION

        h_comp = c*self.beta1

        # This is a point measured from the most compressed point
        # It still uses local coordinates.
        comp_point = np.array([[xec, yec-h_comp]])

        # Map from local to global.
        # The NA point uses global coordinates.
        comp_point = geometry.map_to_global(comp_point, XO, YO, XO, YO, helper_line)

        # Construct compression axis (measured away from xe, ye)
        comp_line = self.construct_line(point_on_line=(comp_point[0][0], comp_point[0][1]), alpha=alpha)

        # Isolate concrete boundary and rebar locally to this function only
        # - Readability
        concrete_boundary = copy.copy(self.concrete_boundary)

        # Make a polygon from the concrete boundary
        concrete_polygon = Polygon(concrete_boundary[['x', 'y']].values)

        # If h_comp beyond the section height, we return the whole concrete shape
        if h_comp >= self.h_max:
            return [np.array(concrete_polygon.exterior.coords.xy).T], na_line

        # If compresison height is zero, return nothing
        if h_comp == 0:
            return [], na_line

        # Get compressed part from the concrete if all the section is crossed by
        # the compression axis
        # Xe, Ye: indicator of the most compressed point, so it knows which part
        # is in compressed
        compressed_parts = geometry.get_compressed_parts(concrete_polygon, comp_line, (Xec, Yec))

        return compressed_parts, na_line

    def first_and_last_node_same(self, points: np.ndarray = None):
        """
        Check if the first and last node are the same
        """
        first_point = points[0]
        last_point = points[-1]
        return (first_point == last_point).all()

    def _secant_method(self, Pn_target: float, alpha: float = 0,
                      c_min: float = -100, c_max: float = 100,
                      tol: float = 1e-1, max_iter: float = 100):
        """
        Find combination of Mnx and Mny, such that Pn is close to Pn_target
            by utilizing secant method.
        """
        iter_count = 0
        c1 = c_min
        c2 = c_max
        c3 = 0
        self.Pn_target = Pn_target
        Pn1 = self.compute_Pn(c1, alpha)

        # Increase tolerance if Pn_target is close to Po.
        if Pn_target > 0.95*self.Po:
            tol = 3

        # Maximum number of iterations
        max_iter = 100

        # Initialize values to track best c
        best_diff = 10000
        best_c = 0

        # Iterate until convergence
        while iter_count < max_iter:
            iter_count += 1
            error = False

            try:
                Pn2 = self.compute_Pn(c2, alpha)
            except ValueError:
                error = True

            # Find the best c that will be used later
            # if it's not converging.
            # Usually this condition will often happen
            # near full compression.
            if abs(Pn2-Pn_target) < best_diff:
                best_diff = abs(Pn2-Pn_target)
                best_c = c2

            # Avoid division by zero
            if abs(Pn2 - Pn1) < 1e-6:
                error = True

            # If c1 and c2 are out of range (getting bigger for each iteration),
            # it indicates that the method is not converging.
            # The possibility is initial c1 and c2 are wrong.
            # Hence, we need to swap c1 and c2 (based on c_min and c_max).
            # We choose 5_000 as limit, since physically it is very unlikely
            # to have a concrete dimension greater than 5_000 mm
            elif abs(c1) > 5_000 or abs(c2) > 5_000:
                error = True

            # If there is an error, we are attempting to swap c1 and c2
            if error:
                # Reset and swap values (opposite of initial value before looping)
                c1 = c_max
                c2 = c_min
                Pn1 = self.compute_Pn(c1, alpha)
                continue

            # Update the new value by using secant formula:
            # h_new = c2 - [f(c2) - target] * (c2 - c1) / (f(c2) - f(c1))
            c3 = c2 - (Pn2 - Pn_target) * (c2 - c1) / (Pn2 - Pn1)

            # Return the result if the iteration has reached the convergence
            if abs(Pn2 - Pn_target) < tol:
                Pn2, Mnx2, Mny2, eps_t = self.compute_Pn_Mn(c2, alpha)
                return Pn2, Mnx2, Mny2, eps_t

            # To prevent c3 from being out of range (less than 0, impossible),
            # bisection method is being used by averaging c_min=0 and c2
            if c3 < 0:
                c3 = c2/2

            # To prevent c3 from being out of range (greater than c_max),
            # bisection method is being used by averaging c_max and c2
            if c3 > c_max:
                c3 = (c_max+c2)/2


            # Update the values for next iteration if there is no error
            # and still not reaching the convergence
            c1, Pn1 = c2, Pn2
            c2 = c3

        # Return the result from the best estimate if not convergence
        return self.compute_Pn_Mn(best_c, alpha)

    def interpolate_PM(self, theta_target: float, results: Capacity, element: list = ["Pn", "Mnx", "Mny"]):
        """
        Build the 2D PM interaction curve for a given bending direction theta.

        Background
        ----------
        The biaxial analysis produces a 3D interaction surface stored as a stack
        of rings — one ring per axial-load level Pn. Each ring is a closed polygon
        of moment capacity points (Mnx, Mny).

        To draw a 2D PM diagram the engineer picks a bending direction theta (e.g.
        0° for pure x-bending, 90° for pure y-bending). The question for each ring
        is: "what is the maximum moment the section can carry in direction theta at
        this axial load?" That moment equals the distance from the origin (0, 0) to
        the point on the ring that lies exactly along the theta direction.

        How the ring is cut
        -------------------
        Imagine standing at the origin and shining a flashlight in direction theta.
        The beam is a ray:  Mnx = M · cos(theta),  Mny = M · sin(theta).
        Where that beam hits the ring polygon is the PM boundary.

        To find the crossing numerically, each ring point is assigned a "side value":

            side = −sin(theta) × Mnx + cos(theta) × Mny

        side = 0 means the point is exactly on the theta ray (or its opposite).
        side > 0 means the point is to the left of the beam.
        side < 0 means the point is to the right of the beam.

        When two consecutive ring points are on opposite sides (one positive, one
        negative), the edge between them crosses the beam. Linear interpolation gives
        the exact crossing location (Mnx_cross, Mny_cross), and then:

            M = cos(theta) × Mnx_cross + sin(theta) × Mny_cross

        gives the moment at that crossing — which is the distance from the origin to
        that point on the ring.

        We keep the largest M (positive side of the PM curve) and smallest M
        (negative side), giving one point on each side per axial-load level.

        Special cases
        -------------
        • Full-tension ring: at the very bottom (full tension, no moment), all ring
          points collapse to nearly one location. The "side" values are all nearly
          zero, so there are no meaningful crossings. The single point is projected
          directly onto the theta direction instead.

        • Ring with no crossing on one side: the section has no capacity in that
          direction at this axial load. That side is set to zero.

        • Asymmetric sections near full tension: for an asymmetric section the
          full-tension point has an inherent moment eccentricity (e.g. the resultant
          lands at Mnx = −28 kNm, not at the origin). For a bending direction theta
          that does not pass through that point (e.g. theta = 90°), the capacity ring
          near full tension never crosses the theta beam. Those rings are excluded from
          the curve entirely. The two arcs are then closed with a V-tip at Mn = 0
          anchored to the Pn of the last excluded ring — this is the deepest axial
          load at which the section has any capacity in the theta direction.

        This works for any section shape — symmetric, asymmetric, L, T, or H.

        Parameters
        ----------
        theta_target : bending direction in degrees (0° = x-axis, 90° = y-axis)
        results      : Capacity object containing the precomputed capacity rings
        element      : three dictionary keys for [axial, x-moment, y-moment];
                       defaults to ["Pn", "Mnx", "Mny"] for nominal capacity and
                       should be changed to ["Pnr", "Mnxr", "Mnyr"] for reduced

        Returns
        -------
        Mn : list[float]  moment along the PM curve [kN·m]
        Pn : list[float]  corresponding axial load [kN]
        """
        Pn_string  = element[0]
        Mnx_string = element[1]
        Mny_string = element[2]

        capacity = results.capacity

        # Remove rings with duplicate Pn values (SNI axial-load cap).
        seen_pn: set = set()
        capacity_dedup = []
        for level in capacity:
            pn_key = round(level[Pn_string][0], 6)
            if pn_key not in seen_pn:
                seen_pn.add(pn_key)
                capacity_dedup.append(level)

        theta_rad = theta_target / 180 * math.pi
        cos_t     = math.cos(theta_rad)
        sin_t     = math.sin(theta_rad)

        def _ring_cut(Mnx_list, Mny_list):
            """
            Find the PM boundary of a single ring by cutting it with the theta beam.
            Returns (pos_M, neg_M): moment on the positive and negative side.
            """
            n = len(Mnx_list)

            # "Side value" for each ring point: positive = left of beam, negative = right.
            # A sign change between two consecutive points means the edge crosses the beam.
            side = [-sin_t * Mnx_list[i] + cos_t * Mny_list[i] for i in range(n)]

            # Special case: full-tension ring where all points are at nearly the same
            # location. All side values are nearly zero, so treat as a single point.
            ON_BEAM_TOLERANCE = 0.001  # kN·m
            if all(abs(s) < ON_BEAM_TOLERANCE for s in side):
                proj = [cos_t * Mnx_list[i] + sin_t * Mny_list[i] for i in range(n)]
                return max(proj), min(proj)

            pos_M = None  # moment on the +theta side (beam going forward)
            neg_M = None  # moment on the -theta side (beam going backward)

            for i in range(n):
                j      = (i + 1) % n   # next ring point (wraps around to 0 at the end)
                side_i = side[i]
                side_j = side[j]

                if side_i * side_j < 0:
                    # Edge crosses the beam: interpolate to find the exact crossing point.
                    # t = 0 is point i, t = 1 is point j; t = side_i/(side_i - side_j)
                    # places us exactly where side = 0.
                    t      = side_i / (side_i - side_j)
                    Mnx_cross = Mnx_list[i] + t * (Mnx_list[j] - Mnx_list[i])
                    Mny_cross = Mny_list[i] + t * (Mny_list[j] - Mny_list[i])
                elif abs(side_i) < ON_BEAM_TOLERANCE and abs(side_j) >= ON_BEAM_TOLERANCE:
                    # Point i sits exactly on the beam — use it directly as the crossing.
                    Mnx_cross = Mnx_list[i]
                    Mny_cross = Mny_list[i]
                else:
                    continue  # edge does not cross the beam

                # Moment at the crossing = distance from origin along the theta direction.
                M = cos_t * Mnx_cross + sin_t * Mny_cross
                if M >= 0:
                    pos_M = M if pos_M is None else max(pos_M, M)
                else:
                    neg_M = M if neg_M is None else min(neg_M, M)

            # Handle rings where one or both sides had no crossing.
            if pos_M is None or neg_M is None:
                # Check whether the ring is degenerate (all points close together).
                # This catches the full-tension ring that wasn't caught by the tolerance
                # check above (e.g. points are close but not all on the beam).
                cx    = sum(Mnx_list) / n
                cy    = sum(Mny_list) / n
                max_spread = max(math.sqrt((Mnx_list[i] - cx) ** 2 +
                                           (Mny_list[i] - cy) ** 2) for i in range(n))
                if max_spread < 1.0:
                    # Degenerate ring (e.g. full-tension ring collapsed to one point).
                    # Project the centroid onto theta. If the projection is effectively
                    # zero, the full-tension point doesn't lie in the theta direction —
                    # return None so the ring is excluded from the curve.
                    M_pt = cx * cos_t + cy * sin_t
                    if abs(M_pt) > ON_BEAM_TOLERANCE:
                        if pos_M is None:
                            pos_M = M_pt
                        if neg_M is None:
                            neg_M = M_pt
                    else:
                        return None, None
                else:
                    if pos_M is None and neg_M is None:
                        # Ring has no capacity in the theta direction at all.
                        # Exclude it so the curve starts at the first ring that does.
                        return None, None
                    # Only one side is missing: the boundary genuinely touches zero.
                    if pos_M is None:
                        pos_M = 0.0
                    if neg_M is None:
                        neg_M = 0.0

            return pos_M, neg_M

        max_proj_list: list = []
        min_proj_list: list = []
        pn_level_list: list = []

        first_included_idx = None   # ring index of the first ring with capacity in theta direction
        last_skipped_pn   = None    # Pn of the deepest ring that had no capacity (used for the V-tip)

        for idx, ring in enumerate(capacity_dedup):
            pn_level = float(ring[Pn_string][0])
            pos_M, neg_M = _ring_cut(
                [float(v) for v in ring[Mnx_string]],
                [float(v) for v in ring[Mny_string]],
            )
            if pos_M is None and neg_M is None:
                # Only track the tip candidate before the first included ring.
                # Skipped rings at the top (full compression) must be ignored.
                if first_included_idx is None:
                    last_skipped_pn = pn_level
                continue
            if first_included_idx is None:
                first_included_idx = idx
            max_proj_list.append(pos_M if pos_M is not None else 0.0)
            min_proj_list.append(neg_M if neg_M is not None else 0.0)
            pn_level_list.append(pn_level)

        # If rings at the bottom were skipped (section has no capacity in this
        # direction near full tension), insert a tip point at Mn=0 at the Pn of
        # the last skipped ring. Both arcs converge here, giving a sharp V bottom
        # instead of a flat horizontal cap.
        if last_skipped_pn is not None and first_included_idx is not None and first_included_idx > 0:
            max_proj_list.insert(0, 0.0)
            min_proj_list.insert(0, 0.0)
            pn_level_list.insert(0, last_skipped_pn)

        # Build closed PM curve: positive arc (ascending Pn) then negative arc
        # (descending Pn), closed back to the starting point.
        Mn = max_proj_list + list(reversed(min_proj_list)) + [max_proj_list[0]]
        Pn = pn_level_list + list(reversed(pn_level_list)) + [pn_level_list[0]]
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
            # Normalize theta
            theta = self.normalize_theta(theta)

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

    @cached_property
    def Po(self):
        """
        Nominal compressive strength, without reduction to consider accidental eccentricity
        """
        fc_ = self.material.fc_
        fy = self.material.fy
        concrete_area = geometry.calculate_concrete_area(self.concrete_boundary[['x', 'y']].values)
        rebar_area = np.sum(geometry.calculate_rebar_area(self.rebar['diameter'].values))

        return (0.85*fc_*(concrete_area-rebar_area) + fy*rebar_area) / 1e3

    def Pcn(self, concrete_compressed_boundary_xy: np.ndarray):
        """
        Calculate concrete nominal compression axial capacity, Pn.
        The result has not considered beta1 yet.

        concrete_boundary_xy: coordinates of concrete boundary (global axis)
        XO, YO: origin to calculate moment length
        """
        Pcn: float = 0

        # Calculate compression area
        compression_area = geometry.calculate_concrete_area(concrete_compressed_boundary_xy)

        # Calculate Pn contributed by concrete compression
        Pcn = 0.85*self.material.fc_ * compression_area / 1e3

        return Pcn

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
        Po = self.Po
        Pn_min = self.Pt_max()

        # Nominal capacity and reduced capacity
        nc = Capacity(False, [])
        rc = Capacity(True, [])

        # 25 evenly spaced levels from 1 KN above full tension to Po full compression.
        number = 25
        Pns_target = list(np.linspace(Pn_min + 1, Po, number))

        # Before calculating the solution in this zone, we need to determine
        # the maximum allowable axial capacity.
        phi_c = sni.compute_phi(reduction_type="compression", eps_t=0, eps_y=self.material.fy / self.material.Es)
        phiPo = phi_c * Po
        Pn_max = sni.compute_Pn_max(Po)
        phiPn_max = phi_c * Pn_max

        for _, Pn_target in enumerate(Pns_target):
            Pn0 = []
            Mnx0 = []
            Mny0 = []
            theta0 = []
            phi0 = []
            eps_t0 = []

            with multiprocessing.Pool(processes=16) as pool:
                # Map the function to the inputs and get results in order
                for Pn, Mnx, Mny, eps_t, theta, phi in pool.starmap(self.find_Mx_My_at_P, [(Pn_target, alpha) for alpha in alphas]):
                    Pn0.append(Pn)
                    Mnx0.append(Mnx)
                    Mny0.append(Mny)
                    eps_t0.append(eps_t)
                    theta0.append(theta)
                    phi0.append(phi)

            nc.capacity.append({
                "Pn": Pn0,
                "Mnx": Mnx0,
                "Mny": Mny0,
                "theta": theta0,
            })
            rc.capacity.append({
                "Pn": [sni.compute_phiPn(phi, Pn, phiPn_max) for phi, Pn in zip(phi0, Pn0)],
                "Mnx": [phi * Mnx for phi, Mnx in zip(phi0,  Mnx0)],
                "Mny": [phi * Mny for phi, Mny in zip(phi0, Mny0)],
                "theta": theta0,
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
    target_thetas = [0, 90]
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
