import numpy as np
import shapely
import shapely.ops
from shapely.geometry import LineString, Point, Polygon
import pandas as pd

def map_to_local(nodes: np.ndarray, x_cen_glob: float, y_cen_glob: float, x_cen_loc: float, y_cen_loc: float, line: LineString):
    """Transform nodes coordinate from global to local,
    with x_cen_glob, y_cen_glob as the center of the global coordinate system
    and x_cen_loc, y_cen_loc as the center of the local coordinate system
    and line as the line of the local coordinate system
    """
    Xa = x_cen_loc - x_cen_glob
    Ya = y_cen_loc - y_cen_glob
    cos_alpha = cos_from_line(line)
    sin_alpha = sin_from_line(line)
    Xi = nodes[:, 0]
    Yi = nodes[:, 1]
    x_ = (Xi - Xa) * cos_alpha + (Ya - Yi) * sin_alpha
    y_ = (Yi - Ya) * cos_alpha + (Xi - Xa) * sin_alpha
    return np.column_stack((x_, y_))

def map_to_global(local_nodes: np.ndarray, x_cen_glob: float, y_cen_glob: float, x_cen_loc: float, y_cen_loc: float, line: LineString):
    """Transform nodes coordinate from local to global,
    with x_cen_glob, y_cen_glob as the center of the global coordinate system
    and x_cen_loc, y_cen_loc as the center of the local coordinate system
    and line as the line of the local coordinate system
    """
    Xa = x_cen_loc - x_cen_glob
    Ya = y_cen_loc - y_cen_glob
    x = local_nodes[:, 0]
    y = local_nodes[:, 1]
    cos_alpha = cos_from_line(line)
    sin_alpha = sin_from_line(line)
    Xi = Xa + x*cos_alpha + y*sin_alpha
    Yi = Ya - x*sin_alpha + y*cos_alpha
    return np.column_stack((Xi, Yi))

def sin_from_line(line: LineString):
    x1 = line.coords[0][0]
    y1 = line.coords[0][1]
    x2 = line.coords[1][0]
    y2 = line.coords[1][1]
    r = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
    return (y1 - y2) / r

def cos_from_line(line: LineString):
    x1 = line.coords[0][0]
    y1 = line.coords[0][1]
    x2 = line.coords[1][0]
    y2 = line.coords[1][1]
    r = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
    return (x2 - x1) / r


def find_farthest_node(nodes: np.ndarray):
    """Find the farthest node from the neutral axis (perpendicular to the neutral axis)
    """
    # sort nodes based on y (second column)
    sorting_indices = nodes[:, 1].argsort()
    sorted = nodes[sorting_indices]

    # return farthest nodes:
    # - the lowest (usually in tension or negative value) and
    # - the highest (usually in compression or positive value)
    return sorted[0], sorted[-1]


def calculate_rebar_area(rebar_dia: np.ndarray):
    """Calculate rebar area
    """
    return np.pi * rebar_dia ** 2 / 4

def calculate_concrete_area(concrete_boundary: np.ndarray):
    """Calculate concrete area
    """
    area = calculate_signed_area(concrete_boundary)
    return abs(area)

def calculate_signed_area(concrete_boundary: np.ndarray):
    """Calculate signed area
    """
    n = len(concrete_boundary)
    area = 0
    for i in range(n):
        j = (i + 1) % n
        area += concrete_boundary[i, 0] * concrete_boundary[j, 1]
        area -= concrete_boundary[j, 0] * concrete_boundary[i, 1]
    area = area / 2
    return area

def polygon_oritentation(concrete_boundary: np.ndarray):
    """
    Check if the polygon is clockwise or counterclockwise using signed area method
    """
    area = calculate_signed_area(concrete_boundary)
    if area < 0:
        return "clockwise"
    elif area > 0:
        return "counterclockwise"
    else:
        raise ValueError("Polygon is degenerate")

def make_clockwise(concrete_boundary: pd.DataFrame):
    concrete_boundary_np = concrete_boundary[['x', 'y']].to_numpy(copy=True)
    if polygon_oritentation(concrete_boundary_np) == "clockwise":
        return concrete_boundary
    else:
        return concrete_boundary[::-1].reset_index(drop=True)

def make_counterclockwise(concrete_boundary: pd.DataFrame):
    concrete_boundary_np = concrete_boundary[['x', 'y']].to_numpy(copy=True)
    if polygon_oritentation(concrete_boundary_np) == "counterclockwise":
        return concrete_boundary
    else:
        return concrete_boundary[::-1].reset_index(drop=True)

def get_compressed_parts(polygon: Polygon, line: LineString, XeYe: tuple[float, float]):
    """
    Get the compression part of the polygon by the line.
    XeYe is an indicator node that is on compression part.
    """
    Xe = XeYe[0]
    Ye = XeYe[1]
    # Cut the polygon into 2 polygons by the line
    # At this point, we don't know which part is the compression part
    parts = shapely.ops.split(polygon, line)

    # If no compression part is found, return only small square (1x1 unit)
    if len(parts.geoms) <= 1:
        return np.array([[[0, 0],
                         [0, 1],
                         [1, 1],
                         [1, 0]]])

    # compr_index holds the index of the polygon part
    # where XeYe is located 
    compr_index = None
    # side_map will holds the left parts and right parts from cutting line.
    # Later, we will determine which part is the compression part.
    side_map = {
        "left": [],
        "right": []
    }
    # We create virtual line by offsetting the cutting line (neutral axis)
    line_offset = line.offset_curve(-0.00001)
    # For each split part, check if the part is on the left or right side
    # by checking if the part intersects with the virtual line
    for i, geom in enumerate(parts.geoms):
        if geom.intersects(line_offset):
            side_map["right"].append(i)
        else:
            side_map["left"].append(i)

        # Check if this part contains node XeYe.
        # If yes, then this part is the compression part
        for _, node in enumerate(geom.exterior.coords):
            if abs(Xe - node[0]) < 0.01 and abs(Ye - node[1]) < 0.01:
                compr_index = i

    # Check compr_index is in side_map["left"] or side_map["right"]
    # Then, we will get all indexes of the compression part
    if compr_index in side_map["left"]:
        compression_index = side_map["left"]
    elif compr_index in side_map["right"]:
        compression_index = side_map["right"]
    else:
        print(f"Error: compr_index: {compr_index} not in side_map: {side_map}")

    # Get the coordinates of compression part
    # The format: x is the first row, y is the second row
    compression_parts = []
    for i in compression_index:
        compression_part = np.array(parts.geoms[i].exterior.coords.xy)
        # Transpose so x is the first column, y is the second
        compression_part = compression_part.T
        compression_parts.append(compression_part)

    return compression_parts

def combine(concrete_boundary: pd.DataFrame, openings: list[pd.DataFrame]):
    """
    Combine concrete boundary and openings as one polygon
    """
    if len(openings) == 0:
        return concrete_boundary

    combined: pd.DataFrame = pd.DataFrame(columns=['x', 'y'])
    smallest_distance = 9999999
    node_at_boundary = 0
    node_at_opening = 0

    # Connect each opening to the concrete boundary
    for i, opening in enumerate(openings):
        # Get the smallest distance between the opening and the concrete boundary
        for row in opening.itertuples():
            x = row.x
            y = row.y

            # Get the distance between the node and the concrete boundary
            distance = np.sqrt((concrete_boundary['x'] - x)**2 + (concrete_boundary['y'] - y)**2)
            min_distance = distance.min()
            min_distance_node_index = distance.argmin()

            # Update the minimum distance and index
            if min_distance < smallest_distance:
                smallest_distance = min_distance
                node_at_boundary = min_distance_node_index
                node_at_opening = row.Index

        # Shift opening nodes, so that the selected node at opening become the first node in the array
        # pd.concat([opening.iloc[node_at_opening:], opening.iloc[:node_at_opening]], axis=0, ignore_index=True)
        opening = pd.concat([opening.iloc[node_at_opening:], opening.iloc[:node_at_opening]]).reset_index(drop=True)

        # Get the node index at boundary
        node_boundary = int(node_at_boundary)

        # Determine if the increment should be positive or negative
        x = concrete_boundary['x'].iloc[node_boundary]
        y1 = concrete_boundary['y'].iloc[node_boundary]
        y2 = opening['y'].iloc[0]
        opening_y = opening['y'].iloc[0]
        increment = y1 > y2

        # Two additional lines are needed to connect concrete boundary and opening,
        # one is the line from opening, the other is the line from concrete boundary.
        # But, these 2 lines can' share the same node.
        # Hence, we will make a small gap between these 2 lines.
        dx = 0.1
        if increment:
            x = x + dx
            opening_x = opening['x'].iloc[0] + dx
        else:
            x = x - dx
            opening_x = opening['x'].iloc[0] - dx

        # Combine the concrete boundary and openings
        combined = pd.concat([concrete_boundary[:node_boundary+1], opening], axis=0, ignore_index=True)
        combined['x'] = combined['x'].astype('float64')
        combined = pd.concat([combined, pd.DataFrame({'x': [opening_x], 'y': [opening_y]})], axis=0, ignore_index=True)
        combined = pd.concat([combined, pd.DataFrame({'x': [x], 'y': [y1]})], axis=0, ignore_index=True)
        combined = pd.concat([combined, concrete_boundary[node_boundary+1:]], axis=0, ignore_index=True)

        # Update concrete boundary
        concrete_boundary = combined

        # Reset smallest distance
        smallest_distance = 9999999

    return combined