import numpy as np
import shapely
import shapely.ops
from shapely.geometry import LineString, Point, Polygon

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
    n = len(concrete_boundary)
    area = 0
    for i in range(n):
        j = (i + 1) % n
        area += concrete_boundary[i, 0] * concrete_boundary[j, 1]
        area -= concrete_boundary[j, 0] * concrete_boundary[i, 1]
    area = abs(area) / 2
    return area


def get_compression_part(polygon: Polygon, line: LineString, XeYe: tuple[float, float]):
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
    if len(parts.geoms) == 1:
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
    for _, i in enumerate(compression_index):
        compression_part = np.array(parts.geoms[i].exterior.coords.xy)
        # Transpose so x is the first column, y is the second
        compression_part = compression_part.T
        compression_parts.append(compression_part)

    return compression_parts
