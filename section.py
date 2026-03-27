import numpy as np
import pandas as pd


class Section:
    def __init__(self):
        pass

    def create_rectangle(self, width: float, height: float):
        return pd.DataFrame(
            {"x": [0, width, width, 0, 0], "y": [0, 0, height, height, 0]}
        )

    def create_circle(self, diameter: float):
        theta = np.linspace(0, 2 * np.pi, 30, endpoint=False)
        x = diameter / 2 * np.cos(theta)
        y = diameter / 2 * np.sin(theta)

        return pd.DataFrame({"x": x, "y": y})

    def create_l_shape(self, shape_type: str, a: float, b: float, c: float, d: float):
        if shape_type == "L-shape-1":
            x = [0, a, a, c, c, 0, 0]
            y = [0, 0, d, d, b, b, 0]
        elif shape_type == "L-shape-2":
            x = [0, a, a, a - c, a - c, 0, 0]
            y = [0, 0, b, b, d, d, 0]
        elif shape_type == "L-shape-3":
            x = [a - c, a, a, 0, 0, a - c, a - c]
            y = [0, 0, b, b, b - d, b - d, 0]
        elif shape_type == "L-shape-4":
            x = [0, c, c, a, a, 0, 0]
            y = [0, 0, b - d, b - d, b, b, 0]
        else:
            return None

        return pd.DataFrame({"x": x, "y": y})
