import numpy as np
import pandas as pd


class Section:
    def __init__(self):
        pass

    def create_rectangle(self, width: float, height: float):
        return pd.DataFrame(
            {"x": [0, width, width, 0, 0], "y": [0, 0, height, height, 0]}
        )

    def create_circle(self, diameter: float, x: float = 0, y: float = 0, n: int = 30):
        theta = np.linspace(0, 2 * np.pi, n, endpoint=False)
        x = diameter / 2 * np.cos(theta) + x
        y = diameter / 2 * np.sin(theta) + y

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


class Reinforcement:
    def __init__(self, diameter: float):
        self.diameter = diameter

    def _distribute_rebars(self, p1, p2, n):
        if n < 2:
            return []
        xs = np.linspace(p1[0], p2[0], n)
        ys = np.linspace(p1[1], p2[1], n)
        return list(zip(xs, ys))

    def build_rectangle_rebar(
        self,
        width: float,
        height: float,
        n_left: int,
        n_top: int,
        n_right: int,
        n_bottom: int,
        cover: float = 40,
    ):
        p_bl = (cover, cover)
        p_br = (width - cover, cover)
        p_tr = (width - cover, height - cover)
        p_tl = (cover, height - cover)

        rebars = []
        rebars.extend(self._distribute_rebars(p_bl, p_br, n_bottom))
        rebars.extend(self._distribute_rebars(p_br, p_tr, n_right))
        rebars.extend(self._distribute_rebars(p_tr, p_tl, n_top))
        rebars.extend(self._distribute_rebars(p_tl, p_bl, n_left))

        df = pd.DataFrame(rebars, columns=["x", "y"])
        df = df.drop_duplicates().reset_index(drop=True)
        df["diameter"] = self.diameter
        return df

    def build_circle_rebar(self, diameter: float, n_rebar: int, cover: float = 40):
        radius = diameter / 2 - cover
        theta = np.linspace(0, 2 * np.pi, n_rebar, endpoint=False)
        x = radius * np.sin(theta)
        y = radius * np.cos(theta)
        return pd.DataFrame({"x": x, "y": y, "diameter": [self.diameter] * n_rebar})

    def build_lshape_rebar(
        self,
        shape_type: str,
        a: float,
        b: float,
        c: float,
        d: float,
        n_a: int,
        n_b: int,
        n_c: int,
        n_d: int,
        n_e: int,
        n_f: int,
        cover: float = 40,
    ):
        if shape_type == "L-shape-1":
            p = [
                (cover, cover),
                (a - cover, cover),
                (a - cover, d - cover),
                (c - cover, d - cover),
                (c - cover, b - cover),
                (cover, b - cover),
            ]
            n_list = [n_a, n_d, n_f, n_e, n_c, n_b]
        elif shape_type == "L-shape-2":
            p = [
                (cover, cover),
                (a - cover, cover),
                (a - cover, b - cover),
                (a - c + cover, b - cover),
                (a - c + cover, d - cover),
                (cover, d - cover),
            ]
            n_list = [n_a, n_b, n_c, n_e, n_f, n_d]
        elif shape_type == "L-shape-3":
            p = [
                (a - c + cover, cover),
                (a - cover, cover),
                (a - cover, b - cover),
                (cover, b - cover),
                (cover, b - d + cover),
                (a - c + cover, b - d + cover),
            ]
            n_list = [n_c, n_b, n_a, n_d, n_f, n_e]
        elif shape_type == "L-shape-4":
            p = [
                (cover, cover),
                (c - cover, cover),
                (c - cover, b - d + cover),
                (a - cover, b - d + cover),
                (a - cover, b - cover),
                (cover, b - cover),
            ]
            n_list = [n_c, n_e, n_f, n_d, n_a, n_b]
        else:
            return None

        rebars = []
        for i in range(6):
            rebars.extend(self._distribute_rebars(p[i], p[(i + 1) % 6], n_list[i]))

        df = pd.DataFrame(rebars, columns=["x", "y"])
        df = df.drop_duplicates().reset_index(drop=True)
        df["diameter"] = self.diameter
        return df
