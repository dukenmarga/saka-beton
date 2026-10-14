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


    def create_t_shape(self, shape_type: str, a: float, b: float, c: float, d: float):
        """
        a: overall width, b: overall height,
        c: web thickness, d: flange thickness
        """
        e = (a - c) / 2
        f = (b - c) / 2
        if shape_type == "T-shape-1":
            x = [a, 0, 0, e, e, e + c, e + c, a, a]
            y = [b, b, b - d, b - d, 0, 0, b - d, b - d, b]
        elif shape_type == "T-shape-2":
            x = [a, a, a - d, a - d, 0, 0, a - d, a - d, a]
            y = [0, b, b, f + c, f + c, f, f, 0, 0]
        elif shape_type == "T-shape-3":
            x = [0, a, a, e + c, e + c, e, e, 0, 0]
            y = [0, 0, d, d, b, b, d, d, 0]
        elif shape_type == "T-shape-4":
            x = [0, 0, d, d, a, a, d, d, 0]
            y = [b, 0, 0, f, f, f + c, f + c, b, b]
        else:
            return None

        return pd.DataFrame({"x": x, "y": y})

    def create_hollow_circle(
        self, outer_diameter: float, inner_diameter: float, n: int = 30
    ):
        return self.create_circle(outer_diameter, n=n), [
            self.create_circle(inner_diameter, n=n)
        ]

    def create_hollow_rectangle(self, width: float, height: float, thickness: float):
        t = thickness
        # Both loops are left open (no repeated closing node), like create_circle.
        # A closed loop makes combine() leave a zero-length edge in the merged
        # boundary, which the mesher rejects.
        boundary = pd.DataFrame({"x": [0, width, width, 0], "y": [0, 0, height, height]})
        opening = pd.DataFrame(
            {
                "x": [t, width - t, width - t, t],
                "y": [t, t, height - t, height - t],
            }
        )
        return boundary, [opening]


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

    def build_tshape_rebar(
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
        n_g: int,
        n_h: int,
        cover: float = 40,
    ):
        e = (a - c) / 2
        f = (b - c) / 2
        if shape_type == "T-shape-1":
            p = [
                (a - cover, b - cover),
                (cover, b - cover),
                (cover, b - d + cover),
                (e + cover, b - d + cover),
                (e + cover, cover),
                (e + c - cover, cover),
                (e + c - cover, b - d + cover),
                (a - cover, b - d + cover),
            ]
        elif shape_type == "T-shape-2":
            p = [
                (a - cover, cover),
                (a - cover, b - cover),
                (a - d + cover, b - cover),
                (a - d + cover, f + c - cover),
                (cover, f + c - cover),
                (cover, f + cover),
                (a - d + cover, f + cover),
                (a - d + cover, cover),
            ]
        elif shape_type == "T-shape-3":
            p = [
                (cover, cover),
                (a - cover, cover),
                (a - cover, d - cover),
                (e + c - cover, d - cover),
                (e + c - cover, b - cover),
                (e + cover, b - cover),
                (e + cover, d - cover),
                (cover, d - cover),
            ]
        elif shape_type == "T-shape-4":
            p = [
                (cover, b - cover),
                (cover, cover),
                (d - cover, cover),
                (d - cover, f + cover),
                (a - cover, f + cover),
                (a - cover, f + c - cover),
                (d - cover, f + c - cover),
                (d - cover, b - cover),
            ]
        else:
            return None

        n_list = [n_a, n_b, n_c, n_d, n_e, n_f, n_g, n_h]

        rebars = []
        for i in range(8):
            rebars.extend(self._distribute_rebars(p[i], p[(i + 1) % 8], n_list[i]))

        df = pd.DataFrame(rebars, columns=["x", "y"])
        df = df.drop_duplicates().reset_index(drop=True)
        df["diameter"] = self.diameter
        return df

    def build_hollow_rectangle_rebar(
        self,
        width: float,
        height: float,
        thickness: float,
        n_left: int,
        n_top: int,
        n_right: int,
        n_bottom: int,
        cover: float = 40,
    ):
        outer = self.build_rectangle_rebar(
            width, height, n_left, n_top, n_right, n_bottom, cover=cover
        )
        # Inner layer sits at cover from the opening face
        inner = self.build_rectangle_rebar(
            width - 2 * thickness,
            height - 2 * thickness,
            n_left,
            n_top,
            n_right,
            n_bottom,
            cover=cover,
        )
        inner["x"] += thickness
        inner["y"] += thickness

        df = pd.concat([outer, inner], ignore_index=True)
        df = df.drop_duplicates().reset_index(drop=True)
        return df
