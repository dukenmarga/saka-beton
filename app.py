# app.py
from flask import Flask, render_template, request, send_from_directory
import subprocess
import os
from biaxial_analysis import ConcreteInteractionDiagram
import pandas as pd
import pickle
import uuid
import sys
from section import Section, Reinforcement
from version import __version__

app = Flask(__name__)


@app.route("/")
def index():
    license_content = ""
    if os.path.exists("LICENSE"):
        with open("LICENSE", "r") as f:
            license_content = f.read()
    return render_template("index.html", license_content=license_content, version=__version__)


@app.route("/images/<path:filename>")
def serve_image(filename):
    return send_from_directory("images", filename)


@app.route("/analyse")
def analyse():
    # Get data from input form
    section_shape = request.args.get("section", "rectangle")
    fc = float(request.args.get("fc", "30"))
    fy = float(request.args.get("fy", "400"))
    traverse = request.args.get("transversal", "ties")

    # Reinforcement
    rebar_dia = float(request.args.get("rebar_dia", "16"))
    cover = float(request.args.get("cover", "40"))
    reinf = Reinforcement(rebar_dia)

    section = Section()
    openings = []
    if section_shape == "rectangle":
        width = float(request.args.get("width", "300"))
        height = float(request.args.get("height", "500"))
        concrete_boundary = section.create_rectangle(width, height)
        section_info = {"shape": section_shape, "width": width, "height": height}
        # Reinforcement
        n_bottom = int(request.args.get("n_bottom", "3"))
        n_top = int(request.args.get("n_top", "3"))
        n_left = int(request.args.get("n_left", "3"))
        n_right = int(request.args.get("n_right", "3"))
        rebar = reinf.build_rectangle_rebar(
            width, height, n_left, n_top, n_right, n_bottom, cover=cover
        )
    elif section_shape == "circle":
        diameter = float(request.args.get("diameter", "400"))
        concrete_boundary = section.create_circle(diameter)
        section_info = {"shape": section_shape, "diameter": diameter}
        # Reinforcement
        n_total = int(request.args.get("n_total", "8"))
        rebar = reinf.build_circle_rebar(diameter, n_total, cover=cover)
    elif section_shape.startswith("L-shape"):
        a = float(request.args.get("dim_a", "500"))
        b = float(request.args.get("dim_b", "500"))
        c = float(request.args.get("dim_c", "150"))
        d = float(request.args.get("dim_d", "150"))
        concrete_boundary = section.create_l_shape(section_shape, a, b, c, d)
        section_info = {"shape": section_shape, "a": a, "b": b, "c": c, "d": d}
        # Reinforcement
        n_a = int(request.args.get("n_a", "3"))
        n_b = int(request.args.get("n_b", "3"))
        n_c = int(request.args.get("n_c", "3"))
        n_d = int(request.args.get("n_d", "2"))
        n_e = int(request.args.get("n_e", "2"))
        n_f = int(request.args.get("n_f", "2"))
        rebar = reinf.build_lshape_rebar(
            section_shape, a, b, c, d, n_a, n_b, n_c, n_d, n_e, n_f, cover=cover
        )
    elif section_shape.startswith("T-shape"):
        a = float(request.args.get("t_dim_a", "600"))
        b = float(request.args.get("t_dim_b", "600"))
        c = float(request.args.get("t_dim_c", "200"))
        d = float(request.args.get("t_dim_d", "200"))
        concrete_boundary = section.create_t_shape(section_shape, a, b, c, d)
        section_info = {"shape": section_shape, "a": a, "b": b, "c": c, "d": d}
        # Reinforcement
        n_a = int(request.args.get("tn_a", "5"))
        n_b = int(request.args.get("tn_b", "3"))
        n_c = int(request.args.get("tn_c", "3"))
        n_d = int(request.args.get("tn_d", "4"))
        n_e = int(request.args.get("tn_e", "3"))
        n_f = int(request.args.get("tn_f", "4"))
        n_g = int(request.args.get("tn_g", "3"))
        n_h = int(request.args.get("tn_h", "3"))
        rebar = reinf.build_tshape_rebar(
            section_shape, a, b, c, d,
            n_a, n_b, n_c, n_d, n_e, n_f, n_g, n_h, cover=cover
        )
    elif section_shape == "hollow-circle":
        outer_dia = float(request.args.get("outer_dia", "600"))
        inner_dia = float(request.args.get("inner_dia", "400"))
        concrete_boundary, openings = section.create_hollow_circle(outer_dia, inner_dia)
        section_info = {"shape": section_shape, "outer_diameter": outer_dia, "inner_diameter": inner_dia}
        # Reinforcement
        n_ring = int(request.args.get("n_ring", "12"))
        rebar = reinf.build_circle_rebar(outer_dia, n_ring, cover=cover)
    elif section_shape == "hollow-rectangle":
        box_width = float(request.args.get("box_width", "800"))
        box_height = float(request.args.get("box_height", "1000"))
        box_thickness = float(request.args.get("box_thickness", "200"))
        concrete_boundary, openings = section.create_hollow_rectangle(
            box_width, box_height, box_thickness
        )
        section_info = {"shape": section_shape, "width": box_width, "height": box_height, "thickness": box_thickness}
        # Reinforcement
        n_bottom = int(request.args.get("box_n_bottom", "3"))
        n_top = int(request.args.get("box_n_top", "3"))
        n_left = int(request.args.get("box_n_left", "4"))
        n_right = int(request.args.get("box_n_right", "4"))
        rebar = reinf.build_hollow_rectangle_rebar(
            box_width, box_height, box_thickness,
            n_left, n_top, n_right, n_bottom, cover=cover
        )
    else:
        raise ValueError(f"Unknown section shape: {section_shape}")

    section_info["cover"] = cover

    # Parse loads from textarea
    loads_text = request.args.get("loads", "")
    pu_list, mux_list, muy_list = [], [], []
    if loads_text.strip():
        for line in loads_text.strip().split("\n"):
            try:
                parts = [p.strip() for p in line.replace(",", " ").split() if p.strip()]
                if len(parts) >= 3:
                    pu_list.append(float(parts[0]))
                    mux_list.append(float(parts[1]))
                    muy_list.append(float(parts[2]))
            except ValueError:
                continue

    # Use default loads if none provided
    if not pu_list:
        pu_list, mux_list, muy_list = (
            [1500, 1100, -300],
            [120, 210, 200],
            [170, 150, 400],
        )

    solution = ConcreteInteractionDiagram(
        concrete_boundary=concrete_boundary,
        rebar=rebar,
        openings=openings,
        section=section_info,
        fc_=fc,
        epsilon_0=0.002,
        epsilon_u=0.0038,
        gamma=0.15,
        fy=fy,
        Es=200000,
        eps_sh=0.015,
        eps_su=0.05,
        transverse=traverse,
        loads={
            "Pu": pu_list,
            "Mux": mux_list,
            "Muy": muy_list,
        },
    )

    # Dump solution to file
    tmp_dir = "tmp"
    os.makedirs(tmp_dir, exist_ok=True)
    random_name = f"{tmp_dir}/{str(uuid.uuid4())}"
    with open(random_name, "wb") as f:
        pickle.dump(solution, f)

    try:
        # Run subprocess
        script_path = os.path.join(os.path.dirname(__file__), "executor.py")
        result = subprocess.run(
            [sys.executable, script_path, random_name],
            capture_output=True,
            text=True,
            timeout=1000,
        )

        if result.returncode == 0:
            return ""
        else:
            return f"<div style='color:red; border:1px solid red; padding:15px; border-radius:8px; text-align:center;'><h2>Error:</h2><pre style='text-align:left;'>{result.stderr}</pre></div>"
    except Exception as e:
        return f"<div style='color:red; border:1px solid red; padding:15px; border-radius:8px; text-align:center;'><h2>Terjadi kesalahan:</h2><p>{str(e)}</p></div>"
    finally:
        # Remove temporary file
        os.remove(random_name)



if __name__ == "__main__":
    app.run(debug=True)
