import pandas as pd
from biaxial_analysis import ConcreteInteractionDiagram, run
from section import Section, Reinforcement

def example4():
    section = Section()
    reinforcement = Reinforcement(16)
    return ConcreteInteractionDiagram(
        concrete_boundary = section.create_l_shape("L-shape-3", 500, 500, 150, 150),
        rebar = reinforcement.build_lshape_rebar(
            "L-shape-3", 500, 500, 150, 150, 6, 6, 2, 2, 5, 5, 40
        ),
        fc_ = 27.58,
        epsilon_0 = 0.002,
        epsilon_u = 0.0038,
        gamma = 0.15,
        fy = 413.69,
        Es = 200000,
        eps_sh = 0.015,
        eps_su = 0.05,
        transverse="ties",
        loads = {
            "Pu": [1500, 1100, -300],
            "Mux": [120, 210, 200],
            "Muy": [170, 150, 400],
        }
    )

if __name__ == "__main__":
    run(example4())