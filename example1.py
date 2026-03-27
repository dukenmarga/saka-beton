import pandas as pd
from biaxial_analysis import ConcreteInteractionDiagram, run

def example1():
    return ConcreteInteractionDiagram(
        concrete_boundary = pd.DataFrame({
            'x': [0, 0, 250, 250],
            'y': [0, 400, 400, 0]
        }),
        rebar = pd.DataFrame({
            'diameter': [16, 16, 16, 22],
            'x': [50, 200, 50, 200],
            'y': [50, 50, 350, 350]
        }),
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
    run(example1())