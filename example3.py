import pandas as pd
from biaxial_analysis import ConcreteInteractionDiagram, run

def example3():
    return ConcreteInteractionDiagram(
        concrete_boundary = pd.DataFrame({
            'x': [150, 150, 300, 300, -300, -300, -150, -150, -300, -300, 300, 300],
            'y': [-200, 200, 200, 350, 350, 200, 200, -200, -200, -350, -350, -200]
        }),
        rebar = pd.DataFrame({
            'diameter': [16, 16, 16, 16, 16, 16, 16, 16],
            'x': [-250, -100, 100, 250, -250, -100, 100, 250],
            'y': [300, 300, 300, 300, -300, -300, -300, -300],
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
    run(example3())