import numpy as np

class Material:
    def __init__(self, fc_, fy, fc__, fsu, epsilon_0, epsilon_u, gamma, Es, eps_sh, eps_su):
        self.fc_ = fc_
        self.fy = fy
        self.fc__ = 0.85 * fc_ if fc__ is None else fc__
        self.fsu = 1.1 * fy if fsu is None else fsu
        self.epsilon_0 = epsilon_0
        self.epsilon_u = epsilon_u
        self.gamma = gamma
        self.Es = Es
        self.eps_sh = eps_sh
        self.eps_su = eps_su

    def concrete_stress_strain(self, eps_c: float):
        fc__ = self.fc__
        gamma = self.gamma
        eps_0 = self.epsilon_0
        eps_u = self.epsilon_u

        if eps_c < 0:
            return 0
        # Ascending curve
        if eps_c <= eps_0:
            A = 2 * eps_c / eps_0
            B = (eps_c / eps_0) ** 2
            return fc__ * (A - B)
        # Descending line
        A = eps_c - eps_0
        B = eps_u - eps_0
        return fc__ * (1 - gamma * (A / B))

    def concrete_stress_strain_ufunc(self, eps_c: float):
        return np.vectorize(self.concrete_stress_strain, otypes=[float])(eps_c)

    # Calculate stress-strain curve for steel rebar
    # Developed from paper:
    # Stress-strain model for grade 275 reinforcing steel with cyclic loading
    # K. J. Thompson and R. Park (1978)
    def steel_rebar_stress_strain(self, eps_s: float):
        Es = self.Es
        fy = self.fy
        fsu = self.fsu
        eps_sh = self.eps_sh
        eps_su = self.eps_su

        sign = 1
        eps_y = fy / Es
        if eps_s < 0:
            sign = -1
            eps_s = -eps_s

        # Region 1 (Elastic)
        # fs = Es * eps_s
        if eps_s <= eps_y:
            return Es * eps_s * sign
        # Region 2 (Plateau)
        # fs = fy
        if eps_s <= eps_sh:
            return fy * sign

        return fy * sign

        # Note: Strain hardening region is not included at the moment.
        # Region 3 (Strain hardening)
        if eps_s <= eps_su:
            q = eps_su - eps_sh
            A = (fsu / fy) * (30 * q + 1)**2 - 60 * q - 1
            B = 15 * q**2
            m = A / B
            R = m * (eps_s - eps_sh) + 2
            S = 60 * (eps_s - eps_sh) + 2
            T = (eps_s - eps_sh) * (60 - m)
            U = 2 * (30 * q + 1)**2
            return fy * (R / S + T / U) * sign
        # Region 4 (Failing)
        return fsu*sign

    def steel_rebar_stress_strain_ufunc(self, eps_s: float):
        return np.vectorize(self.steel_rebar_stress_strain, otypes=[float])(eps_s)
