from typing import Literal
import numpy as np

class SNI_2847_2019:
    def __init__(self, transverse: Literal["ties", "spiral"] = "ties"):
        self.transverse = transverse

    def compute_phi(
        self,
        reduction_type: Literal[
            "compression",  # pure compression
            "moment_axial", # combination of moment and axial
            "tension",      # pure tension
        ],
        eps_t: float = 0,
        eps_y: float = 0,
    ) -> float:
        """Set phi, Strength reduction"""
        if (0.005 - eps_y) == 0:
            return 0
        
        transversal_reinforcement_type = self.transverse

        # Checking if force is axial, moment, or both combined
        if reduction_type in ["compression", "moment_axial"]:
            # Compression-controlled
            # spiral (ACI 318M-14 Table 21.2.2(a))
            if eps_t <= eps_y and transversal_reinforcement_type == "spiral":
                phi = 0.75
            # other (ACI 318M-14 Table 21.2.2(b))
            elif eps_t <= eps_y and transversal_reinforcement_type in ["other", "ties"]:
                phi = 0.65

            # Tension-controlled
            # spiral and other (ACI 318M-14 Table 21.2.2(e, f))
            elif eps_t >= 0.005 and transversal_reinforcement_type in [
                "spiral",
                "other",
                "ties",
            ]:
                phi = 0.9

            # Transition: eps_y < eps_t < 0.005
            # spiral (ACI 318M-14 Table 21.2.2(c))
            elif transversal_reinforcement_type == "spiral":
                phi = 0.75 + 0.15 * (eps_t - eps_y) / (0.005 - eps_y)
            # other (ACI 318M-14 Table 21.2.2(d))
            elif transversal_reinforcement_type in ["other", "ties"]:
                phi = 0.65 + 0.25 * (eps_t - eps_y) / (0.005 - eps_y)

            # as failover if transversal_reinforcement_type is not defined or wrong
            else:
                phi = 0.65

        elif reduction_type == "tension":
            phi = 0.9

        # failover if reduction_type is not defined
        else:
            print("reduction_type is not defined")
            phi = 0.65

        return phi

    def compute_phi_vector(
            self,
            reduction_type: Literal[
                "compression",
                "moment_axial",
                "tension",
            ],
            eps_t: float = 0,
            eps_y: float = 0,
        ):
        return np.vectorize(self.compute_phi)(reduction_type, eps_t, eps_y)


    def compute_phiPn(
        self, phi: float = 0, Pn: float = 0, phiPn_max: float = 99_999_999
    ) -> float:
        """Compute reduced compression axial capacity of column"""

        phiPn = Pn * phi

        if phiPn > phiPn_max:
            phiPn = phiPn_max
        return phiPn

    def compute_phiPn_max(
        self,
        phiPn: float = 0,
    ) -> float:
        """Compute pure compression axial capacity of column,
        but reduced to anticipate accidental eccentricity, KN.
        """

        # Compute constant for maximum axial capacity
        # See Table 22.4.21 ACI 318_14
        if self.transverse == "ties":
            K = 0.8
        elif self.transverse in ["spiral", "other"]:
            K = 0.85
        else:
            K = 0.0

        phiPn = K * phiPn

        return phiPn
