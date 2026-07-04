from orbitalsim.simulation import calculate_accelerations
from orbitalsim.universal_constants import UniversalConstants

class SymplecticEulerIntegrator:

    @staticmethod
    def step(bodies, positions, velocities):
        accelerations = calculate_accelerations(bodies, positions)
        velocities += accelerations * UniversalConstants.dt
        positions += velocities * UniversalConstants.dt
        return positions, velocities