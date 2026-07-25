from orbitalsim.simulation import calculate_accelerations
from orbitalsim.universal_constants import UniversalConstants

class SymplecticEulerIntegrator:

    @staticmethod
    def step(bodies, positions, velocities, dt):
        accelerations = calculate_accelerations(bodies, positions)
        velocities += accelerations * dt
        positions += velocities * dt
        return positions, velocities