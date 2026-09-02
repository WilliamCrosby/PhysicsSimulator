from orbitalsim.simulation import calculate_accelerations
from orbitalsim.universal_constants import UniversalConstants

class SymplecticEulerIntegrator:

    @staticmethod
    def step(positions, velocities, GMs, dt):
        accelerations = calculate_accelerations(positions, GMs)
        velocities += accelerations * dt
        positions += velocities * dt
        return positions, velocities