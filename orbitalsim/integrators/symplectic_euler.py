from orbitalsim.simulation import calculate_accelerations
from orbitalsim.universal_constants import UniversalConstants

class SymplecticEulerIntegrator:

    @staticmethod
    def step(positions, velocities, masses, dt):
        accelerations = calculate_accelerations(positions, masses)
        velocities += accelerations * dt
        positions += velocities * dt
        return positions, velocities