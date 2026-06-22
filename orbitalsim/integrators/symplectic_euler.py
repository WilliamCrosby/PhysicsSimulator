from orbitalsim.universal_constants import UniversalConstants


class SymplecticEulerIntegrator:

    @staticmethod
    def step(bodies, positions, velocities, accelerations):
        velocities += accelerations * UniversalConstants.dt
        positions += velocities * UniversalConstants.dt
        return positions, velocities