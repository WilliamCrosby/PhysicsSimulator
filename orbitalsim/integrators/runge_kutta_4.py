from orbitalsim.universal_constants import UniversalConstants
from orbitalsim.simulation import calculate_accelerations

class RungeKutta4Integrator:

    @staticmethod
    def step(bodies, positions, velocities):

        k1_pos = velocities
        k1_vel = calculate_accelerations(bodies, positions)

        k2_pos = velocities + 0.5 * UniversalConstants.dt * k1_vel
        k2_vel = calculate_accelerations(bodies, positions + 0.5 * UniversalConstants.dt * k1_pos)

        k3_pos = velocities + 0.5 * UniversalConstants.dt * k2_vel
        k3_vel = calculate_accelerations(bodies, positions + 0.5 * UniversalConstants.dt * k2_pos)

        k4_pos = velocities + UniversalConstants.dt * k3_vel
        k4_vel = calculate_accelerations(bodies, positions + UniversalConstants.dt * k3_pos)

        positions += UniversalConstants.dt / 6 * (k1_pos + 2 * k2_pos + 2 * k3_pos + k4_pos)
        velocities += UniversalConstants.dt / 6 * (k1_vel + 2 * k2_vel + 2 * k3_vel + k4_vel)

        return positions, velocities