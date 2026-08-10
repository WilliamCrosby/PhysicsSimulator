from orbitalsim.universal_constants import UniversalConstants
from orbitalsim.simulation import calculate_accelerations

class RungeKutta4Integrator:

    @staticmethod
    def step(positions, velocities, masses, dt):

        k1_pos = velocities
        k1_vel = calculate_accelerations(positions, masses)

        k2_pos = velocities + 0.5 * dt * k1_vel
        k2_vel = calculate_accelerations(positions + 0.5 * dt * k1_pos, masses)

        k3_pos = velocities + 0.5 * dt * k2_vel
        k3_vel = calculate_accelerations(positions + 0.5 * dt * k2_pos, masses)

        k4_pos = velocities + dt * k3_vel
        k4_vel = calculate_accelerations(positions + dt * k3_pos, masses)

        positions += dt / 6 * (k1_pos + 2 * k2_pos + 2 * k3_pos + k4_pos)
        velocities += dt / 6 * (k1_vel + 2 * k2_vel + 2 * k3_vel + k4_vel)

        return positions, velocities