import numpy as np

from orbitalsim.calculations import mechanical_energy_calculator, momentum_calculator, update_extrema
from orbitalsim.celestial_body import CelestialBody
from orbitalsim.universal_constants import UniversalConstants

def calculate_accelerations(positions, GMs):
    accelerations = np.zeros_like(positions)
    n = len(GMs)
    G = UniversalConstants.G

    for i in range(n):
        position_i = positions[i]
        acceleration_i = accelerations[i]

        for j in range(n):
            if i == j:
                continue

            r = positions[j] - position_i
            dist_sq = np.dot(r, r)

            acceleration_i += (GMs[j] * r / (dist_sq * np.sqrt(dist_sq)))

        accelerations[i] = acceleration_i

    return accelerations

class NBodySimulation:

    def __init__(self, integrator, bodies: list[CelestialBody]):
        self.integrator = integrator
        self.bodies = bodies

    def simulator(self, steps, dt, do_extrema_calculations=False):

        bodies = self.bodies
        integrator = self.integrator

        positions = np.array([b.position for b in bodies], dtype = float)
        velocities = np.array([b.velocity for b in bodies], dtype = float)
        GMs = np.array([b.GM for b in bodies], dtype = float)
        masses = np.array([b.mass for b in bodies], dtype = float)

        trajectories = np.empty((steps + 1, len(bodies), positions.shape[1]), dtype = float) # for visuals later
        trajectories[0] = positions

        saved_velocities = np.empty((steps + 1, len(bodies), velocities.shape[1]), dtype = float)
        saved_velocities[0] = velocities

        energy = mechanical_energy_calculator(masses, positions, velocities)
        maximum_energy = energy
        minimum_energy = energy

        lin_momentum, ang_momentum = momentum_calculator(masses, positions, velocities)
        maximum_linear_momentum = np.linalg.norm(lin_momentum)
        minimum_linear_momentum = np.linalg.norm(lin_momentum)
        maximum_angular_momentum = np.linalg.norm(ang_momentum)
        minimum_angular_momentum = np.linalg.norm(ang_momentum)


        # primary loop
        for step in range(steps):

            integrator.step(positions, velocities, GMs, dt)

            trajectories[step + 1] = positions
            saved_velocities[step + 1] = velocities

            if do_extrema_calculations:
                maximum_energy, minimum_energy, maximum_linear_momentum, minimum_linear_momentum, maximum_angular_momentum, minimum_angular_momentum = update_extrema(masses, positions, velocities, maximum_energy, minimum_energy, maximum_linear_momentum, minimum_linear_momentum, maximum_angular_momentum, minimum_angular_momentum)


        for i, body in enumerate(bodies):
            body.position = positions[i]
            body.velocity = velocities[i]

        if do_extrema_calculations:
            return trajectories, saved_velocities, positions, maximum_energy, minimum_energy, maximum_linear_momentum, minimum_linear_momentum, maximum_angular_momentum, minimum_angular_momentum
        else:
            return trajectories, saved_velocities, positions
