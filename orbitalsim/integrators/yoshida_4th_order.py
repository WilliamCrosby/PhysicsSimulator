import numpy as np
from orbitalsim.simulation import calculate_accelerations
from orbitalsim.integrators.velocity_verlet import VelocityVerletIntegrator

class Yoshida4thOrderIntegrator:

    # currently designed to work with velocity verlet, can work with any symplectic time-reversible integrator (i think)
    @staticmethod
    def step(bodies, positions, velocities, dt):

        w1 = 1.0 / (2.0 - 2.0**(1.0/3.0))
        w0 = -2.0**(1.0/3.0) / (2.0 - 2.0**(1.0/3.0))

        positions, velocities = VelocityVerletIntegrator.step(bodies, positions, velocities, dt * w1)
        positions, velocities = VelocityVerletIntegrator.step(bodies, positions, velocities, dt * w0)
        positions, velocities = VelocityVerletIntegrator.step(bodies, positions, velocities, dt * w1)

        return positions, velocities