from orbitalsim.diagnostic_calculators import *
from orbitalsim.universal_constants import *

def __init__(self, integrator, bodies):
    self.integrator = integrator
    self.bodies = bodies

def CalculateAccelerations(bodies, positions):
    accelerations = np.zeros_like(positions)

    for i in range(len(bodies)):
        for j in range(len(bodies)):
            if i == j:
                continue

            r = positions[j] - positions[i]

            dist_sq = np.dot(r, r)
            dist = np.sqrt(dist_sq)

            accelerations[i] += (UniversalConstants.G * bodies[j].mass * r / dist ** 3)

    return accelerations


class NBodySimulation:


    def __init__(self, integrator, bodies: list[CelestialBody]):
        self.integrator = integrator
        self.bodies = bodies

    def Simulator(self, steps):

        positions = np.array([b.position.copy() for b in self.bodies], dtype = float)
        velocities = np.array([b.velocity.copy() for b in self.bodies], dtype = float)

        trajectories = [positions.copy()] # for visuals later

        energy = MechanicalEnergyCalculator(self.bodies, positions, velocities)
        maximumEnergy = energy
        minimumEnergy = energy

        linMomentum, angMomentum = MomentumCalculator(self.bodies, positions, velocities)
        maximumLinearMomentum = np.linalg.norm(linMomentum)
        minimumLinearMomentum = np.linalg.norm(linMomentum)
        maximumAngularMomentum = np.linalg.norm(angMomentum)
        minimumAngularMomentum = np.linalg.norm(angMomentum)


        # primary loop
        for step in range(steps):

            accelerations = CalculateAccelerations(self.bodies, positions)

            self.integrator.Step(self.bodies, positions, velocities, accelerations)

            trajectories.append(positions.copy())

            maximumEnergy, minimumEnergy, maximumLinearMomentum, minimumLinearMomentum, maximumAngularMomentum, minimumAngularMomentum = UpdateExtrema(self.bodies, positions, velocities, maximumEnergy, minimumEnergy, maximumLinearMomentum, minimumLinearMomentum, maximumAngularMomentum, minimumAngularMomentum)


        for i, body in enumerate(self.bodies):
                body.positions = positions[i]
                body.velocities = velocities[i]

        return np.array(trajectories), maximumEnergy, minimumEnergy, maximumLinearMomentum, minimumLinearMomentum, maximumAngularMomentum, minimumAngularMomentum
