import numpy as np
from orbitalsim.celestial_body import CelestialBody
from orbitalsim.universal_constants import UniversalConstants


def MechanicalEnergyCalculator(bodies:list[CelestialBody], positions, velocities):

    gravitationalPotentialEnergies = []
    kineticEnergies = []

    totalKineticEnergy = 0.0
    totalGravitationalPotentialEnergy = 0.0

    for i in range(len(bodies)):

        kineticEnergies.append(0.5 * bodies[i].mass * np.linalg.norm(velocities[i])**2)

        GPE = 0.0

        for j in range(i + 1, len(bodies)):

            if j == i:
                continue

            GPE += (-1 * UniversalConstants.G * bodies[j].mass * bodies[i].mass / np.linalg.norm(positions[i] - positions[j]))

        gravitationalPotentialEnergies.append(GPE)

    for i in range(len(gravitationalPotentialEnergies)):
        totalGravitationalPotentialEnergy += gravitationalPotentialEnergies[i]
    for j in range (len(kineticEnergies)):
        totalKineticEnergy += kineticEnergies[j]

    return totalKineticEnergy + totalGravitationalPotentialEnergy


def MomentumCalculator(bodies:list[CelestialBody], positions, velocities):

    totalLinearMomentum = np.zeros(3)
    totalAngularMomentum = np.zeros(3)

    for i in range(len(bodies)):
        p = bodies[i].mass * velocities[i]

        totalLinearMomentum += p
        totalAngularMomentum += np.cross(positions[i], p)

    return totalLinearMomentum, totalAngularMomentum


def UpdateExtrema(bodies, positions, velocities, maximumEnergy, minimumEnergy, maximumLinearMomentum, minimumLinearMomentum, maximumAngularMomentum, minimumAngularMomentum):
    if MechanicalEnergyCalculator(bodies, positions, velocities) > maximumEnergy:
        maximumEnergy = MechanicalEnergyCalculator(bodies, positions, velocities)
    if MechanicalEnergyCalculator(bodies, positions, velocities) < minimumEnergy:
        minimumEnergy = MechanicalEnergyCalculator(bodies, positions, velocities)
    if np.linalg.norm(MomentumCalculator(bodies, positions, velocities)[0]) > maximumLinearMomentum:
        maximumLinearMomentum = np.linalg.norm(MomentumCalculator(bodies, positions, velocities)[0])
    if np.linalg.norm(MomentumCalculator(bodies, positions, velocities)[0]) < minimumLinearMomentum:
        minimumLinearMomentum = np.linalg.norm(MomentumCalculator(bodies, positions, velocities)[0])
    if np.linalg.norm(MomentumCalculator(bodies, positions, velocities)[1]) > maximumAngularMomentum:
        maximumAngularMomentum = np.linalg.norm(MomentumCalculator(bodies, positions, velocities)[1])
    if np.linalg.norm(MomentumCalculator(bodies, positions, velocities)[1]) < minimumAngularMomentum:
        minimumAngularMomentum = np.linalg.norm(MomentumCalculator(bodies, positions, velocities)[1])

    return maximumEnergy, minimumEnergy, maximumLinearMomentum, minimumLinearMomentum, maximumAngularMomentum, minimumAngularMomentum
