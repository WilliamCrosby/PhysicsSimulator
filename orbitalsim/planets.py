import numpy as np
from orbitalsim.celestial_body import CelestialBody
from orbitalsim.universal_constants import UniversalConstants

class Planets:

    @staticmethod
    def CircularVelocities(r: float):
        return np.sqrt(UniversalConstants.G * 1.989e30 / r)

    mercury = CelestialBody(
        3.301e23,
        np.array([0.39 * UniversalConstants.AU, 0.0, 0.0]),
        np.array([0.0, CircularVelocities(0.39 * UniversalConstants.AU), 0.0]),
        name="Mercury"
    )

    venus = CelestialBody(
        4.867e24,
        np.array([0.72 * UniversalConstants.AU, 0.0, 0.0]),
        np.array([0.0, CircularVelocities(0.72 * UniversalConstants.AU), 0.0]),
        name="Venus"
    )

    earth = CelestialBody(
        5.972e24,
        np.array([UniversalConstants.AU, 0.0, 0.0]),
        np.array([0.0, CircularVelocities(UniversalConstants.AU), 0.0]),
        name="Earth"
    )

    mars = CelestialBody(
        6.417e23,
        np.array([1.52 * UniversalConstants.AU, 0.0, 0.0]),
        np.array([0.0, CircularVelocities(1.52 * UniversalConstants.AU), 0.0]),
        name="Mars"
    )

    jupiter = CelestialBody(
        1.898e27,
        np.array([5.20 * UniversalConstants.AU, 0.0, 0.0]),
        np.array([0.0, CircularVelocities(5.20 * UniversalConstants.AU), 0.0]),
        name="Jupiter"
    )

    saturn = CelestialBody(
        5.683e26,
        np.array([9.58 * UniversalConstants.AU, 0.0, 0.0]),
        np.array([0.0, CircularVelocities(9.58 * UniversalConstants.AU), 0.0]),
        name="Saturn"
    )

    uranus = CelestialBody(
        8.681e25,
        np.array([19.2 * UniversalConstants.AU, 0.0, 0.0]),
        np.array([0.0, CircularVelocities(19.2 * UniversalConstants.AU), 0.0]),
        name="Uranus"
    )

    neptune = CelestialBody(
        1.024e26,
        np.array([30.05 * UniversalConstants.AU, 0.0, 0.0]),
        np.array([0.0, CircularVelocities(30.05 * UniversalConstants.AU), 0.0]),
        name="Neptune"
    )

    allPlanets = [mercury, venus, earth, mars, jupiter, saturn, uranus, neptune]
