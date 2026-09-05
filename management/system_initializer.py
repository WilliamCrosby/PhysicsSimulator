from orbitalsim.celestial_body import CelestialBody

class SystemInitializer:

    @staticmethod
    def adjust_to_barycenter(bodies: list[CelestialBody]) -> None:
        total_mass = sum(b.mass for b in bodies)
        center_of_mass = sum(b.mass * b.position for b in bodies) / total_mass
        center_of_mass_velocity = sum(b.mass * b.velocity for b in bodies) / total_mass

        for b in bodies:
            b.position -= center_of_mass
            b.velocity -= center_of_mass_velocity
