class CelestialBody:

    def __init__(self, mass, position, velocity, name, horizons_id=None, barycenter_id=None):
        self.mass = mass
        self.position = position
        self.velocity = velocity
        self.name = name
        self.horizons_id = horizons_id
        self.barycenter_id = barycenter_id
