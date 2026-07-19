class CelestialBody:

    def __init__(self, mass, position, velocity, name, horizons_id=-1):
        self.mass = mass
        self.position = position
        self.velocity = velocity
        self.name = name
        self.horizons_id = horizons_id
