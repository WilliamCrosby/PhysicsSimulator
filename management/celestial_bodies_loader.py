import json
import numpy as np
from orbitalsim.celestial_body import CelestialBody
from orbitalsim.universal_constants import UniversalConstants


class CelestialBodyLoader:

    GM_UNITS_TO_SI = {
        "m3/s2": 1.0,
        "km3/s2": 1e9,
    }

    @classmethod
    def _gm_to_si(cls, gm: float, unit: str) -> float:
        try:
            return gm * cls.GM_UNITS_TO_SI[unit]
        except KeyError as exc:
            raise ValueError(f"Unsupported GM unit: {unit}") from exc

    @staticmethod
    def _compute_circular_velocity(r: float, central_gm: float) -> float:
        return np.sqrt(central_gm / r)

    @classmethod
    def from_dict(cls, data: dict, central_gm: float = 1.32712440018e20) -> CelestialBody:
        name = data["name"]
        gm = cls._gm_to_si(data["gm"], data.get("gm_unit", "m3/s2"))
        mass = gm / UniversalConstants.G
        horizons_id = data.get("horizons_id")
        barycenter_id = data.get("barycenter_id")

        if "position_au" in data:
            pos = np.array(data["position_au"], dtype=float) * UniversalConstants.AU_TO_M
        else:
            pos = np.array(data["position"], dtype=float)

        if data.get("circular_orbit", False):
            r = np.linalg.norm(pos)
            v_mag = cls._compute_circular_velocity(r, central_gm)
            vel = np.array([0.0, v_mag, 0.0], dtype=float)
        else:
            vel = np.array(data["velocity"], dtype=float)

        return CelestialBody(
            GM=gm,
            mass=mass,
            position=pos,
            velocity=vel,
            name=name,
            horizons_id=horizons_id,
            barycenter_id=barycenter_id
        )

    @classmethod
    def load_from_json(cls, json_path: str) -> list[CelestialBody]:
        with open(json_path, 'r', encoding='utf-8') as f:
            raw_bodies = json.load(f)

        sun_data = next((b for b in raw_bodies if b["name"] == "Sun"), None)
        central_gm = (
            cls._gm_to_si(sun_data["gm"], sun_data.get("gm_unit", "m3/s2"))
            if sun_data
            else 1.32712440018e20
        )

        return [cls.from_dict(b_data, central_gm=central_gm) for b_data in raw_bodies]
