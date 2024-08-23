#!/usr/bin/env python3
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import solar_system_ephemeris, EarthLocation, Angle, SkyCoord
from astropy.coordinates import get_body_barycentric, get_body, get_moon
import datetime
# Where are we observing from?
latitude = Angle("-26:41:46.0", unit=u.deg)
longitude = Angle("116:38:13.0", unit=u.deg)
observing_location = EarthLocation(lat=latitude, lon=longitude)

# Get the current local time
t = Time(datetime.datetime.now(datetime.UTC), scale='utc', location=observing_location)
t.format = 'iso'
print("Solar System objects as of %s" %(t))
with solar_system_ephemeris.set('builtin'):
    for planet in ["sun", "moon", "mercury", "venus", "mars", "jupiter", "saturn", "uranus", "neptune"]:
        planet_sc = get_body(planet, t, observing_location)
        dir_str = planet_sc.to_string(style='hmsdms')
        print("%8s : %s" %(planet, dir_str))
