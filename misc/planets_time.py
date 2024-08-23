#!/usr/bin/env python3
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import solar_system_ephemeris, EarthLocation, Angle, SkyCoord
from astropy.coordinates import get_body_barycentric, get_body, get_moon
import sys

if len(sys.argv) != 3:
    sys.exit("Usage %s YYYY-MM-DD hh:mm:ss" %(sys.argv[0]))

# Time is in format 2021-05-26 06:06:00
# Where are we observing from?
latitude = Angle("-26:41:46.0", unit=u.deg)
longitude = Angle("116:38:13.0", unit=u.deg)
observing_location = EarthLocation(lat=latitude, lon=longitude)

# When are we observing?
t = Time("%s %s" %(sys.argv[1], sys.argv[2]), format='iso', scale='utc')

print("Solar System objects as of %s" %(t))
with solar_system_ephemeris.set('builtin'):
    for planet in ["sun", "moon", "mercury", "venus", "mars", "jupiter", "saturn", "uranus", "neptune"]:
        planet_sc = get_body(planet, t, observing_location)
        dir_str = planet_sc.to_string(style='hmsdms')
        print("%7s : %s" %(planet, dir_str))
