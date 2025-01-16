# Plots the orbit given a TLE (3LE)

import matplotlib.pyplot as plt

from Orbit import Orbit
import SolarSystem
from SolarSystem import Bodies

# Some examples (with orbits updated March 2023):

## The International Space Station
TLE = """ISS (ZARYA)             
1 25544U 98067A   23066.65796090  .00063326  00000+0  11315-2 0  9991
2 25544  51.6402 111.4629 0006757  46.8215  99.0592 15.49401145386044"""

## A GPS satellite
# TLE = """NAVSTAR 73 (USA 260)    
# 1 40534U 15013A   23065.97345392  .00000030  00000+0  00000+0 0  9995
# 2 40534  53.5531 253.8092 0077821  24.6199 335.7514  2.00565805 58210"""

## GOES-16, a meteorological satellite in geostationary orbit
# TLE = """GOES 16                 
# 1 41866U 16071A   23066.40259302 -.00000249  00000+0  00000+0 0  9996
# 2 41866   0.0911 285.9379 0000672  69.8671 238.7723  1.00272472 23094"""

### Kosmos 2546, in a Molniya orbit
# TLE = """COSMOS 2546             
# 1 45608U 20031A   23063.38511294 -.00000416  00000+0  00000+0 0  9996
# 2 45608  63.1911 195.9777 6604588 266.6958 149.2155  2.00571036 20373"""

# ===========================================

orbit = Orbit(primary=Bodies["Earth"])
orbit.from_TLE(TLE)
plt.figure(figsize=(8,8))
orbit.plot(show_primary=True, projection="3D", units=1e6)
plt.tight_layout()
plt.show()
