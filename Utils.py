# Generic utility functions and constants

import datetime as dtm
from math import pi, sqrt, sin, cos, acos, atan2, degrees
import sys

import numpy as np

# ==============================================================================
# CONSTANTS

# Standard gravitational acceleration (m/s^2)
g0 = 9.80665

# Newton's universal gravitation constant
# 2014 CODATA-recommended value (m^3 kg^-1 s^-2)
GRAV = 6.67408e-11

# Astronomical unit (meters)
AU = 1.495978707e11

# J2000 epoch
J2000 = dtm.datetime(2000,1,1,11,58,56)

# Sidereal year in seconds
YR = 3600*24*365.242190402

# Other constants (in MKS)
KM = 1000

# ==============================================================================

# Parses a string containing a date & time (or just a date if date_only is True)
def parse_date(date_str, date_only=False):

  if date_str.startswith("A.D."):
    date_str = date_str.replace("A.D.", "")

  date_formats = [r"%Y-%m-%d", r"%Y-%b-%d", r"%Y-%B-%d", r"%b %d, %Y", r"%B %d, %Y"]
  time_formats = [r"%H:%M:%S.%f", r"%H:%M:%S", r"%H:%M", r"%H"]

  date_str = date_str.strip().lower()

  if date_only:
    fmts = date_formats
  else:
    fmts = [a + " " + b for a in date_formats for b in time_formats]

  dt = None
  for fmt in fmts:
    try:
      dt = dtm.datetime.strptime(date_str, fmt)
      break
    except:
      pass

  if dt is not None:
    if date_only:
      return dt.date()
    else:
      return dt
  else:
    print("Error: could not parse date/time '{}'".format(date_str))
    sys.exit()

# ==============================================================================
# Date ranges

def datetimes_range(start_date, end_date, N=200):
  td = dtm.timedelta(seconds=((end_date-start_date).total_seconds()/(N-1)))
  dts = []
  for i in range(0,N):
    dts.append(start_date + i*td)
  return dts

# ==============================================================================
# Date conversions

_MJD_epoch = dtm.datetime(1858, 11, 17, 0, 0, 0)
_MJD2000_epoch = dtm.datetime(2000, 1, 1, 12, 0, 0)

# Converts a Modified Julian Date to a Python datetime object
def mjd2dt(mjd):
  return _MJD_epoch + dtm.timedelta(days=float(mjd))

# Converts a Python datetime object to a Modified Julian Date
def dt2mjd(dt):
  return (dt - _MJD_epoch)/dtm.timedelta(days=1)
def dt2mjd2000(dt):
  return (dt - _MJD2000_epoch)/dtm.timedelta(days=1)

# ==============================================================================

# Returns the circular orbit velocity at given altitude (in km) around body
def circ_speed(body, altitude):
  return sqrt(body.mu_self/(body.radius+altitude))

# ==============================================================================
# Unit conversions

def fps_to_mps(fps):
  return fps * 0.3045

def mi_to_km(mi):
  return mi * 1.609

# ==============================================================================
# "Nice" formatters

# Input assumed in meters
def nice_distance(dist):
  if dist < 500e6:
    return "{:,.1f} km".format(dist/1e3)
  else:
    return "{:.3f} AU".format(dist/AU)

# Input assumed in seconds
def nice_time(_time):
  if _time < 2*3600:
    return "{:.1f} min".format(_time/60)
  elif _time < 24*3600:
    return "{:.1f} hours".format(_time/3600)
  elif _time < 365*24*3600:
    return "{:.1f} days".format(_time/(24*3600))
  else:
    return "{:.2f} yr".format(_time/(365.24*24*3600))

# ==============================================================================

# Calculates a unit vector normal to the 2D curve at the given index
# Will try to compute a centered difference. Requires at least two points.
def curve_normal(xs, ys, i0):
  assert len(xs) >= 2
  assert len(ys) >= 2
  x0 = xs[i0]; y0 = ys[i0]
  i1 = max(0, i0-1)
  x1 = xs[i1]; y1 = ys[i1]
  i2 = min(len(xs)-1, i0+1)
  x2 = xs[i2]; y2 = ys[i2]
  if x1 == x2:
    if y1 == y2:
      print("Error: curve normal is undefined")
      return None
    if x1 >= 0:
      nx = (1, 0)
    else:
      ny = (-1, 0)
    return np.array([nx, ny])
  else:
    m = (y2-y1)/(x2-x1)
    s = sqrt(1+m**2)
    nx = m/s
    ny = -1/s
    if nx*x0 < 0:
      nx *= -1
      ny *= -1
    return np.array([nx, ny])

# ==============================================================================

# Does a bisection search on a *sorted* array, returning the element and
# position index if the item is found, or the closest element (and index)
# if it is not.
# If passed, key must be a one-argument function that takes an item from
# the list and returns the key to be compared against the searched item.
def binary_search(items, searched, key=None):
  if key is None:
    key = lambda x: x
  N = len(items)
  imax = N - 1
  if searched <= key(items[0]):
    return (items[0], 0)
  elif searched >= key(items[imax]):
    return (items[imax], imax)
  i1 = 0
  i2 = imax
  while i2 - i1 > 1:
    if key(items[i1]) == searched:
      return (items[i1], i1)
    elif key(items[i2]) == searched:
      return (items[i2], i2)
    else:
      im = (i1+i2)//2
      if key(items[im]) == searched:
        return (items[im], im)
      elif key(items[im]) > searched:
        i2 = im
      else:
        i1 = im
  diff1 = abs(key(items[i1]) - searched)
  diff2 = abs(key(items[i2]) - searched)
  if diff1 <= diff2:
    return (items[i1], i1)
  else:
    return (items[i2], i2)

# ==============================================================================
# Angles

# Clamps an angle in degrees to the range [0, 360)
def clamp_degs(degrees, min_degs=0, max_dex=360):
  while degrees >= max_dex:
    degrees -= 360
  while degrees < min_degs:
    degrees += 360
  return degrees

# Clamps an angle in radians to the range [0, 2*pi)
def clamp_rads(radians):
  while radians >= 2*pi:
    radians -= 2*pi
  while radians < 0:
    radians += 2*pi
  return radians

# ==============================================================================

# Generates points along an ellipse in the xy plane with the given semi-major
# and semi-minor axes, which coincide with the x and y axes, respectively, and
# having its center at the origin.
# Step is in radians
def gen_ellipse(a, b, step=0.01):
  xs = []; ys = []
  theta = 0
  while theta < 2*pi:
    xs.append(a*cos(theta))
    ys.append(b*sin(theta))
    theta += step
  return xs, ys

# Generates points along a hyperbola in the xy plane with the given semi-major
# axis and eccentricity, with its focus at the origin, the axis of symmetry
# along the x-axis and opening to the left, with theta_max the maximum polar
# angle (in radians) measured at the focus from the vertex to the end points.
# Step is in radians
def gen_hyperbola(a, e, theta_max, step=0.01):
  xs = []; ys = []
  theta = -theta_max
  while theta <= theta_max:
    r = a*(e**2-1)/(1+e*cos(theta))
    xs.append(r*cos(theta))
    ys.append(r*sin(theta))
    theta += step
  return np.array(xs), np.array(ys)

# ==============================================================================
# 2D and 3D rotations

# Planar rotation of the given point through angle theta (degrees)
def rotate2D(x, y, theta):
  theta_rad = pi * theta / 180
  xp = x*cos(theta_rad) - y*sin(theta_rad)
  yp = x*sin(theta_rad) + y*cos(theta_rad)
  return xp, yp

# 3D rotation of the given point through angle theta (degrees) around
# axis specified by ux, uy, uz
def rotate3D(x, y, z, theta, ux, uy, uz):
  if (ux**2 + uy**2 + uz**2 - 1) > 1e-3:
    norm = sqrt(ux**2 + uy**2 + uz**2)
    ux /= norm; uy /= norm; uz /= norm
  theta_rad = pi * theta / 180
  C = cos(theta_rad)
  S = sin(theta_rad)
  xp = x*(C+ux**2*(1-C)) + y*(ux*uy*(1-C)-uz*S) + z*(ux*uz*(1-C)+uy*S)
  yp = x*(uy*ux*(1-C)+uz*S) + y*(C+uy**2*(1-C)) + z*(uy*uz*(1-C)-ux*S)
  zp = x*(uz*ux*(1-C)-uy*S) + y*(uz*uy*(1-C)+ux*S) + z*(C+uz**2*(1-C))
  return xp, yp, zp

# 3D elementary rotation around the x-axis
def Rx(x, y, z, theta):
  theta_rad = pi * theta / 180
  C = cos(theta_rad)
  S = sin(theta_rad)
  xp = x
  yp = C*y - S*z
  zp = S*y + C*z
  return xp, yp, zp

# 3D elementary rotation around the y-axis
def Ry(x, y, z, theta):
  theta_rad = pi * theta / 180
  C = cos(theta_rad)
  S = sin(theta_rad)
  xp = C*x + S*z
  yp = y
  zp = -S*x + C*z
  return xp, yp, zp

# 3D elementary rotation around the z-axis
def Rz(x, y, z, theta):
  theta_rad = pi * theta / 180
  C = cos(theta_rad)
  S = sin(theta_rad)
  xp = C*x - S*y
  yp = S*x + C*y
  zp = z
  return xp, yp, zp

# ==============================================================================
# VECTOR AND LINEAR ALGEBRA UTILS

# Projects coordinates (x,y,z) into plane defined by normal vector n
# Returns (x,y) coordinates on the plane
def project_to_plane(x, y, z, n):
  nx = n[0]; ny = n[1]; nz = n[2]
  rn = sqrt(nx**2 + ny**2 + nz**2)
  theta = acos(nz/rn) * 180/pi
  phi = atan2(ny, nx) * 180/pi
  x1, y1, z1 = Rz(x, y, z, -phi)
  x2, y2, z2 = Ry(x1, y1, z1, -theta)
  return x2, y2

# Returns difference between two vectors
def diff(v1, v2):
  return vec_diff(v1, v1)
def vec_diff(vec1, vec2):
  return np.array([(vec1[i]-vec2[i]) for i in range(len(vec1))])

# Returns distance between two vectors, i.e. magn(vec1-vec2)
def distance(vec1, vec2):
  return vec_dist(vec1, vec2)
def vec_dist(vec1, vec2):
  return magnitude(vec_diff(vec1, vec2))
  # return sqrt(sum([(vec1[i]-vec2[i])**2 for i in range(len(vec1))]))

# Returns the squared distance between two vectors, i.e. magn(vec1-vec2)**2
def vec_dist_squared(vec1, vec2):
  return sum([(vec1[i]-vec2[i])**2 for i in range(len(vec1))])

# Squared vector magnitude
def magnsq(vec):
  return sum([x**2 for x in vec])

# Vector magnitude / norm
def magn(vec):
  return magnitude(vec)
def norm(vec):
  return magnitude(vec)
def magnitude(vec):
  return sqrt(magnsq(vec))

# Returns the normal vector in the direction of vec
def normalize(vec):
  m = magnitude(vec)
  return np.array([x/m for x in vec])

# Returns the angle, in degrees, between two vectors, using their dot product
def angle_between(vec1, vec2):
  dotp = sum(vec1[i]*vec2[i] for i in range(len(vec1)))
  return degrees(acos(dotp/(magn(vec1)*magn(vec2))))


# ==============================================================================
