# A class to represent general 3D orbits, with the ability to plot them
# (using matplotlib) or compute position and/or velocity at arbitrary
# times

from datetime import datetime, timedelta
from math import sqrt, pi, sin, cos, tan, atan, atan2, acos
import sys

import numpy as np
from numpy.linalg import norm

import Utils

# ==============================================================================

class Orbit:

  # The only required argument is the primary, which must be a CelestialBody
  # instance.
  # The optional arguments elements or state_vectors specify the orbit:
  #   elements must be a dict containining the elements in the form specified
  #     by from_elements()
  #   state_vectors must be a list/tuple with two lists/tuples containing the
  #     3D cartesian position and velocity
  # If both are given, only elements is used.
  def __init__(self, primary, elements=None, state_vectors=None):

    self.a = None
    self.e = None
    self.i = None
    self.w = None
    self.LAN = None
    self.M0 = None
    self.epoch = None
    self.long_peri = None
    self.L0 = None
    self.b = None
    self.c = None
    self.P = None
    self.n = None
    self.rpe = None
    self.rap = None
    self.vpe = None
    self.vap = None
    self.hap = None
    self.hpe = None
    self.vinf = None
    self.C3 = None
    self.th_inf = None
    self.phi = None
    self.equatorial = None
    self.retrograde = None

    self.fig = None
    self.ax = None
    self.color = None
    self.projection = "2D"
    self.units = 1.0
    self.l_scl = 1.0

    self.primary = primary

    if elements is not None:
      self.from_elements(elements)
    elif state_vectors is not None:
      self.from_state_vectors(*state_vectors)

  # ----------------------------------------------------------------------------

  def __str__(self):
    if self.hpe is None or self.hap is None or self.e is None:
      return "<Orbit around {} (uninitialized)>".format(self.primary.name)
    else:
      if self.e < 1:
        if self.e == 0: kind = "circular"
        elif self.e < 0.08: kind = "near-circular"
        else: kind = "elliptical"
        return "<Orbit around {}, {}, {} x {} alt, {:.1f}° inc, P={}>".format(self.primary.name, kind, Utils.nice_distance(self.hpe), Utils.nice_distance(self.hap), self.i, Utils.nice_time(self.P))
      else:
        return "<Orbit around {}, hyperbolic, hpe={}, vinf={:.1f} km/s, C3={:.1f} km^2/s^2>".format(self.primary.name, Utils.nice_distance(self.hpe), self.vinf/1e3, self.C3/1e6)

  def __repr__(self):
    if self.hpe is None or self.hap is None or self.e is None:
      return "<Orbit, primary={} (uninitialized)>".format(self.primary.name)
    else:
      if self.e < 1:
        if self.e == 0: kind = "circular"
        elif self.e < 0.08: kind = "near-circular"
        else: kind = "elliptical"
        return "<Orbit, primary={}, a={:.6e} m, e={:.6f}, i={:.1f}°, LAN={:.1f}°, w={:.1f}°, M0={:.1f}°, epoch={})>".format(self.primary.name, self.a, self.e, self.i, self.LAN, self.w, self.M0, self.epoch)
      else:
        return "<Orbit, primary={}, a={:.1f} m, e={:.6f}, i={:.1f}°, LAN={:.1f}°, w={:.1f}°, M0={:.1f}°, epoch={}>".format(self.primary.name, self.a, self.e, self.i, self.LAN, self.w, self.M0, self.epoch)

  # ----------------------------------------------------------------------------

  def from_elements(self, elements, to_epoch=None):

    """Defines the orbit from the "standard" set of orbital elements, which
    must be a dictionary of the form
    {'a': a, 'e': e, 'i': i, 'arg': arg, 'LAN': LAN, 'M0': M0, 'epoch': epoch}
    containing the following keys/values:
      'a': semi-major axis (meters)
      'e': eccentricity
      'i': inclination (degrees)
      'arg': argument of periapsis (degrees)
      'LAN': longitude of the ascending node (degrees)
      'M0': mean anomaly at given epoch (degrees)
      'epoch': UTC time at which the elements are given (tz-naive datetime)
    If the optional parameter to_epoch is given (a tz-naive datetime),
    the mean anomaly at epoch will be adjusted to that."""

    if to_epoch is not None:
      epoch = to_epoch
    else:
      epoch = elements['epoch']

    self.epoch = epoch
    self.a = elements['a']
    self.e = elements['e']
    self.i = elements['i']
    self.w = elements['arg']
    self.LAN = elements['LAN']
    self.M0 = elements['M0']

    if to_epoch is not None:
      if e > 1: n = sqrt(-self.primary.mu/self.a^3)
      else: n = sqrt(self.primary.mu/self.a^3)
      self.M0 -= n*self.secs_since_epoch(epoch)

    self.set_derived_params()

  # ----------------------------------------------------------------------------

  def from_state_vectors(self, pos, vel, date, to_epoch=None):

    """Defines the orbit from the 3D state vectors:
    pos: 3D cartesian position vector, in m
    vel: 3D cartesian velocity vector, in m/s
    date: UTC date as which the state is given (tz-naive datetime objects)
    to_epoch: (optional) UTC date to displace the saved epoch to (tz-naive
      datetime object)
    If to_epoch is given, the mean anomaly at epoch will be adjusted to that.
    If 2D vectors for pos and vel are given, a z coordinate of zero will be
    assumed."""

    if to_epoch is None:
      self.epoch = date
    else:
      self.epoch = to_epoch

    pos = np.array(pos)
    vel = np.array(vel)
    if len(pos) == 2: pos = np.append(pos, 0)
    if len(vel) == 2: vel = vel.append(pos, 0)

    r = norm(pos)
    v = norm(vel)

    # Semi-major axis (negative for hyperbolic orbits)
    a = 1/(2/r - v**2/self.primary.mu)

    # Angular momentum vector (specific)
    hvec = np.cross(pos, vel)
    h = norm(hvec)
    if h == 0:
      raise(Exception("Radial orbits (h=0) not implemented!"))

    # Eccentricity vector (points to periapsis)
    evec = np.cross(vel, hvec)/self.primary.mu - pos/r
    e = norm(evec)

    # Inclination
    i = acos(hvec[2]/h)
    if i == -pi: i = pi
    eps = 1e-10
    self.retrograde = (abs(i-pi) < eps)

    # Line of nodes vector (points to ascending node)
    self.equatorial = (abs(i) < eps or abs(i-pi) < eps)
    if self.equatorial:
      nvec = np.array([1,0,0])
    else:
      nvec = np.cross(np.array((0,0,1)), hvec)
    n = norm(nvec)

    # Longitude of the ascending node
    if self.equatorial:
      LAN = 0
    else:
      LAN = acos(nvec[0]/n)
      if nvec[1] < 0:
        LAN = 2*pi - LAN

    # Argument of periapsis
    w = acos(np.dot(nvec, evec)/(n*e))
    if self.equatorial:
      if evec[1] < 0:
        w = 2*pi - w
    elif evec[2] < 0:
      w = 2*pi - w

    # True anomaly
    foo = np.dot(evec, pos)/(e*r)
    if foo > 1 and (1-foo) < 1e-5:
      foo = 1
    elif foo < -1 and (foo+1) < 1e5:
      foo = -1
    nu = acos(foo)
    if np.dot(pos,vel) < 0:
      nu = 2*pi - nu

    # Mean anomaly
    if e < 1:
      E = 2*atan(sqrt((1-e)/(1+e))*tan(nu/2))
      M = E - e*sin(E)
    elif e > 1:
      F = 2*np.arctanh(sqrt((e-1)/(e+1))*tan(nu/2))
      M = e*np.sinh(F) - F
    else:
      raise(Exception("Parabolic orbits not implemented!"))

    # Set mean anomaly at epoch (epoch may be different from date)
    if to_epoch is None:
      M0 = M
    else:
      if e > 1:
        motion = sqrt(-self.primary.mu/self.a^3)
      else:
        motion = sqrt(self.primary.mu/self.a^3)
      M0 = M - motion*self.secs_since_epoch(date)

    # Clamp M0 to [0, 360) if elliptic; leave sign if hyperbolic
    if e < 1:
      M0 = Utils.clamp_degs(np.degrees(M0))
    else:
      M0 = np.degrees(M0)

    self.a = a
    self.e = e
    self.i = np.degrees(i)
    self.w = np.degrees(w)
    self.LAN = np.degrees(LAN)
    self.M0 = M0

    self.set_derived_params()

  # ----------------------------------------------------------------------------

  def from_TLE(self, TLE):

    """Defines the orbit from a TLE (two-line element set) or 3LE

    The argument TLE can be either a string, with newlines characters separating
    the (2 or 3) lines, or a 2- or 3-element tuple with the lines."""

    if isinstance(TLE, str):
      TLE = TLE.split("\n")
    assert len(TLE) in [2,3]

    if len(TLE) == 3:
      line0 = TLE[0]
      line1 = TLE[1]
      line2 = TLE[2]
    else:
      line0 = None
      line1 = TLE[0]
      line2 = TLE[1]

    epoch_year = int(line1[18:20])
    if epoch_year < 57: epoch_year += 2000
    else: epoch_year + 1900
    epoch_day = float(line1[20:32])

    epoch = datetime(epoch_year, 1, 1) + timedelta(days=epoch_day-1)
    inc = float(line2[8:16])
    RAAN = float(line2[17:25])
    ecc = float("0."+line2[26:33])
    w = float(line2[34:42])
    M0 = float(line2[43:51])
    motion = float(line2[52:63])   # in revolutions per day

    P = 86400/motion
    a = (self.primary.mu*(P/(2*pi))**2)**(1.0/3)

    elements = {'a': a, 'e': ecc, 'i': inc, 'LAN': RAAN, 'arg': w, 'M0': M0, 'epoch': epoch}

    self.from_elements(elements)

  # ----------------------------------------------------------------------------

  def set_derived_params(self):
    """"Sets derived orbital parameters from main elements"""

    self.long_peri = Utils.clamp_degs(self.w + self.LAN)
    self.L0 = Utils.clamp_degs(self.M0 + self.long_peri)
    self.E = -self.primary.mu/(2*self.a)

    if self.e < 1:
      
      self.b = self.a*sqrt(1-self.e**2)
      self.c = self.e * self.a
      self.n = sqrt(self.primary.mu/self.a**3)
      self.P = 2*pi/self.n
      self.rpe = (1-self.e)*self.a
      self.rap = (1+self.e)*self.a
      self.hpe = self.rpe - self.primary.radius
      self.hap = self.rap - self.primary.radius
      self.vpe = self.speed_at_radius(self.rpe)
      self.vap = self.speed_at_radius(self.rap)

    elif self.e > 1:
      
      self.b = -self.a*sqrt(self.e**2-1)
      self.c = -self.a*self.e
      self.l = -self.a*(self.e**2-1)
      self.n = sqrt(-self.primary.mu/self.a**3)
      self.P = np.inf
      self.rpe = -self.a*(self.e-1)
      self.rap = np.inf
      self.hpe = self.rpe - self.primary.radius
      self.hap = np.nan
      self.vpe = self.speed_at_radius(self.rpe)
      self.vap = np.nan
      self.vinf = sqrt(-self.primary.mu/self.a)
      self.C3 = self.vinf**2
      self.th_inf = np.degrees(acos(-1/self.e))
      self.phi = 180 - self.th_inf

    else:
      raise(Exception("Parabolic orbits not implemented"))

  def get_elements(self):
    return {"a": self.a, "e": self.e, "i": self.i, "arg": self.w, "LAN": self.LAN, "M0": self.M0, "long_peri": self.long_peri, "L0": self.L0, "epoch": self.epoch}

  # ----------------------------------------------------------------------------

  def set_axes(self, ax):
    self.ax = ax

  def set_projection(self, projection):
    self.projection = projection

  # ----------------------------------------------------------------------------
  # Orbital mechanics utility functions

  # Returns the (fractional) number of seconds since the epoch
  # Negative it before epoch
  def secs_since_epoch(self, date):
    return (date - self.epoch).total_seconds()

  # Speed at given radius (i.e. vis-viva equation)
  def speed_at_radius(self, radius):
    assert self.rpe <= radius <= self.rap
    return sqrt(self.primary.mu*(2/radius - 1/self.a))

  # Radius at given true anomaly
  def radius_at_nu(self, nu):
    return self.a*(1-self.e**2)/(1+self.e*cos(np.radians(nu)))

  # Mean anomaly M at provided date, in degrees
  # date must be a UTC tz-naive datetime object
  def M_at_date(self, date):
    return np.degrees(self.n*self.secs_since_epoch(date)) + self.M0

  # Returns a 2-tuple with the (x,y,z) positions of the periapsis and apoapsis
  # Will return None for the apoapsis if orbit if hyperbolic
  def get_apsides(self):
    return np.array((self.get_periapsis(), self.get_apoapsis()))

  # Returns the (x,y,z) position of the periapsis
  def get_periapsis(self, vel=False):
    xpe, ype, zpe = self.perifocal_to_inertial(self.rpe, 0, 0)
    xpe /= self.units; ype /= self.units; zpe /= self.units
    return xpe, ype, zpe

  # Returns the (x,y,z) position of the apoapsis
  def get_apoapsis(self, vel=False):
    xap, yap, zap = self.perifocal_to_inertial(-self.rap, 0, 0)
    xap /= self.units; yap /= self.units; zap /= self.units

    if vel:
      return xap, yap, zap, vxap, vyap, vzap
    else:
      return xap, yap, zap

  # Change periapsis
  def dv_peri_change(self, new_hpe):
    new_orb = Orbit(self.primary)
    return new_orb.vap - self.vap

  # ----------------------------------------------------------------------------

  # Returns the position and velocity vectors (in the inertial frame) at
  # the arbitrary date given, which must be a tz-naive datetime object
  # specified in the same time system as the epoch
  # Can also return perifocal position and velocity if perifocal=True
  # This is a wrapper
  def posvel_at_date(self, date, perifocal=False, dist_units=None, vel_units=None):

    nu = self.nu_at_date(date)
    pos, vel = self.posvel_at_nu(nu)
    
    if dist_units is not None:
      pos /= dist_units
    if vel_units is not None:
      vel /= vel_units

    return (pos, vel)
  
  # ----------------------------------------------------------------------------

  # Solves the kepler equation for the true anomaly at the given date, which
  # must be a tz-naive datetime object specified in the same time system as
  # the epoch.
  def nu_at_date(self, date):

    # Compute mean anomaly at specified date (in radians)
    M = np.radians(self.M_at_date(date))

    # Solve the Kepler equation and obtain the true anomaly
    if self.e < 1:
      E = self.solveKepler(M, self.e)
      # nu = 2*atan(sqrt((1+self.e)/(1-self.e))*tan(E/2))
      nu = 2*atan2(sqrt(1+self.e)*sin(E/2), sqrt(1-self.e)*cos(E/2))
    else:
      F = self.solveKepler(M, self.e)
      # nu = acos((self.e - np.cosh(F))/(self.e*np.cosh(F)-1))
      # nu = 2*atan(sqrt((self.e+1)/(self.e-1))*np.tanh(F/2))
      nu = 2*atan2(sqrt(self.e+1)*np.sinh(F/2), sqrt(self.e-1)*np.cosh(F/2))
      # if F < 0:
      #   nu = -np.abs(nu)

    return nu

  # ----------------------------------------------------------------------------

  # Returns the position and velocity vectors (in the inertial frame) at
  # the given true anomaly nu, in radians
  # Can also return perifocal position and velocity if perifocal=True
  def posvel_at_nu(self, nu, perifocal=False, dist_units=None, vel_units=None):

    # Perifocal position
    p = self.a*(1-self.e**2)
    r = p/(1 + self.e*cos(nu))
    x = r*cos(nu)
    y = r*sin(nu)
    z = 0

    # Perifocal velocity
    vz = 0
    if self.e < 1:
      # vx = sqrt(self.primary.mu*self.a)/r*(-sin(E))
      # vy = sqrt(self.primary.mu*self.a)/r*(sqrt(1-self.e**2)*cos(E))
      vx = -sqrt(self.primary.mu/p)*sin(nu)
      vy = sqrt(self.primary.mu/p)*(self.e + cos(nu))
    else:
      # vx = sqrt(-self.primary.mu*self.a)/r*(-np.sinh(F))
      # vy = sqrt(-self.primary.mu*self.a)/r*(sqrt(self.e**2-1)*np.cosh(F))
      vx = -sqrt(self.primary.mu/p)*sin(nu)
      vy = sqrt(self.primary.mu/p)*(self.e + cos(nu))

    if perifocal:
      return np.array((x, y, z)), np.array((vx, vy, vz))

    # Transform to inertial reference frame
    x, y, z = self.perifocal_to_inertial(x, y, z)
    vx, vy, vz = self.perifocal_to_inertial(vx, vy, vz)

    # Scale distance and/or velocity if units are provided
    if dist_units is not None:
      x, y, z = x/dist_units, y/dist_units, z/dist_units
    if vel_units is not None:
      vx, vy, vz = vx/vel_units, vy/vel_units, vz/vel_units

    return np.array((x, y, z)), np.array((vx, vy, vz))

  # ----------------------------------------------------------------------------

  # Applies the coordinate transform from perifocal coordinates to
  # the inertial coordinates in which the elements are defined (e.g.
  # equatorial coordinates)
  def perifocal_to_inertial(self, x, y, z):
    if self.equatorial and self.retrograde:
      x, y, z = Utils.Rx(x, y, z, 180)
    x, y, z = Utils.Rz(x, y, z, self.w)
    if not self.equatorial:
      x, y, z = Utils.Rx(x, y, z, self.i)
      x, y, z = Utils.Rz(x, y, z, self.LAN)
    return x, y, z

  # ----------------------------------------------------------------------------

  # Solves Kepler's equation, M = E - e*sin(E), for the eccentric anomaly E,
  # given M and e; if e > 1, solves the hyperbolic Kepler equation instead,
  # M = e*sinh(F) - F, for the hyperbolic eccentric anomaly F (although it's
  # still called E in the code below).
  def solveKepler(self, M, e):

    # Tolerance (absolute) for Newton's method
    tol = 1e-7
    max_iter = 1000

    if (M == 0):
      return 0

    # Starting values from Vallado
    if e < 1:
      if -pi < M < 0 or M > pi:
        E = M - e
      else:
        E = M + e
    elif e > 1:
      if e < 1.6:
        if -pi < M < 0 or M > pi:
          E = M - e
        else:
          E = M + e
      else:
        if e < 3.6 and abs(M) > pi:
          E = M - e*(1 if M >= 0 else -1)
        else:
          E = M/(e-1)
    else:
      raise Exception("Kepler solver for parabolic orbits not implemented!")

    err = tol*2
    iter = 1
    while (err > tol):
      if (e < 1):
        Enew = E - (E-e*sin(E)-M)/(1-e*cos(E))
      else:
        Enew = E - (e*np.sinh(E)-E-M)/(e*np.cosh(E)-1)
      err = abs(Enew - E)
      iter += 1
      if iter > max_iter:
        print("\nKepler solver reached max iterations!")
        print("M=", M)
        print("e=", e)
        print("current E=", Enew)
        print("last E=", E)
        print("err=", err)
        raise Exception("solveKepler() did not converge!")
      E = Enew

    return E

  # ----------------------------------------------------------------------------

  # Imports matplotlib modules for potting
  def _matplotlib_imports(self):
    global plt, Axes3D
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

  # Sets the axes uses fot plotting (will create a new figure/axes if needed)
  def _init_axes(self, projection=None):
    if projection is None:
      projection = self.projection
    print(len(plt.get_fignums()))
    if len(plt.get_fignums()) == 0:
      self.fig = plt.figure()
    else:
      self.fig = plt.gcf()
    if len(self.fig.axes) == 0:
      if self.projection == "3D":
        self.ax = self.fig.add_subplot(111, projection='3d')
      elif self.projection == "2D":
        self.ax = self.fig.add_subplot(111)
    else:
      self.ax = plt.gca()

  # Plots the position of the body at the given date
  # date must be a UTC tz-naive datetime object
  def plot_at_date(self, date, projection=None, ax=None, units=None, opts=None, label=None, label_opts=None):

    self._matplotlib_imports()

    if ax is not None:
      self.ax = ax
    if projection is not None:
      self.projection = projection
    if units is not None:
      self.units = units
    if opts is None:
      opts = {}
    if label_opts is None:
      label_opts = {}

    if self.ax is None:
      self.init_axes()

    if "color" not in opts and self.color is not None:
      opts["color"] = self.color

    if "zorder" not in opts:
      opts["zorder"] = 10

    (x, y, z), (vx, vy, vz) = self.posvel_at_date(date)
    x /= self.units; y /= self.units; z /= self.units

    if self.projection == "3D":
      ln = self.ax.scatter([x],[y],[z], **opts)
    elif self.projection == "2D":
      if "s" not in opts:
        opts["s"] = 30
      ln = self.ax.scatter([x],[y], **opts)
    
    if label is not None:
      r = sqrt(x**2+y**2+z**2)
      nx, ny, nz = x/r, y/r, z/r
      if self.projection == "3D":  
        self.ax.text(x, y, z, label, color=opts["color"])
      elif self.projection == "2D":
        if "xytext" in label_opts:
          xytext = label_opts["xytext"]
          label_opts["textcoords"] = "offset pixels"
        else:
          d = 25
          label_opts["xytext"] = (d*nx, d*ny)
          label_opts["textcoords"] = "offset pixels"
        plt.annotate(label, xy=(x,y), ha="center", va="center", color=opts["color"], **label_opts)
        # bbox = dict(fc="black", alpha=0.7, pad=0, lw=0)
        # plt.annotate(label, xy=(x,y), xytext=(13, 0), textcoords="offset pixels", ha="left", va="center", color=opts["color"], bbox=bbox, **label_opts)

    return ln

  # ----------------------------------------------------------------------------

  # Generates the orbit data
  def gen_orbit(self, units=None, plot_step=0.01, rmax=None):

    if units is not None:
      self.units = units

    # Generate points along trajectory in perifocal coords
    xs = []; ys = []; zs = []
    if self.e < 1:
      xs1, ys1 = Utils.gen_ellipse(self.a, self.b, step=plot_step)
      xs1 = [x-self.c for x in xs1]
    elif self.e > 1:
      if rmax is None:
        rmax = abs(self.a)*20
      theta_max = np.arccos((self.l/rmax-1)/self.e)
      xs1, ys1 = Utils.gen_hyperbola(abs(self.a), self.e, theta_max, step=plot_step)
    
    # Transform to standard coords
    for i in range(len(xs1)):
      x, y, z = self.perifocal_to_inertial(xs1[i], ys1[i], 0)
      x /= self.units; y /= self.units; z /= self.units
      xs.append(x); ys.append(y); zs.append(z)

    return np.array(xs), np.array(ys), np.array(zs)

  # ----------------------------------------------------------------------------

  # Plots the orbit
  # Can optionally provide the axes
  def plot(self, ax=None, projection=None, normal=None, show_apsides=False, show_line_apsides=False, show_nodes=False, show_axes=False, show_primary=False, units=None, plot_step=0.01, rmax=None, dark=True, opts=None):

    self._matplotlib_imports()

    # if self.e > 1:
    #   print("Can't plot hyperbolic orbits yet!")
    #   return

    # For hyperbolic orbits, compute a decent maximum radius if not given
    if rmax is None:
      rmax = (-self.a)*20

    if opts is None:
      opts = {}

    if ax is not None:
      self.ax = ax
    if projection is not None:
      self.projection = projection.upper()
    if units is not None:
      self.units = units

    if self.projection == "2D":
      if normal is None:
        normal = np.array([0, 0, 1])

    if self.ax is None:
      self._init_axes()

    if "color" not in opts and self.color is not None:
      opts["color"] = self.color

    # Generate orbit data points
    if self.e < 1:
      xs, ys, zs = self.gen_orbit(units=units, plot_step=plot_step)
    elif self.e > 1:
      xs, ys, zs = self.gen_orbit(units=units, plot_step=plot_step, rmax=rmax)

    # Plot orbit
    if self.projection == "3D":
      ln, = self.ax.plot(xs, ys, zs, **opts)
    elif self.projection == "2D":
      ln, = self.ax.plot(xs, ys, **opts)

    if self.e < 1:
      scale = self.rap*1.2 / self.units
    elif self.e > 1:
      scale = rmax*1.2 / self.units
    self.ax.set_xlim(-scale, scale)
    self.ax.set_ylim(-scale, scale)
    if self.projection == "3D":
      self.ax.set_zlim(-scale, scale)

    # Set equal aspect ratio
    if self.projection == "3D":
      self.ax.set_box_aspect((1, 1, 1))
    elif self.projection == "2D":
      self.ax.set_aspect('equal')

    # Primary as a circle/sphere
    if show_primary:

      R = self.primary.radius/self.units

      if self.projection == "3D":
        
        u = np.linspace(0, 2 * np.pi, 40)
        v = np.linspace(0, np.pi, 30)
        x = R * np.outer(np.cos(u), np.sin(v))
        y = R * np.outer(np.sin(u), np.sin(v))
        z = R * np.outer(np.ones(np.size(u)), np.cos(v))
        self.ax.plot_surface(x, y, z, color="gray", alpha=0.7)

      elif self.projection == "2D":
        
        import matplotlib.patches as mpatches
        circle = mpatches.Circle((0,0), R, fill=True, color="gray", zorder=10)
        plt.gca().add_artist(circle)

    # else:
    #
    #   # Just at mark at the center of gravity
    #   if self.projection == "3D":
    #     self.ax.plot([0], [0], [0], "x", color="gray")
    #   elif self.projection == "2D":
    #     self.ax.plot([0], [0], "x", color="gray")

    # Show apsides as ticks and/or line of apsides
    if show_apsides or show_line_apsides:

      color = "gray" if dark else "k"

      # Eyeball tick_len by order of magnitude of orbit size
      # TODO: use pixels instead of physical coordinates for ticks!
      if self.a <= 1e9:
        tick_len = 1e6
      elif self.a <= 5 * Utils.AU:
        tick_len = 0.03 * Utils.AU
      else:
        print(self.a)
        tick_len = Utils.AU

      xpe, ype, zpe = self.get_periapsis()
      xap, yap, zap = self.get_apoapsis()
      if self.projection == "3D":
        self.ax.scatter([xpe], [ype], [zpe], s=10, color=opts["color"])
        self.ax.scatter([xap], [yap], [zap], s=10, color=opts["color"])
        if show_line_apsides:
          self.ax.plot([xpe,xap], [ype,yap], [zpe,zap], ls=(0,(3,3)), color=opts["color"])
      elif self.projection == "2D":
        self.plot_tick(nu=0, tick_len=tick_len, opts=opts)
        self.plot_tick(nu=pi, tick_len=tick_len, opts=opts)
        # self.ax.scatter([xpe], [ype], s=10, color=opts["color"])
        # self.ax.scatter([xap], [yap], s=10, color=opts["color"])
        if show_line_apsides:
          self.ax.plot([xpe,xap], [ype,yap], ls=(0,(3,3)), color=opts["color"])

    # Show nodes and line of nodes
    if show_nodes:

      r1 = self.radius_at_nu(-self.w)
      xno1 = r1
      yno1 = 0
      xno1, yno1 = rotate2D(xno1, yno1, self.LAN)
      r2 = self.radius_at_nu(180-self.w)
      xno2 = -r2
      yno2 = 0
      xno2, yno2 = rotate2D(xno2, yno2, self.LAN)
      xno1 /= self.units; yno1 /= self.units
      xno2 /= self.units; yno2 /= self.units
      self.ax.plot([xno1], [yno1], "bo")
      self.ax.plot([xno2], [yno2], "bo")
      self.ax.plot([xno1,xno2], [yno1,yno2], "b", ls=(3,(3,3)))

    # Show coordinates axes
    if show_axes:

      self.ax.plot([0,self.a/self.units], [0,0], "r")
      self.ax.plot([0,0], [0,self.a/self.units], "g")
      if self.projection == "3D":
        self.ax.plot([0,0], [0,0], [0,self.a/self.units], "b")

    # Return value: for animation
    return ln, (xs, ys, zs)

  # ----------------------------------------------------------------------------

  # Plots a tick normal to the orbit at the given date
  # Alternatively, if either nu, or both pos & vel, are given, use those
  def plot_tick(self, date=None, nu=None, pos=None, vel=None, tick_len=None, which="both", ax=None, opts=None):

    self._matplotlib_imports()

    if opts is None:
      opts = {}
    if ax is not None:
      self.ax = ax
    if self.ax is None:
      self.init_axes()
    if "color" not in opts:
      opts["color"] = self.color

    valid_which_2D = ["both", "in", "out"]
    valid_which_3d = ["both", "normal", "binormal"]

    which = which.lower()
    if self.projection == "2D" and which not in valid_which_2D:
      raise(Exception("Error in plot_tick(): 'which' must be one of {}".format(valid_which_2D)))
    elif self.projection == "3D" and which not in valid_which_3d:
      raise(Exception("Error in plot_tick(): 'which' must be one of {}".format(valid_which_3d)))

    if date is not None:
      pos, vel = self.posvel_at_date(date)
    elif nu is not None:
      pos, vel = self.posvel_at_nu(nu)
    elif pos is not None and vel is not None:
      pass
    else:
      print("Error in Orbit.plot_tick(): either provide 'date', 'nu', or both 'pos' and 'vel'.")

    hvec = np.cross(pos, vel)
    hnorm = Utils.normalize(hvec)
    nvec = np.cross(vel, hvec)
    normal = Utils.normalize(nvec)

    if tick_len is None:
      r = Utils.magn(pos)
      tick_len = r/10

    if self.projection == "2D":

      if which == "both":
        dp1 = tick_len/2 * normal
        dp2 = tick_len/2 * normal
      elif which == "in":
        dp1 = tick_len/2 * normal
        dp2 = 0
      elif which == "out":
        dp1 = 0
        dp2 = tick_len/2 * normal
      p1 = (pos - dp1)/self.units
      p2 = (pos + dp2)/self.units
      self.ax.plot([p1[0], p2[0]], [p1[1], p2[1]], **opts)

    elif self.projection == "3D":

      p1 = (pos - tick_len/2 * normal)/self.units
      p2 = (pos + tick_len/2 * normal)/self.units
      self.ax.plot([p1[0], p2[0]], [p1[1], p2[1]], [p1[2], p2[2]], **opts)

      if which.lower() in ["both", "binormal"]:
        p3 = (pos - tick_len/2 * hnorm)/self.units
        p4 = (pos + tick_len/2 * hnorm)/self.units
        self.ax.plot([p3[0], p4[0]], [p3[1], p4[1]], [p3[2], p4[2]], **opts)

    if self.projection == "2D":
      return p1, p2

# ==============================================================================
# The following are utility functions designed to be called without directly
# building an Orbit object

# Computes a trajectory from a list of osculating orbital elements
# The list of elements must be *sorted*, with each item of the form
#   (date, a, e, i, arg, LAN, M)
# where date is a tz-naive datetime object, a is in meters and angles
# are in degrees. start_date and end_date must be tz-naive datetime objects.
# max_dist is the maximum distance between successive points; the timestep
# will be adaptively adjusted to remain below that (but not too much below)
# Returned is a list of (date,x,y,z,vx,vy,vz) tuples
def compute_trajectory(primary, elements_list, start_date, end_date, max_dist):
  idx_max = len(elements_list) - 1
  elements, idx = search(elements_list, start_date, key=lambda x: x['epoch'])
  orbit = Orbit(primary, elements)
  pos, vel = orbit.posvel_at_date(elements['epoch'])
  data = [(start_date, pos, vel)]
  date = start_date
  time_step = timedelta(minutes=1)
  thresh_hi = max_dist**2
  thresh_lo = (max_dist/2)**2
  while date <= end_date:
    accept = False
    while not accept:
      new_date = date + time_step
      elements = elements_list[idx]
      # print(new_date, idx, elements['epoch'], elements['e'])
      orbit = Orbit(primary, elements)
      pos, vel = orbit.posvel_at_date(new_date)
      while idx < idx_max and new_date >= elements_list[idx+1]['epoch']:
        idx += 1
      lastpos = data[-1][1]
      last_dist2 = vec_dist_squared(pos, lastpos)
      # print("dist=", sqrt(last_dist2)/1e3)
      if last_dist2 > thresh_hi:
        time_step /= 1.25
        # print("Decreased timestep to ", time_step.total_seconds())
        if time_step.total_seconds() < 1:
          accept = True
      elif last_dist2 < thresh_lo:
        time_step *= 1.25
        # print("Increased timestep to ", time_step.total_seconds())
        accept = True
      else:
        accept = True
    # print("Accepted with dt:", time_step.total_seconds())
    data.append((new_date, pos, vel))
    date = new_date
  return data

# Returns the body-centric inertial position and velocity at date computed
# from a *sorted* list of osculating elements, which must be of the form
#   (date, a, e, i, arg, LAN, M)
# where date is a tz-naive datetime object, a is in meters and angles
# are in degrees.
def posvel_at_date(primary, elements_list, date):
  elset = search(elements_list, date, key=lambda x: x[0])
  o = Orbit(primary, elements=elset)
  return o.posvel_at_date(date)
