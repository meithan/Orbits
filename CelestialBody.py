# A class to represent celestial bodies, with some wrappers into the Orbit class

import Utils

# A celestial body
class CelestialBody:

  def __init__(self, name, mu=None, radius=None, rot_period=None, orbit=None):
    self.name = name
    self.mu = mu
    if mu is not None:
      self.mass = Utils.GRAV/self.mu
    self.radius = radius
    self.rot_period = rot_period
    self.orbit = orbit
    self.color = None

  def get_elements(self):
    return self.orbit.get_elements()

  def plot_at_date(self, date_time, ax=None, projection=None, units=None,
  opts=None, label=None, label_opts=None):
    if self.color is not None:
      if opts is None: opts = {}
      opts["color"] = self.color
    sc = self.orbit.plot_at_date(date_time, ax=ax, projection=projection, units=units,
    opts=opts, label=label, label_opts=label_opts)
    return sc

  def plot_orbit(self, projection=None, ax=None, show_apsides=False,  show_line_apsides=False, show_nodes=False, show_axes=False, show_primary=False, units=None, dark=True, opts=None):
    if self.color is not None:
      if opts is None: opts = {}
      opts["color"] = self.color
    ln = self.orbit.plot(projection=projection, ax=ax, show_apsides=show_apsides,
    show_nodes=show_nodes, show_axes=show_axes, show_primary=show_primary,
    opts=opts, units=units)
    return ln

  def posvel_at_date(self, date_time, dist_units=None, vel_units=None):
    return self.orbit.posvel_at_date(date_time, dist_units=dist_units, vel_units=vel_units)

  def __repr__(self):
    return "<CelestialBody '%s'>" % self.name
