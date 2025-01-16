# Functions to compute transfer between celestial bodies.
#
# Main functions are:
# compute_porkchop(): creates (and plots) a porkchop plot,
# pykep_lambert_transfer(): computes a single Lambert solution transfer
#
# This program relies heavily on PyKEP.

import copy
import datetime as dtm
from math import sqrt, asin, pi, degrees
import os
import sys

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import matplotlib.ticker as mticker
import numpy as np

import pykep
import pykep.planet

import Utils

# ==============================================================================

# Converts a Python datetime to a pykep epoch object
def dt2pykep(dt):
  return pykep.epoch_from_string(dt.strftime("%Y-%m-%d %H:%M:%S"))

# J2000 epoch as a pykep epoch
pykep_J2000 = pykep.epoch_from_string('2000-01-01 00:00:00')

def nice_tof(t):
  if t > 500:
    return t/365.24, "yr", "years"
  else:
    return t, "d", "days"

# ==============================================================================

def compute_porkchop(orig_body, dest_body, departures, tofs=None, arrivals=None, dep_alt=200, arr_alt=200, plot_var="dep_C3", vmin=None, vmax=None, plot=True, save_fname=None):
  """ Computes a porkchop plot by repeatedly solving Lambert's problem.
  Arguments:
    orig_body: origin celestial body, a pykep planet
    dest_body: destination celesial body, a pykep planet
    departures: a (min_date, max_date, num_dates) tuple with the range and
      number of departure dates to compute.
    tofs: a (min_tof, max_tof, num_tofs) tuple with the range of times-of-flight
      to compute. Either this or 'arrivals' must be provided.
    arrivals: a (min_date, max_date, num_dates) tuple with the range and number
      of departure dates to compute. Either this or 'tofs' must be provided.
    dep_alt: altitude of departure circular parking orbit, in km.
    arr_alt: altitude of arrival circular orbit, in km.
    plot_var: variable to plot in the porkchop; see valid_plot_vars for
      recognized options.
    vmin, vmax: minimum and maximum values for the colormap; set to None for
      automatic limits. If vmax is given, values above this will be set to gray.
    plot: set to False to skip plot, an only return the arrays.
    save_fname: if given, will save the plot to this filename or filepath.
    mu: the gravitational parameter of the primary, in m^3/s^2; defaults to Sun's.
    All dates must be timezone-naive datetime objects with UTC time.
  Returns:
    An dict with the following keys containinng numpy arrrays, all of shape
    (num_deps x num_deps or num_tofs):
      dep_C3: launch characteristic energy (C3, km²/s²)
      dep_vinf: launch hyperbolic excess speed (vinf, km/s)
      ejec_dv: ejection Δv from initial parking orbit (km/s)
      arr_C3: arrival C3 (km²/s²)
      arr_vinf: arrival vinf (km/s)
      cap_dv: minimum capture Δv (km/s)
      low_dv: low-orbit insertion Δv (km/s)
  """

  # Variables to be computed (and possibly plotted)
  variables = {
    "dep_C3": dict(name="Launch C3", short_name="Dep C3", units="km²/s²"),
    "arr_C3": dict(name="Arrival C3", short_name="Arr C3", units="km²/s²"),
    "dep_vinf": dict(name="Departure vinf", short_name="Dep vinf", units="km/s"),
    "arr_vinf": dict(name="Arrival vinf", short_name="Arr vinf", units="km/s"),
    "ejec_dv": dict(name="Ejection Δv", short_name="Ejec Δv", units="km/s"),
    "cap_dv": dict(name="Capture Δv", short_name="Cap Δv", units="km/s"),
    "low_dv": dict(name="Low-orbit Δv", short_name="Low Δv", units="km/s"),
    "tot_dv": dict(name="Total Δv", short_name="Tot Δv", units="km/s"),
  }

  # Check plot_var
  if plot_var not in variables:
    print("Error: plot_var must be one of:")
    print(variables.keys())
    return None

  orig_name = orig_body.name.title()
  dest_name = dest_body.name.title()

  # Range of dates (by mjd)
  min_dep_date, max_dep_date, num_deps = departures
  min_dep_mjd = dt2pykep(min_dep_date).mjd
  max_dep_mjd = dt2pykep(max_dep_date).mjd
  dep_mjds = np.linspace(min_dep_mjd, max_dep_mjd, num_deps)
  dep_dates = [Utils.mjd2dt(x) for x in dep_mjds]

  # Range of tofs
  if tofs is not None:
    mode = "tof"
    min_tof, max_tof, num_tofs = tofs
    tofs = np.linspace(min_tof, max_tof, num_tofs)
    num_ys = num_tofs
    ys = tofs
    print("y-axis is time-of-flight")

  # Range of arrival dates
  elif arrivals is not None:
    mode = "arrival"
    min_arr_date, max_arr_date, num_arrs = arrivals
    min_arr_mjd = dt2pykep(min_arr_date).mjd
    max_arr_mjd = dt2pykep(max_arr_date).mjd
    arr_mjds = np.linspace(min_arr_mjd, max_arr_mjd, num_arrs)
    arr_dates = [Utils.mjd2dt(x) for x in arr_mjds]
    num_ys = num_arrs
    ys = arr_mjds
    print("y-axis is arrivale date")

  else:
    print("You must provide either tofs=(min_tof, max_tof, num_tofs) --or-- arr_dates=(min_arr, max_arr, num_arrs)")
    return None

  # Allocate data arrays
  result = {
    "dep_C3": np.zeros((num_deps, num_ys)),
    "dep_vinf": np.zeros((num_deps, num_ys)),
    "ejec_dv": np.zeros((num_deps, num_ys)),
    "arr_C3": np.zeros((num_deps, num_ys)),
    "arr_vinf": np.zeros((num_deps, num_ys)),
    "cap_dv": np.zeros((num_deps, num_ys)),
    "low_dv": np.zeros((num_deps, num_ys)),
    "tot_dv": np.zeros((num_deps, num_ys)),
  }

  # The missing variable (tof or arr_date) can be computed from others,
  # so it is not saved
  result["dep_dates"] = dep_dates
  if mode == "tof":
    result["tofs"] = tofs
  elif mode == "arrival":
    result["arr_dates"] = arr_dates

  # Best transfer for each variable
  result["best"] = {}
  for var in result.keys():
    result["best"][var] = None

  print("\nComputing porkchop plot ...")
  for i, dep_mjd in enumerate(dep_mjds):

    t1 = pykep.epoch(dep_mjd, "mjd")

    for j, y in enumerate(ys):

      if mode == "tof":
        tof = y
        arr_mjd = dep_mjd + tof*pykep.DAY2SEC
      elif mode == "arrival":
        arr_mjd = y
        tof = arr_mjd - dep_mjd

      t2 = pykep.epoch(t1.mjd2000 + tof)
      ro1, vo1 = orig_body.eph(t1)
      rd1, vd1 = dest_body.eph(t1)
      ro2, vo2 = orig_body.eph(t2)
      rd2, vd2 = dest_body.eph(t2)

      # Lambert solver
      vinf1, vinf2, DLA, ejec_dv, cap_dv, low_dv, lam = pykep_lambert_transfer(orig_body, dest_body, t1, tof, dep_alt, arr_alt)

      result["dep_vinf"][i,j] = vinf1/1e3
      result["arr_vinf"][i,j] = vinf2/1e3
      result["dep_C3"][i,j] = vinf1**2/1e6
      result["arr_C3"][i,j] = vinf2**2/1e6
      result["ejec_dv"][i,j] = ejec_dv/1e3
      result["cap_dv"][i,j] = cap_dv/1e3
      result["low_dv"][i,j] = low_dv/1e3
      result["tot_dv"][i,j] = (ejec_dv+cap_dv)/1e3

      for var in variables.keys():
        if result["best"][var] is None or result[var][i,j] < result["best"][var][0]:
          result["best"][var] = (result[var][i,j], (i,j))

    if i == 0 or (i+1) % (num_deps//10) == 0:
      print(Utils.mjd2dt(dep_mjd).strftime("%Y-%m-%d"), "({}/{})".format(i+1, num_deps))

  var_name = variables[plot_var]["name"]
  var_short_name = variables[plot_var]["short_name"]
  var_units = variables[plot_var]["units"]

  print(f"\nBest transfer for {var_name}:")
  value = result["best"][plot_var][0]
  i,j = result["best"][plot_var][1]
  dep_mjd = dep_mjds[i]
  dep_dt = Utils.mjd2dt(dep_mjd)
  if mode == "tof":
    tof = tofs[j]
    arr_mjd = dep_mjds[j] + tof
  elif mode == "arrival":
    arr_mjd = arr_mjds[j]
    tof = arr_mjds[j] - dep_mjds[i]
  arr_dt = Utils.mjd2dt(arr_mjd)
  t, u, uu = nice_tof(tof)
  dep_C3 = result["dep_C3"][i,j]
  arr_C3 = result["arr_C3"][i,j]
  ejec_dv = result["ejec_dv"][i,j]
  cap_dv = result["cap_dv"][i,j]
  low_dv = result["low_dv"][i,j]
  tot_dv = result["tot_dv"][i,j]
  print(f"{var_name}: {value:.1f} {var_units}")
  print(f"Departing {dep_dt.strftime('%Y-%m-%d')}")
  print(f"Arriving {arr_dt.strftime('%Y-%m-%d')}")
  print(f"Time of flight: {t:.1f} {uu}")
  print(f"Departure C3: {dep_C3:.1f} km²/s²")
  print(f"Arrival C3: {arr_C3:.1f} km²/s²")
  print(f"Departure Δv: {ejec_dv:.1f} km/s")
  print(f"Capture Δv: {cap_dv:.1f} km/s")
  print(f"Low-orbit Δv: {low_dv:.1f} km/s")
  print(f"Total Δv: {tot_dv:.1f} km/s")

  # ------------------------------------
  # Plot porkchop
  if plot:

    data = result[plot_var]

    if vmin is None: vmin = data.min()
    if vmax is None:
      if var_units == "km/s":
        vmax = 10
      elif var_units == "km²/s²":
        vmax = 200

    xlabel = "Departure date"

    if mode == "tof":
      t, tof_units, tof_units_label = nice_tof(min(tofs))
      if tof_units == "yr":
        tofs /= 365.24
      ylabel = f"Time of flight [{tof_units_label}]"
    elif mode == "arrival":
      ylabel = "Arrival date"
      tof_units = "d"

    fig = plt.figure(figsize=(7.5,6))
    ax = plt.gca()

    # Compute departure dates
    if mode == "arrival":
      ys = arr_dates

    # data[data > vmax] = np.nan

    cmap = copy.copy(cm.get_cmap("jet"))
    cmap.set_over("gray")

    # im = plt.pcolormesh(dep_mjds, tofs, data, cmap=cmap, norm=mcolors.LogNorm(vmin=vmin, vmax=vmax), shading="flat")
    im = plt.pcolormesh(dep_dates, ys, data.T, cmap=cmap, vmin=vmin, vmax=vmax, shading="auto")

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    cb = plt.colorbar(im, fraction=0.046, pad=0.02)
    cb.set_label(f"{var_name} [{var_units}]")

    title = f"{orig_name} → {dest_name} transfer (direct) • {var_name}"
    plt.title(title)

    plt.gca().xaxis.set_major_locator(mdates.MonthLocator())

    # cts = plt.contour(dep_mjds, tofs, data, colors="black", levels=[6, 7, 8, 9, 10, 15, 20, 30, 40, 50], alpha=1.0)
    # plt.clabel(cts, colors="black", fmt='%i')

    def format_coord(x, y):
      min_x, max_x = min(dep_dates), max(dep_dates)
      min_y, max_y = min(ys), max(ys)
      i = int((mdates.num2date(x).replace(tzinfo=None)-min_x)/(max_x-min_x)*(len(dep_dates)-1))
      if mode == "tof":
        tof = y
        arr_dt = mdates.num2date(x) + dtm.timedelta(days=y)
        j = int((y-min_y)/(max_y-min_y)*(len(ys)-1))
      elif mode == "arrival":
        tof = y - x
        arr_dt = mdates.num2date(y)
        j = int((mdates.num2date(y).replace(tzinfo=None)-min_y)/(max_y-min_y)*(len(ys)-1))
      dep_str = mdates.num2date(x).strftime('%Y-%m-%d')
      arr_str = arr_dt.strftime('%Y-%m-%d')
      if mode == "arrival":
        tof, tof_units, _ = nice_tof(tof)
      tof_str = f"{tof:.1f} {tof_units}"
      return f"Dep: {dep_str}, Arr: {arr_str}, ToF: {tof_str}, {var_short_name}: {data[i,j]:.1f} {var_units}"
    plt.gca().format_coord = format_coord

    plt.tight_layout()

    if save_fname is not None:
      plt.savefig(save_fname)
      print("Saved", save_fname)
    plt.show()

  # ------------------------------------

  return result

# ==============================================================================

def pykep_lambert_transfer(orig, dest, t1, tof, h1, h2=None):
  """Computes a direct Hohmann-type transfer departing from a circular

  Arguments:
    orig: the origin, a pykep.planet instance
    dest: the destination, a pykep.planet instance
    t1: the time of departure (a Python datetime or a pykep time)
    tof: time of flight, in days
    h1: altitude of departing circular orbit, km
    h2: altitude of arrival circular orbit, km (optional)

  Returns:
    The two hyperbolic excess speeds (vinfs), the declination of the launch asymptote, the ejection and minimum capture Δv's, and the lambert object returned by pykep (which can be used to generate the trajectory)

  Note: the delta-v's assume the launch/arrival hyperbolas are coplanar with the initial/final orbits.

  Requires PyKEP
  """

  if orig.mu_central_body != dest.mu_central_body:
    print("Planets {} and {} have differrent mu_central_body, but they must orbit the same primary.".format(orig.name, dest.name))
    return None
  mu = orig.mu_central_body

  if type(t1) in [dtm.datetime, dtm.date]:
    _t1 = dt2pykep(t1)
  else:
    _t1 = t1

  _t2 = pykep.epoch(_t1.mjd2000 + tof)
  rO1,vO1 = orig.eph(_t1)
  rO2,vO2 = orig.eph(_t2)
  rD1,vD1 = dest.eph(_t1)
  rD2,vD2 = dest.eph(_t2)

  if h2 is None:
    h2 = h1

  # Lambert solver
  lam = pykep.lambert_problem(r1=rO1, r2=rD2, tof=tof*pykep.DAY2SEC, mu=mu, cw=False, max_revs=0)

  # # Extract lowest-energy solution from Lambert solver
  # vels1 = lam.get_v1()
  # vels2 = lam.get_v2()
  # lowest = None
  # best_i = None
  # for i in range(len(lam.get_v1())):
  #   v1 = vels1[i]
  #   v1mag = magn(v1)
  #   if lowest is None or v1mag < lowest:
  #     lowest = v1mag
  #     best_i = i
  # v1 = lam.get_v1()[best_i]
  # v2 = lam.get_v2()[best_i]

  # Extract solution from Lambert solver (0 revs case)
  v1 = lam.get_v1()[0]
  v2 = lam.get_v2()[0]
  vinf1_vec = Utils.vec_diff(v1, vO1)
  vinf2_vec = Utils.vec_diff(vD2, v2)
  vinf1 = Utils.magn(vinf1_vec)
  vinf2 = Utils.magn(vinf2_vec)

  if not hasattr(orig, "SOI"):
    orig.SOI = float("inf")
  if not hasattr(dest, "SOI"):
    dest.SOI = float("inf")

  # Compute departure quantities
  # The deltav's assume that the launch and arrival hyperbolas are coplanar 
  # with the initial/final orbits, which may not necessarily be true.
  E1 = 0.5*vinf1**2 - orig.mu_self/orig.SOI
  r1 = orig.radius + h1*1e3
  vper1 = sqrt(2*(E1+orig.mu_self/r1))
  vcirc1 = Utils.circ_speed(orig, h1*1e3)
  ejec_dv = vper1 - vcirc1
  DLA = degrees(asin(vinf1_vec[2]/vinf1))

  # Compute arrival quantities
  # Same comment above applies
  E2 = 0.5*vinf2**2 - dest.mu_self/dest.SOI
  r2 = dest.radius + h2*1e3
  vper2 = sqrt(2*(E2+dest.mu_self/r2))
  vcap2 = sqrt(2*dest.mu_self/r2)    # Speed at periapsis for E=0
  vcirc2 = Utils.circ_speed(dest, h2*1e3)
  cap_dv = vper2 - vcap2
  low_dv = vper2 - vcirc2

  return vinf1, vinf2, DLA, ejec_dv, cap_dv, low_dv, lam

# ==============================================================================

# Generates the trajectory for the given pykep lambert problem solution
# Returned is a list of (dt, r, v) tuples, where dt is the total time
# elapsed since the initial position, and r and v are the 3D position
# and velocity vectors
def pykep_lambert_trajectory(lam, N=200):

  r = lam.get_r1()
  v = lam.get_v1()[0]
  T = lam.get_tof()
  mu = lam.get_mu()
  dt = T / (N - 1)

  traj = [(0, r, v)]
  for i in range(1,N):
    r, v = pykep.propagate_lagrangian(r, v, dt, mu)
    dt_tot = dt*i
    traj.append((dt_tot,r,v))

  return traj

# ==============================================================================

# Generates the trajectory for the given pykep planet object between
# the given dates (python datetimes)
# Returned is a list of (t, r, v) tuples, where t are Python datetime
# objects and r and v are the 3D position  and velocity vectors
def pykep_planet_trajectory(planet, date1, date2, N=200):

  t0 = dt2mjd(date1)
  T = (date2-date1)/dtm.timedelta(days=1)
  dt = T / (N - 1)

  r, v = planet.eph(pykep.epoch(t0, "mjd"))
  traj = [(date1, r, v)]
  for i in range(1,N):
    r, v = planet.eph(pykep.epoch(t0 + i*dt, "mjd"))
    t = date1 + dtm.timedelta(days=i*dt)
    traj.append((t, r, v))

  return traj

# ==============================================================================

# Returns the position and velocity vectors of planet at given date
def pykep_planet_posvel(planet, date, units=None):
  t0 = dt2mjd(date)
  r, v = planet.eph(pykep.epoch(t0, "mjd"))
  if units is not None:
    r = tuple(x/units for x in r)
    v = tuple(x/units for x in v)
  return r, v

# ==============================================================================

if __name__ == "__main__":

  # Example: compute Earth to Mars transfer
  import numpy as np

  # Specify departure date and time of flight (days)
  departure_date = dtm.datetime(2022, 9, 8)
  time_of_flight = 200

  Earth = pykep.planet.jpl_lp("earth")
  Mars = pykep.planet.jpl_lp("mars")

  # # Add orbital periods, in seconds
  # Earth.P = Earth.compute_period(Utils.pykep_J2000)
  # Mars.P = Mars.compute_period(Utils.pykep_J2000)

  # Add SOIs
  Earth.SOI = 9.24e8
  Mars.SOI = 5.76e8

  lam, vinf1, vinf2, DLA, deltav1, deltav2, deltavcap = pykep_lambert_transfer(Earth, Mars, departure_date, time_of_flight, 200, 200)

  dep_C3 = (np.linalg.norm(vinf1)/1e3)**2
  print(departure_date.date(), dep_C3)

  # Find best transfer given that time_of_flight and dates range
  # start_date = dtm.datetime(2022, 9, 1)
  # end_date = dtm.datetime(2022, 10, 15)
  # departure_date = start_date
  # while departure_date <= end_date:
  #
  #   lam, vinf1, vinf2, DLA, deltav1, deltav2, deltavcap = pykep_lambert_transfer(Earth, Mars, departure_date, time_of_flight, 200, 200)
  #
  #   dep_C3 = (np.linalg.norm(vinf1)/1e3)**2
  #   print(departure_date.date(), dep_C3)
  #
  #   departure_date += dtm.timedelta(days=1)
