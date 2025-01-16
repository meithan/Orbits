# Contains CelestialBody definitions of the bodies of the solar system (currently all planets plus Pluto and Ceres). These can be accessed through the Bodies dict. The show_planets() method (called when this program is run directly) plots the positions of bodies of the solar system at a given date (defaults to now).
import csv
import datetime as dtm
import os
import sys

from CelestialBody import CelestialBody
from Orbit import Orbit
import Utils

# ==============================================================================

AU = Utils.AU
J2000 = Utils.J2000

def parse_float(x):
  try:
    v = float(x)
  except:
    return None
  return v

def sanitize(s):
  if s == "":
    return None
  else:
    return s

# ------------------------------------------------------------------------------

# ALL VALUES GIVEN IN MKS

# Solar system bodies
# Data read from CSV file solar_system_bodies.csv
# The orbital elements of the planets were taken from:
# http://www.met.rdg.ac.uk/~ross/Astronomy/Planets.html
# Do note that this source lists the longitude of the perihelion, ~omega, and
# the mean longitude at epoch, L0, instead of the more commonly used argument
# of periapsis, omega, and mean anomaly at epoch, M0, as is usual. The tradi-
# tional elements can be obtained through the relations:
#   omega = ~omega - LAN
#   M0 = L0 - ~omega
# The epoch for all elements is J2000. These elements should not be used as-is
# outside of the time interval 1800 AD - 2050 AD.
# Ceres elements from JPL Solar System Dynamics page.
package_dir = os.path.dirname(os.path.abspath(__file__))
fpath = os.path.join(package_dir, "solar_system_bodies.csv")
bodies_props = {}
with open(fpath) as f:

  f.readline()
  reader = csv.reader(f)
  for row in reader:

    row = [sanitize(x) for x in row]
    props = {}
    props["name"] = row[0].strip('"')
    props["mu"] = parse_float(row[1])
    props["radius"] = parse_float(row[2])
    props["color"] = row[3]
    props["primary"] = row[4]
    props["a"] = parse_float(row[5])
    props["e"] = parse_float(row[6])
    props["i"] = parse_float(row[7])
    props["L0"] = parse_float(row[8])
    props["~omega"] = parse_float(row[9])
    props["LAN"] = parse_float(row[10])

    bodies_props[props["name"]] = props

# Create bodies first
Bodies = {}
for name in bodies_props:
  Bodies[name] = CelestialBody(name)

# Fill in body properties
for name in bodies_props:
  body = Bodies[name]
  props = bodies_props[name]
  if props["mu"] is not None: body.mu = props["mu"]
  if props["radius"] is not None: body.radius = props["radius"] * 1e3
  if props["color"] is not None: body.color = props["color"]
  if props["primary"] is not None:
    primary = Bodies[props["primary"]]
    if all([props[x] is not None for x in ["a", "e", "i", "L0", "~omega", "LAN"]]):
      body.orbit = Orbit(primary=primary)
      a = props["a"]*AU if primary.name == "Sun" else props["a"]*1e3
      arg_peri = props["~omega"] - props["LAN"]
      M0 = props["L0"] - props["~omega"]
      body.orbit.from_elements({'a': a, 'e': props["e"], 'i': props["i"], 'arg': arg_peri, 'LAN': props["LAN"], 'M0': M0, 'epoch': J2000})

# Add SOIs to planets, in meters
Bodies["Earth"].SOI = 9.24e8
Bodies["Mars"].SOI = 5.76e8
Bodies["Venus"].SOI = 6.16e8

# Separate main objects
Sun = Bodies["Sun"]
Mercury = Bodies["Mercury"]
Venus = Bodies["Venus"]
Earth = Bodies["Earth"]
Moon = Bodies["Moon"]
Mars = Bodies["Mars"]
Ceres = Bodies["Ceres"]
Jupiter = Bodies["Jupiter"]
Saturn = Bodies["Saturn"]
Uranus = Bodies["Uranus"]
Neptune = Bodies["Neptune"]
Pluto = Bodies["Pluto"]

astro_symbols = {
  "Sun": "☉",
  "Mercury": "☿",
  "Venus": "♀",
  # "Earth": "♁",
  "Earth": "⊕",
  "Mars": "♂",
  "Ceres": "⚳",
  "Jupiter": "♃",
  "Saturn": "♄",
  "Uranus": "♅",
  "Neptune": "♆",
  "Pluto": "♇"
}

# ==============================================================================

# Plots the current positions of the planets
def show_planets(datetime=None, projection="2D", bodies=None, dark=True, ax=None, show=True, show_date=True, noaxes=False):

  import datetime as dtm
  import sys
  import matplotlib as mpl
  import matplotlib.pyplot as plt
  from mpl_toolkits.mplot3d import Axes3D

  # -----------------------------------
  # Check arguments

  # Date and time to plot
  if datetime is None:
    datetime = dtm.datetime.now()

  # projection
  projection = projection.upper()
  proj_valid_values = ["2D", "3D"]
  if projection not in proj_valid_values:
    print("'projection' must be one of", proj_valid_values)
    sys.exit()

  # Select celestial bodies to plot
  valid_bodies =  ["Mercury", "Venus", "Earth", "Mars", "Ceres", "Jupiter", "Saturn","Uranus", "Neptune", "Pluto"]
  if bodies is None or bodies == [] or bodies == "":
    bodies = valid_bodies
  elif type(bodies) == str:
    if bodies.capitalize() in valid_bodies:
      bodies = [bodies]
    elif bodies.lower() == "all":
      bodies = valid_bodies
    elif bodies.lower() == "inner":
      bodies = ["Mercury", "Venus", "Earth", "Mars"]
    elif bodies.lower() == "outer":
      bodies = ["Jupiter", "Saturn","Uranus", "Neptune", "Pluto"]
    else:
      print("Error: Unrecognized bodies argument: {}".format(bodies))
      return
  elif type(bodies) in [list, tuple]:
    new_bodies = []
    for b in bodies:
      bb = b.capitalize()
      if bb in valid_bodies:
        new_bodies.append(bb)
      else:
        print("Warning: {} not a valid body name; ignoring".format(b))

  # Use dark style
  if dark:
    plt.style.use('dark_background')
  else:
    plt.style.use('default')

  # ------------------------------------
  # Plot

  if ax is None:
    fig = plt.figure(figsize=(10,10))
    if projection == "3D":
      ax = fig.add_subplot(111, projection='3d')
    elif projection == "2D":
      ax = fig.add_subplot(111)
      ax.set_aspect('equal')

  # Bodies
  for i,name in enumerate(bodies):
    body = Bodies[name]
    # show_axes = True if i == 0 else False
    show_axes = False
    body.plot_orbit(ax=ax, projection=projection, show_axes=show_axes, units=AU, show_apsides=True, dark=True)
    if body.name in astro_symbols:
      label = astro_symbols[body.name]
    else:
      label = body.name
    body.plot_at_date(datetime, opts=dict(s=50), label=label, label_opts=dict(fontsize=18))
    # TEMP: custom labels
    # names = {"Mercury": "Mercurio", "Venus": "Venus", "Earth": "Tierra", "Mars": "Marte"}
    # dxdy = {"Mercury": (20,0), "Venus": (-75,-25), "Earth": (-50,-30), "Mars": (20,0)}
    # label = names[body.name]
    # label = body.name
    # dx, dy = dxdy[body.name]
    # body.plot_at_date(datetime, opts=dict(s=200))
    # pos, vel = body.posvel_at_date(datetime)
    # pos = pos/Utils.AU
    # print(label, pos)
    # plt.annotate(label, xy=(pos[0], pos[1]), xytext=(dx, dy), textcoords="offset pixels", color=body.color, fontsize=18, va="center")

  if projection == "2D":
    max_rap = None
    for body in bodies:
      rap = Bodies[body].orbit.rap/Utils.AU
      if max_rap is None or rap > max_rap:
        max_rap = rap
    L = max_rap*1.05
    plt.xlim(-L, L)
    plt.ylim(-L, L)

  # for name in v:
  #   planet = Bodies[name]
  #   print("%s %.3f x %.3f AU, %.1f°" % (name.ljust(8), planet.orbit.rpe/AU, planet.orbit.rap/AU, planet.orbit.i))

  # Plot Sun
  sun_color = "yellow"
  if projection == "2D":
    plt.scatter([0], [0], color=sun_color, zorder=10)
    plt.annotate(astro_symbols["Sun"], xy=(0,0), xytext=(0, -25), textcoords="offset pixels", ha="center", va="center", fontsize=18, color=sun_color)
  elif projection == "3D":
    ax.scatter([0], [0], [0], color=sun_color, zorder=10)
  
  # Shown date
  if show_date:
    date_str = datetime.strftime("%Y-%m-%d %H:%M %z")
    if projection == "2D":
      plt.text(0.01, 0.99, date_str, transform=ax.transAxes, ha="left", va="top")
    elif projection == "3D":
      pass

  # Grid & panes
  if projection == "2D":
    plt.grid(ls="-", color="0.1" if dark else "0.9")
  elif projection == "3D":
    pane_color = (0, 0, 0, 0)
    ax.xaxis.set_pane_color(pane_color)
    ax.yaxis.set_pane_color(pane_color)
    ax.zaxis.set_pane_color(pane_color)

  # Axis labels
  plt.xlabel("AU")
  plt.ylabel("AU")

  # Hide axes
  if noaxes:
    ax.grid(False)
    plt.subplots_adjust(top=0.99, bottom=0.05, left=0.075, right=0.975, hspace=0.2, wspace=0.2)
    plt.axis("off")

  plt.tight_layout()

  if show:
    plt.show()

# ==============================================================================

if __name__ == "__main__":

  import argparse
  import datetime as dtm

  parser = argparse.ArgumentParser(
    description = "Plots the positions of the solar system planets.",
  )
  parser.add_argument("datetime", nargs='?', default=None, help="An UTC date, and possibly time, in YYYY-MM-DD or YYYY-MM-DD HH:MM:SS format. Defaults to current UTC date and time.")
  parser.add_argument("--inner", action="store_true", help="Show only inner planets")
  parser.add_argument("--outer", action="store_true", help="Show only outer planets")
  parser.add_argument("--all", action="store_true", help="Show all major bodies (the default if no flag is passed).")
  parser.add_argument("--light", action="store_true", default=False, help="Plot with a light background instead of the default dark backgorund.")
  parser.add_argument("--projection", nargs=1, action="store", help="Specify the projectionection for the plot: either '2D' or '3D'. Default is '2D'.")
  parser.add_argument("--noaxes", action="store_true", help="Hide axes and grid lines.")

  args = parser.parse_args()

  if args.datetime is None or args.datetime.lower() == "now":
    # datetime = dtm.datetime.now()
    datetime = dtm.datetime.now(dtm.timezone.utc).replace(tzinfo=None)
  else:
    if len(args.datetime.strip()) == 10:
      datetime = dtm.datetime.strptime(args.datetime, "%Y-%m-%d")
    else:
      datetime = dtm.datetime.strptime(args.datetime, "%Y-%m-%d %H:%M:%S")

  if args.inner:
    bodies = "inner"
  elif args.outer:
    bodies = "outer"
  else:
    bodies = "all"

  if args.projection is not None:
    projection = args.projection[0].lower()
  else:
    projection = "2d"

  dark = not args.light

  # bodies = ["Mercury", "Venus", "Earth", "Mars", "Jupiter", "Uranus"]

  show_planets(datetime=datetime, bodies=bodies, dark=dark, projection=projection, noaxes=args.noaxes)

  # HACK -- Output state vectors
  # dt = dtm.datetime.utcnow()
  # for name in ["Mercury", "Venus", "Earth", "Mars", "Ceres", "Jupiter", "Saturn","Uranus", "Neptune", "Pluto"]:
  #   body = Bodies[name]
  #   pos, vel = body.posvel_at_date(dt)
  #   # print(name, pos, vel)

