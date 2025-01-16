# Module to compute positions, ground tracks and passes of 
# Earth Satellites from their TLEs

# Uses Skyfield

import datetime as dtm
from zoneinfo import ZoneInfo
from math import sqrt, sin, cos, radians, degrees, asin, pi
import requests
import sys

# Requirements: skyfield
# For plotting groundtracks: matplotlib and cartopy
import skyfield.api as sf

UTC = ZoneInfo("UTC")

# ==============================================================================

# Download latest TLEs from Celestrak
# norad_cat_ids = ["40697", "42063"]
# for norad_cat_id in norad_cat_ids:
#   celestrak_url = "https://www.celestrak.com/satcat/tle.php?CATNR={}".format(norad_cat_id)
#   print(celestrak_url)
#   r = requests.get(celestrak_url)
#   if r.status_code != 200:
#     print("Error retrieving norad_cat_id {}".format(norad_cat_id))
#     sys.exit()
#   lines = r.text.split("\n")
#   sat = EarthSatellite(lines[1], lines[2], lines[0], ts)
#   satellites.append(sat)

# ==============================================================================

# Loads a satellite from a TLE (either a file or the two/three lines)
# Accepts both 2-line and 3-line TLEs
def from_TLE(fname=None, line1=None, line2=None, name=None):
  
  if line1 is not None and line2 is not None:
    sat = sf.EarthSatellite(line1=line1, line2=line2, name=name)
  
  elif fname is not None:
    with open(fname) as f:
      lines = f.readlines()
    if len(lines) == 2:
      sat = sf.EarthSatellite(line1=lines[0], line2=lines[1])
    elif len(lines) == 3:
      sat = sf.EarthSatellite(line1=lines[1], line2=lines[2], name=lines[0])
    else:
      print("File must contain a single 2-line or 3-line TLE")
      return None
  else:
    print("Satellite.from_TLE error: Either supply fname or line1 & line2")
    return None
  
  sat.period = 60 * 2*pi / sat.model.no
  
  return sat

# ==============================================================================

def subpoint_now(sat):
  return subpoint_at_dt(sat, datetime.utcnow())

def subpoint_at_dt(sat, dt):
  ts = sf.load.timescale(builtin=True)
  t = ts.from_datetime(dt)
  geocentric = sat.at(t)
  subpoint = geocentric.subpoint()
  lat, lon = sf.wgs84.latlon_of(geocentric)
  height = sf.wgs84.height_of(geocentric)
  return lat.degrees, lon.degrees, height.km

# ==============================================================================

# Computes successive subpoints between dt1 and dt2
# Step in seconds
def calc_ground_track(sat, dt1, dt2, step=1):
  dt = dt1
  track = []
  while dt <= dt2:
    lat, lon, alt = subpoint_at_dt(sat, dt)
    track.append((dt, lat, lon))
    dt += dtm.timedelta(seconds=step)
  return track

# ==============================================================================

# Plots the ground track of a satellite on a world map
# Will plot revs orbits (default=1) centered on the TLE epoch if t1 and t2
# are None; if only t1 is given, will plot revs orbits centered on t1; and
# if both t1 and t2 are given, will plot between these two times.
def plot_ground_track(sat, revs=1, t1=None, t2=None, show_at_epoch=False, proj=None, ax=None, plot_opts={}, step=1):

  import matplotlib.pyplot as plt
  import cartopy.crs as ccrs
  from cartopy.feature.nightshade import Nightshade
  import cartopy.feature as cf

  epoch_dt = sat.epoch.utc_datetime()
  
  if t1 is None and t2 is None:
    # Plot revs orbits centered on epoch
    delta = dtm.timedelta(seconds=revs*sat.period/2)
    t1 = epoch_dt - delta
    t2 = epoch_dt + delta
  if t1 is not None and t2 is None:
    # Plot revs orbits centered on t1
    tc = t1
    delta = dtm.timedelta(seconds=revs*sat.period/2)
    t1 = tc - delta
    t2 = tc + delta

  if proj is None:
    proj = ccrs.Robinson()
    # proj = ccrs.PlateCarree()

  if ax is None:
    plt.figure(figsize=(10,5))
    ax = plt.axes(projection=proj)

  if "color" not in plot_opts:
    plot_opts["color"] = "k"

  # Plot ground track
  track = calc_ground_track(sat, t1, t2, step=step)
  dts, lats, lngs = zip(*track)
  plt.plot(lngs, lats, transform=ccrs.Geodetic(), **plot_opts)

  # Show position at epoch if requested
  if show_at_epoch:
    lat, lng, _ = subpoint_at_dt(sat, epoch_dt)
    plt.scatter([lng], [lat], transform=ccrs.Geodetic(), color="k", s=20)

  # ax.set_global()
  # ax.stock_img()
  # ax.gridlines(ls=":")
  # ax.set_axisbelow(True)
  # ax.add_feature(Nightshade(epoch_dt, alpha=0.2))
  # ax.add_feature(cf.COASTLINE, edgecolor='gray', lw=0.5)
  # ax.add_feature(cf.BORDERS, edgecolor='gray', lw=0.5)

  # plt.title("{} • Epoch: {} UTC".format(sat.name, epoch_dt.strftime("%Y-%m-%d %H:%M:%S")))

  # plt.tight_layout()
  # plt.show()

# # Compute minimum angle above horizon for satellite to be within swath
# # Width (full) of the swath
# def min_altitude_swath(swath, orbit_alt):

#   beta = (swath/2)/R
#   Rph = R + orbit_alt
#   a = sqrt(R**2 + Rph**2 - 2*R*Rph*cos(beta))
#   alpha = pi - asin(Rph/a * sin(beta))
#   theta = alpha - pi/2
#   return degrees(theta)

# ==============================================================================

# Find satellite passes near a geographic location
# 'only' allows to filter passes by day/night based on situation at
# pass culmination, and must be one of the following strings:
#   "day", "daylight": pass occurs during the day at location
#   "night": pass occurs during the night at location
#   "twilight": pass occurs during twilight
#   "visible": pass occurs during twilight or night, and the satellite is lit
def find_passes(satellite, latitude, longitude, min_datetime, max_datetime, min_elevation=10, min_distance=None, only=None, direction="both", verbose=False):

  import cartopy.geodesic

  direction = direction.lower()

  location = sf.Topos(latitude_degrees=latitude, longitude_degrees=longitude)
  
  if min_datetime.tzinfo is None:
    min_datetime = min_datetime.replace(tzinfo=UTC)
  if max_datetime.tzinfo is None:
   max_datetime = max_datetime.replace(tzinfo=UTC)

  ts = sf.load.timescale(builtin=True)
  t1 = ts.from_datetime(min_datetime)
  t2 = ts.from_datetime(max_datetime)
  
  if only is not None:
    import skyfield.almanac as almanac
    eph = sf.load('de421.bsp')
    f = almanac.dark_twilight_day(eph, location)

  if verbose:
    print("\n>> {}".format(satellite.name))

  tts, evts = satellite.find_events(location, t1, t2, altitude_degrees=min_elevation)

  passes = []
  for i in range(len(tts)-2):
    if evts[i] == 0 and evts[i+1] == 1 and evts[i+2] == 2:

      tt_start = tts[i]
      tt_culm = tts[i+1]
      tt_end = tts[i+2]

      dt_start = tt_start.utc_datetime()
      dt_culm = tt_culm.utc_datetime()
      dt_end = tt_end.utc_datetime()

      if direction in ["ascending", "asc", "descending", "desc"]:
        
        lat1, _, _ = subpoint_at_dt(satellite, dt_start)
        lat2, _, _ = subpoint_at_dt(satellite, dt_end)

        if direction in ["ascending", "asc"]:
          if lat2 < lat1:
            continue
        elif direction in ["descending", "desc"]:
          if lat2 > lat1:
            continue

      # Filter by day/night and visibility
      if only is not None:   
        only = only.lower()
        code = int(f(tt_culm))
        if only in ["day", "daylight", "daytime"] and code != 4:
          continue
        elif only in ["night", "nighttime"] and code != 0:
          continue
        elif only in ["twilight"] and code not in [1, 2, 3]:
          continue
        elif only in ["visible"]:
          if code == 4:
            continue
          sunlit = satellite.at(tt_culm).is_sunlit(eph)
          if not sunlit:
            continue

      # Filter by minimum ground track distance
      # Currently just checks (ground) distance at culmination
      if min_distance is not None:
        lat1, lon1, _ = subpoint_at_dt(satellite, dt_culm)
        geoid = cartopy.geodesic.Geodesic()
        result = geoid.inverse((lon1, lat1), (longitude, latitude))
        min_dist = result[0][0]/1e3
        if min_dist > min_distance:
          continue
      
      passes.append((dt_start, dt_culm, dt_end))
      
      if verbose:
        # lat1, lon1, alt1 = subpoint_at_dt(satellite, dt_start)
        # lat2, lon2, alt2 = subpoint_at_dt(satellite, dt_end)        
        # print(dt_culm.strftime("%Y-%m-%d %H:%M:%S"), lat1, lon1, lat2, lon2)
        print("{}, {}, {}, {:.1f} km".format(dt_start.strftime("%Y-%m-%d %H:%M:%S"), dt_culm.strftime("%Y-%m-%d %H:%M:%S"), dt_end.strftime("%Y-%m-%d %H:%M:%S"), min_dist/1e3))
  
  if verbose:
    if len(passes) == 0:
      print("No passes found in given time range")
    else:
      print("{} pass{} found".format(len(passes), "es" if len(passes) > 1 else ""))
  
  return passes

# ==============================================================================

# Plots the groundtrack of a satellite pass
def plot_pass(sat, dt_start, dt_culm, dt_end, ax=None, crs=None, step=10, ticks_interval=None, plot_opts={}, tzone=None):

  import cartopy
  import cartopy.crs as ccrs
  import matplotlib.pyplot as plt  
  
  PlateCarree = ccrs.PlateCarree()
  Geodetic = ccrs.Geodetic()

  if crs is None:
    crs = PlateCarree

  if ax is None:
    plt.figure(figsize=(12,8))
    ax = plt.axes(projection=crs)
    ax.coastlines(resolution='50m')
    ax.add_feature(cartopy.feature.OCEAN)
    ax.add_feature(cartopy.feature.LAND, edgecolor='0.3', facecolor="tan")
    extent = [-106, -92, 15, 25]
    ax.set_extent(extent, crs=PlateCarree)
    # ax.set_global()
    ax.gridlines(linestyle=":", color="0.5", draw_labels=True)

  # Compute ground track
  track = calc_ground_track(sat, dt_start, dt_end, step)
  dts, lats, lons = zip(*track)

  # t_start = ts.from_datetime(dt_start)
  # t_culm = ts.from_datetime(dt_culm)
  # t_end = ts.from_datetime(dt_end)

  # Plot pass ground track
  lines = plt.plot(lons, lats, transform=Geodetic, **plot_opts)

  # Plot arrow indicating motion direction at culmination
  lat1, lon1, _ = subpoint_at_dt(sat, dt_culm)
  lat2, lon2, _ = subpoint_at_dt(sat, dt_culm+dtm.timedelta(seconds=1))
  # plt.scatter([lon], [lat], transform=Geodetic)
  PlateCarree_mpl = PlateCarree._as_mpl_transform(ax)
  plt.annotate("", xy=(lon2, lat2), xytext=(lon1, lat1), arrowprops=dict(color=lines[0].get_color(), headlength=8, headwidth=8), xycoords=PlateCarree_mpl, textcoords=PlateCarree_mpl)

  # Plot ticks, if requested
  if ticks_interval is not None:
    delta = dtm.timedelta(minutes=ticks_interval)
    _dt = dt_start.replace(second=0, microsecond=0) + delta
    while _dt <= dt_end:
      lat, lon, alt = subpoint_at_dt(sat, _dt)
      plt.scatter([lon], [lat], transform=Geodetic, s=10, **plot_opts)
      if tzone is not None:
        dt_str = _dt.astimezone(tzone).strftime("%H:%M")
      else:
        dt_str = _dt.strftime("%H:%M")
      plt.annotate(dt_str, xy=(lon, lat), xytext=(6, 0), textcoords="offset pixels", xycoords=PlateCarree_mpl, va="center")
      _dt += delta


# ==============================================================================

# Returns the minimum distance along the surface of Earth required 
# for a satellite at the given altitude in km to be alpha degrees
# above the horizon
def altitude_radius(altitude, alpha):
  R = 6372.0
  alpha_rad = radians(alpha)
  gamma = asin(R/(R+altitude)*cos(alpha_rad))
  theta = pi/2 - (alpha_rad + gamma)
  return theta * R

# ==============================================================================