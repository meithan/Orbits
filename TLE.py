# Ddecods orbital data given in Two-Line Element (TLE) format

from datetime import datetime, timedelta
from math import pi, sqrt
import string
import sys

# ====================================================================

# Constants, all in MKS
SECS_PER_DAY = 86400
EARTH_MU = 3.986004418e14
EARTH_R = 6378.137e3

# ====================================================================

# Computes a TLE line checksum
def line_checksum(line):
  s = 0
  for c in line[:68]:
    if c in string.digits:
      s += int(c)
    elif c == "-":
      s += 1
  return (s % 10)

# Checks the validity of the TLE, reporting any problems found
def validate(line1, line2):

  for n,line in ((1,line1), (2,line2)):

    if line[0] != str(n):
      return False, "Line %i should start with %i (found: %s)" % (n, n, line[0])

    if len(line) != 69:
      return False, "Line %i has incorrect length (%i, must be 69)" % (n, len(line))

    computed_sum = str(line_checksum(line))
    expected_sum = line[68]
    if computed_sum != expected_sum:
      return False, "Checksum mismatch in line 1 (computed: %s, expected: %s)" % (computed_sum, expected_sum)

  if line1[2:7] != line2[2:7]:
    return False, "Line 1's sat number (%s) does not match line 2's (%s)" % (line1[2:7], line2[2:7])

  return True, ""

# Converts text year and fractional day to a datetime object
def compute_epoch(epoch_year, epoch_day):
  year = int(epoch_year)
  if year < 57:
    year += 2000
  else:
    year += 1900
  days = float(epoch_day) - 1
  return datetime(year, 1, 1) + timedelta(days=days)

# Parses implied-decimal exponential notation
def parse_exp(s):
  if s[0] == "-":
    sgn = -1
  else:
    sgn = +1
  f = float("0." + s[1:-2].strip() + "e" + s[-2:])
  return sgn*f

# Prints the error message
def show_error_msg(TLE):
  print("Bad input: TLE must be either a string with 2 or 3 lines separated by a newline character, or a list with either 2 or 3 elements; the 3-lines forms assumes the first line contains the title")
  print("Given input:")
  print(TLE)

# ===========================

# Takes a TLE and outputs a dict with the decoded satellite info and orbital
# elements.
# The input TLE must be either:
#  1) A string containing either 2 or 3 lines separated by newline characters
#  2) A list or tuple of either 2 or 3 strings, each a single line of the TLE
# The 3-lines versions assume that the first line contains the title (usually
# the common name of the satellite).
# The output is a dict with decoded orbital elemens and rate, indexed by the
# following keys and meanings:
#  "catnum": NORAD/USSPACECOM satellite catalog number (string)
#  "class": Classification (U=Unclassifief, C=Classified, S=Secret) (string)
#  "intldes": International designator (COSPAR ID) (string)
#  "epoch": epoch of orbital elements (tz-naive UTC datetime object)
#  "ndot": first derivative of mean motion in revs/day (float)
#  "nddot": second derivative of mean motion (float)
#  "bstar": drag term (stared balistic coef.) or radiation pressure coef. (float)
# "n", "i", "e", "arg", "LAN", "M0", "P", "a", "rpe", "rap", "hpe", "hap", "revs"
def decode(TLE, check=True, verbose=False, debug=False):

  # Check input, separate parts
  if type(TLE) == str:
    if TLE.count("\n") not in [1, 2]:
      show_error_msg(TLE)
      return None
    else:
      tokens = TLE.split("\n")
  elif type(TLE) == list:
    if len(TLE) not in [2, 3]:
      show_error_msg(TLE)
      return None
    else:
      tokens = TLE
  else:
    show_error_msg(TLE)
    return None
  if len(tokens) == 2:
    TLE_title = ""
    TLE_line1 = tokens[0].strip("\n")
    TLE_line2 = tokens[1].strip("\n")
  else:
    TLE_title = tokens[0].strip("\n")
    TLE_line1 = tokens[1].strip("\n")
    TLE_line2 = tokens[2].strip("\n")

  # Validate TLE
  if check:
    valid, msg = validate(TLE_line1, TLE_line2)
    if valid:
      if verbose: print("TLE valid")
    else:
      print("Invalid TLE!")
      print(msg)
      print(TLE)

  # Parse line 1
  satellite_number = TLE_line1[2:7]
  classification = TLE_line1[7:8]
  intl_designator_raw = TLE_line1[9:17]
  epoch_year = TLE_line1[18:20]
  epoch_day = TLE_line1[20:32]
  ndot_over_2 = TLE_line1[33:43]
  nddot_over_6 = TLE_line1[44:52]
  bstar_drag = TLE_line1[53:61]
  elset_number = TLE_line1[64:68]

  # Parse line 2
  TLE_inc = TLE_line2[8:16]
  TLE_RAAN = TLE_line2[17:25]
  TLE_ecc = TLE_line2[26:33]
  TLE_arg = TLE_line2[34:42]
  TLE_M0 = TLE_line2[43:51]
  TLE_n = TLE_line2[52:63]
  TLE_revs = TLE_line2[63:68]

  if TLE_title == "":
    TLE_title = satellite_number

  if debug:
    print("Parsed TLE:")
    print("NORAD catalog no: '%s'" % satellite_number)
    print("Classification:   '%s'" % classification)
    print("Intl designator:  '%s'" % intl_designator_raw)
    print("Epoch year:       '%s'" % epoch_year)
    print("Epoch day:        '%s'" % epoch_day)
    print("ndot/2:           '%s'" % ndot_over_2.strip())
    print("nddot/6:          '%s'" % nddot_over_6.strip())
    print("bstar drag:       '%s'" % bstar_drag.strip())
    print("elset no:         '%s'" % elset_number.strip())
    print("inc:              '%s'" % TLE_inc)
    print("raan:             '%s'" % TLE_RAAN)
    print("ecc:              '%s'" % TLE_ecc)
    print("arg:              '%s'" % TLE_arg)
    print("M0:               '%s'" % TLE_M0)
    print("n:                '%s'" % TLE_n)
    print("revs:             '%s'" % TLE_revs)

  yr = int(intl_designator_raw[:2])
  yr += 2000 if yr < 57 else 1900
  intl_designator = "{}-{}".format(yr, intl_designator_raw[2:])

  # Compute orbital elements
  epoch = compute_epoch(epoch_year, epoch_day)
  ndot = 2 * float(ndot_over_2)
  nddot = 6 * parse_exp(nddot_over_6)
  bstar = parse_exp(bstar_drag)
  n = float(TLE_n)
  inc = float(TLE_inc)
  ecc = float("0." + TLE_ecc)
  arg = float(TLE_arg)
  raan = float(TLE_RAAN)
  M0 = float(TLE_M0)
  revs = int(TLE_revs)
  P = 1/n * SECS_PER_DAY
  a = (EARTH_MU * (P/(2*pi))**2)**(1./3.)
  rpe = a*(1-ecc)
  rap = a*(1+ecc)
  hpe = rpe - EARTH_R
  hap = rap - EARTH_R
  vpe = sqrt((1+ecc)*EARTH_MU/((1-ecc)*a))
  vap = sqrt((1-ecc)*EARTH_MU/((1+ecc)*a))

  result = {"name": TLE_title, "epoch": epoch, "ndot": ndot, "nddot": nddot, "bstar": bstar, "n": n, "i": inc, "e": ecc, "arg": arg, "LAN": raan, "M0": M0, "P": P, "a": a, "rpe": rpe, "rap": rap, "hpe": hpe, "hap": hap}

  if debug or verbose:
    print()
    if TLE_title != "":
      print("Spacecraft name:  %s" % TLE_title)
    print("NORAD catalog no: %s" % satellite_number)
    print("Intl designator:  %s" % intl_designator)
    print("\nOrbital Elements")
    print("Epoch:         %s" % epoch)
    print("Apogee:        {:,.3f} km".format(hap/1e3))
    print("Perigee:       {:,.3f} km".format(hpe/1e3))
    print("Period:        {:,.3f} minutes".format(P/60))
    print("Inclination:   {:.4f}°".format(inc))
    print("Eccentricty:   {:.7f}".format(ecc))
    print("Arg. perigee:  {:.4f}°".format(arg))
    print("RA asc. node:  {:.4f}°".format(raan))
    print("MA at epoch:   {:.4f}°".format(M0))
    print("Mean motion:   {:.4f} rev/day".format(n))
    print("Perigee speed: {:.3f} km/s".format(vpe/1e3))
    print("Apogee speed:  {:.3f} km/s".format(vap/1e3))
    print("ndot:          {:.4e} rev/day^2".format(ndot))
    print("nddot:         {:.4e} rev/day^3".format(nddot))
    print("bstar:         {:.4e} 1/RE".format(bstar))

  return result

# ===============================

if __name__ == "__main__":

  # TLE = \
#"""1 43175U 18012B   18094.46550631 -.00000083  00000-0  00000+0 0  9995
#2 43175  12.1384 210.8310 5567305 258.1238 218.1208  1.16549561  1007"""

  # Jason-3
  #TLE = """1 41240U 16002A   16017.84863520 -.00000062  00000-0  00000+0 0  9999
  #2 41240  66.0451 114.8369 0016419 234.2524 227.7053 12.88437011    01"""

  # Beresheet, 6 march 2019
#   TLE = \
# """Beresheet
# 1 44049U 19009B   19064.00000000 -.00001457  00000-0  00000+0 0  9997
# 2 44049  27.6566   9.9443 9013856 181.3433 301.0220  0.46403215    91"""

# Sentinel 2A
#  TLE = \
#"""SENTINEL-2A
#1 40697U 15028A   19285.17258184 -.00000021  00000-0  86979-5 0  9991
#2 40697  98.5703 357.9797 0001109  84.2547 275.8768 14.30816136224811"""

# ISS on 31 may 2020
  TLE = """ISS (ZARYA)
1 25544U 98067A   20151.61686127  .00000168  00000-0  11087-4 0  9992
2 25544  51.6444  75.4313 0002297  11.5525  50.1151 15.49398617229298"""

  TLE = """COSMOS 1408
1 13552U 82092A   21319.03826954  .00002024  00000-0  69413-4 0  9995
2 13552  82.5637 123.6906 0018570 108.1104 252.2161 15.29390138142807"""

#   TLE = """ISS (ZARYA)
# 1 25544U 98067A   21319.86359406  .00001699  00000-0  39549-4 0  9996
# 2 25544  51.6442 312.7312 0004574 201.3017 242.9752 15.48581038312108"""

  decode(TLE, debug=True)
