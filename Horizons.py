# Functions to parse JPL Horizons data file

from datetime import datetime, timedelta
import sys
import numpy as np

import Orbit
import Utils

# ==============================================================================

# Loads a JPL Hozirons CSV data file, parsing the info into a dict
# If dates are TDB, this will convert them to UTC (if delta-T not in CSV file,
# will use a reasonable default)
def read_Horizons_csv(fname, convert_numpy=True, convert_TDB_UTD=False):

  import csv
  if convert_numpy:
    import numpy as np

  def parse(values):
    count = 0
    result = []
    for x in values:
      if x.strip() == "":
        count += 1
        result.append(f"unknown{count}")
      else:
        result.append(x.strip())
    return result

  # Find header and read available fields
  f = open(fname)
  start_line = None
  rev_date = None
  for no,line in enumerate(f):
    if no < 10 and "Revised" in line:
      rev_date = Utils.parse_date(line[10:23], date_only=True)
    if line.startswith("$$SOE"):
      start_line = no
      break
  if start_line is None:
    print("Couldn't find $$SOE in {}".format(fname))
    sys.exit()
  f.seek(0)
  for i in range(start_line-2):
    f.readline()
  reader = csv.reader(f)
  header = parse(next(reader))
  # print(header)

  header_idxs = [k for k in range(len(header)) if header[k].strip() != ""]
  dataset = {}
  for k in header_idxs:
    dataset[header[k]] = []

  gave_warning_deltaT = False
  for row in reader:

    row = [x.strip() for x in row]
    if row[0].startswith("*") or row[0] == "$$SOE":
      continue
    elif row[0] == "$$EOE":
      break

    for k in header_idxs:

      if header[k] in ["Calendar Date (TDB)", "Date__(UT)__HR:MN"]:
        
        if "Date" not in dataset:
          dataset["Date"] = []
        dataset[header[k]].append(row[k])

        date = Utils.parse_date(row[k])

        if convert_TDB_UTD:

          dataset[header[k]].append(row[k])
          
          # Convert TDB to UTC         
          if "delta-T" in header:
            TDB_minus_UTC = float(row[header.index("delta-T")])
          else:
            # Full formula is:
            # TDB = UTC + leap_secs + 32.184 + 0.001658 sin (g) + 0.000014 sin (2g) seconds
            # where g = 357.53 + 0.9856003 ( JD – 2451545.0 ) degrees
            # Here we ignore the small time-varying (periodic) component
            leap_secs = 37   # As of 2022
            def_deltaT = leap_secs + 32.184
            if not gave_warning_deltaT:
              print(f"Warning: delta-T not present in file; using default value of {def_deltaT} s")
              gave_warning_deltaT = True
            TDB_minus_UTC = def_deltaT
          
          date -= timedelta(seconds=TDB_minus_UTC)        
          dataset["Date"].append(date)
        
        else:

          dataset["Date"].append(date)

      elif header[k] == "Date__(UT)__HR:MN":

        dataset[header[k]].append(row[k])
        date = Utils.parse_date(row[0])
        
        if "Date" not in dataset:
          dataset["Date"] = []
        dataset["Date"].append(date)        

      else:
        
        try:
          val = float(row[k])
        except:
          val = row[k]
        
        dataset[header[k]].append(val)

  # dataset["dates"] = dataset["Calendar Date (TDB)"]

  if "Date" not in dataset:
    print("Error: couldn't find a Date field")
    sys.exit()

  if convert_numpy:
    for key in dataset:
      if isinstance(dataset[key][0], (float, int)):
        dataset[key] = np.array(dataset[key])

  if rev_date is not None:
    dataset["rev_date"] = rev_date

  return dataset

# ==============================================================================

# Loads JPL Horizons data from file and interpolates values of any variable at
# arbitrary dates (in the range of the data)
class HorizonsData:

  # ----------------------------------------------------------------------------

  def __init__(self, horizons_fname=None, length_scale=None):

    self.data = None
    self.idx = None
    self.vars = []

    if horizons_fname is not None: 
      self.load_data(horizons_fname, length_scale)

  def __getitem__(self, key):
    return getattr(self, key)

  def __repr__(self):
    s = "<HorizonData object, fields: "
    s += ", ".join(f"'{v}'" for v in self.vars)
    s += ">"
    return s

  # ----------------------------------------------------------------------------

  def load_data(self, horizons_fname, length_scale=None):
    dataset = read_Horizons_csv(horizons_fname)
    for key in dataset:
      setattr(self, key, dataset[key])
      self.vars.append(key)
    self.dates = self.Date
    self.N = len(self.dates)
    self.min_date = self.dates[0]
    self.max_date = self.dates[-1]
    if length_scale is not None:
      if hasattr(self, 'X'): self.X /= length_scale
      if hasattr(self, 'X'): self.Y /= length_scale
      if hasattr(self, 'Z'): self.Z /= length_scale

  # ----------------------------------------------------------------------------

  # Interpolates a variable at the given date (which must be in range)
  # varnames can be either a string with a single variable name, or a list
  # of strings to retrieve multiple variables at the same date
  def value_at_date(self, varnames, date):

    if date < self.min_date:
      raise ValueError("Given date ({}) before min data date ({})".format(date, self.min_date))
      return None

    if date > self.max_date:
      raise ValueError("Given date ({}) after max data date ({})".format(date, self.max_date))
      return None

    if type(varnames) == str:
      varnames = [varnames]

    for var in varnames:
      if not hasattr(self, var):
        print(self)
        err_msg = "Dataset has no variable '{}'".format(var)
        raise Exception(err_msg)
        # return None

    # Find index of bracketing dates
    self.idx = self.find_idx(date)

    # Interpolate values
    i1 = self.idx
    i2 = self.idx + 1
    dt1 = self.dates[i1]
    dt2 = self.dates[i2]
    dt2_dt1 = (dt2-dt1).total_seconds()
    d_dt1 = (date-dt1).total_seconds()
    values = []
    for var in varnames:
      y1 = getattr(self, var)[i1]
      y2 = getattr(self, var)[i2]
      m = (y2-y1)/dt2_dt1
      y = y1 + m*d_dt1
      values.append(y)

    if len(values) == 1:
      return values[0]
    else:
      return values

  # ----------------------------------------------------------------------------

  def pos_at_date(self, date, use_numpy=True):
    x, y, z = self.value_at_date(['X', 'Y', 'Z'], date)
    return np.array((x, y, z)) if use_numpy else (z, y, z)

  def vel_at_date(self, date, use_numpy=True):
    vx, vy, vz =  self.value_at_date(['VX', 'VY', 'VZ'], date)
    return np.array((vx, vy, vz)) if use_numpy else (vz, vy, vz)

  def posvel_at_date(self, date, use_numpy=True):
    x, y, z, vx, vy, vz = self.value_at_date(['X', 'Y', 'Z', 'VX', 'VY', 'VZ'], date)
    return (np.array((x, y, z)), np.array((vx, vy, vz))) if use_numpy else (x, y, z, vx, vy, vz)

  def orbit_at_date(self, primary, date):
    A, EC, IN, W, OM, MA = self.value_at_date(['A', 'EC', 'IN', 'W', 'OM', 'MA'], date)
    elements = {'a': A*1e3, 'e': EC, 'i': IN, 'arg': W, 'LAN': OM, 'M0': MA, 'epoch': date}
    return Orbit(primary, elements)

  # ----------------------------------------------------------------------------

  # Finds in the index that brackets (idx and idx+1) a given date
  def find_idx(self, date):

    # Try last index or the one after that
    if self.idx is not None:
      if self.dates[self.idx] <= date <= self.dates[self.idx+1]:
        return self.idx
      elif (self.idx != self.N-2) and (self.dates[self.idx+1] <= date <= self.dates[self.idx+2]):
        return self.idx + 1

    # Replace with binary search, man!
    self.idx = None
    for i in range(self.N-1):
      if self.dates[i] <= date <= self.dates[i+1]:
        self.idx = i
        break

    return self.idx

# ==============================================================================

# ---DEPRECATED---
# The above code is more general
#
# Loads data from a JPL Horizons output file
# columns are the names of the columns present in the file fname, which
# must be in CSV format. The header (object page) will be skipped.
# Returned is a list of dicts with the elements labeled as in columns.
def load_JPL_data(fname, columns):
  column_names = {'JDTDB': 'jd', 'Calendar Date': 'date', 'Calendar Date (TDB)': 'date', 'delta-T': 'TDBUT', 'X': 'x', 'Y': 'y', 'Z': 'z', 'VX': 'vx', 'VY': 'vy', 'VZ': 'vz', 'EC': 'e', "QR": 'rpe', "IN": 'i', "OM": 'LAN', 'W': 'arg', 'Tp': "tp", 'N': 'n', 'MA': 'M0', 'TA': 'nu', 'A': 'a', 'AD': 'rap', 'PR': 'P', "LT": "LT", "RG": "RG", "RR": "RR"}
  data = []
  reading_data = False
  f = open(fname)
  while True:
    line = f.readline()
    if line == "": break
    line = line.strip()
    if line == "$$SOE":
      reading_data = True
      continue
    elif line == "$$EOE":
      break
    if reading_data:
      values = {}
      tokens = line.split(",")
      for i,colname in enumerate(columns):
        if "Calendar Date" in colname:
          value = datetime.strptime(tokens[i][6:], "%Y-%b-%d %H:%M:%S.%f")
        else:
          value = float(tokens[i])
          if colname in ['QR', 'A', 'AD']:
            value *= 1000
        if colname in column_names:
          colname = column_names[colname]
        values[colname] = value
      data.append(values)
  f.close()
  return data
