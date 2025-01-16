# A simple module to retrieve TLEs from Space-Track.org

import datetime as dtm
import html
import requests
import sys

import TLE

# ==============================================================================

# Downloads TLEs from Space-Track.org by NORAD catalog ID
# If only the ID is given, will download the latest TLE
# If either min_date or max_date are given, will download all TLEs in range
# If download_all is set to True, will download all available TLEs
def download_TLE(norad_cat_id, credentials={}, out_fname="tle.txt", download_all=False, min_date=None, max_date=None):

  # Load credentials
  username = None; password = None
  if "username" in credentials and "password" in credentials:
    username = credentials["username"]
    password = credentials["password"]
  else:
    try:
      from space_track_creds import username, password
    except:
      pass
  if username in [None, ""]  or password in [None, ""]:
    print("Erorr: Space-Track login credentials couldn't be loadad. Either fill them in space_track_creds.py or provide them through the 'credentials' argument.")
    sys.exit()

  # Generate query
  query_url = "https://www.space-track.org/basicspacedata/query/class/"

  if download_all or min_date is not None or max_date is not None:
    query_url += "gp_history/"
  else:
    query_url += "gp/"

  if not download_all and (min_date is not None or max_date is not None):
    parts = []
    if min_date is not None:
      parts.append(f">{min_date.strftime('%Y-%m-%d')}")
    if max_date is not None:
      parts.append(f"<{max_date.strftime('%Y-%m-%d')}")
    query_url += f"EPOCH/{','.join(parts)}/"

  query_url += f"NORAD_CAT_ID/{norad_cat_id}/"
  query_url += "orderby/EPOCH DESC/format/tle"

  # Authenticate
  print("\nAuthenticating with Space-Track.org ...", end=" ")
  auth_url = "https://www.space-track.org/ajaxauth/login"
  login_data = {'identity': username, 'password': password}
  session = requests.Session()
  response = session.post(auth_url, data=login_data)
  if response.status_code != 200:
    print("Error!")
    print(response)
    return
  else:
    print("Done")
    # print("Headers:", response.headers)
    # print("Content:", response.content)
    # print("Cookies:", response.cookies)
    # print("Session cookies:", session.cookies)

  # Run query
  print("\nRunning query ...")
  print(query_url)
  response = session.get(query_url)
  if response.status_code != 200:
    print("Request returned error!")
    print("Status code:", response.status_code)
    return

  # Write TLE to file
  with open(out_fname, "w") as fout:
    fout.write(response.text)
  print("\nWrote", out_fname)

# ==============================================================================

if __name__ == "__main__":

  # Example: download latest TLE for the iSS

  # ISS NORAD ID
  norad_cat_id = 25544

  download_TLE(norad_cat_id)

