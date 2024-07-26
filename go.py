#!/usr/bin/python

import csv
import numpy as np
from pyproj import Proj
from scipy.spatial import distance

CSV_FILENAME = "counties.csv"
OUTPUT_FILENAME = "distances.csv"

# returns (headers, data) tuple
def read_data():
  with open(CSV_FILENAME, newline='') as csvfile:
    r = csv.reader(csvfile)
    headers = next(r, None) # pull the first row as headers
    rawdata = [row for row in r] # read everything else
    return (headers, np.array(rawdata))


def lats_longs(dat):
  lat_longs = np.array(dat[:, (3,2)]).astype('float')
  return (list(lat_longs[:, 0]), list(lat_longs[:, 1]))

def get_xys(lats, longs):
  # https://proj.org/en/9.4/operations/projections/gs50.html
  p = Proj("+proj=gs50")

  [xs, ys] = p(longs, lats)
  return list(zip(xs, ys))


def go():
  print("Reading data")
  _, d = read_data()
  codes = d[:, 0]
  lats, longs = lats_longs(d)

  print("Performing map projection")
  xys = get_xys(lats, longs)

  print("Computing distances")
  dists = distance.cdist(xys, xys)
  
  print("Writing output")
  with open(OUTPUT_FILENAME, "w") as out:
    out.write(','.join(codes) + "\n")
    for rn, row in enumerate(dists):
      print(f"{rn} / {len(dists)}")
      for i, e in enumerate(row.astype('str')):
        out.write(e)
        if i != len(row)-1:
          out.write(",")
      out.write("\n")


if __name__ == "__main__":
  go()
