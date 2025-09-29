#!/usr/bin/env python3
from __future__ import annotations

"""
File: find_nearest_cbsa_disttance.py
Author(s): Amin Bemanian
Date: 07/07/24
Description: GPT assisted script which finds the closest CBSA distance
between states as an alternative to Euclidean centroid distance
"""

"""
cbsa_nearest_distance_by_state_fullnames_sym.py
==============================================

Outputs:
1. **cbsa_points.csv** – every CBSA with centroid lat/lon and the set of
   full state names it spans.

2. **state_cbsa_dist.tsv** – symmetric table: two rows
   for each unordered state pair (state1, state2) and (state2, state1),
   with full state names and the minimum CBSA‑to‑CBSA distance (km).

Data source: 2024 CBSA Gazetteer
https://www2.census.gov/geo/docs/maps-data/data/gazetteer/2024_Gazetteer/2024_Gaz_cbsa_national.zip
"""
import io, math, zipfile, re
from itertools import combinations

import certifi, pandas as pd, requests, urllib3
from requests.exceptions import SSLError

urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

GAZ_URL = (
    "https://www2.census.gov/geo/docs/maps-data/data/gazetteer/"
    "2024_Gazetteer/2024_Gaz_cbsa_national.zip"
)

# State abbreviation → full name
ABBR_TO_FULL = {
    "AL":"Alabama","AK":"Alaska","AZ":"Arizona","AR":"Arkansas","CA":"California",
    "CO":"Colorado","CT":"Connecticut","DE":"Delaware","DC":"District of Columbia",
    "FL":"Florida","GA":"Georgia","HI":"Hawaii","ID":"Idaho","IL":"Illinois",
    "IN":"Indiana","IA":"Iowa","KS":"Kansas","KY":"Kentucky","LA":"Louisiana",
    "ME":"Maine","MD":"Maryland","MA":"Massachusetts","MI":"Michigan",
    "MN":"Minnesota","MS":"Mississippi","MO":"Missouri","MT":"Montana",
    "NE":"Nebraska","NV":"Nevada","NH":"New Hampshire","NJ":"New Jersey",
    "NM":"New Mexico","NY":"New York","NC":"North Carolina","ND":"North Dakota",
    "OH":"Ohio","OK":"Oklahoma","OR":"Oregon","PA":"Pennsylvania","RI":"Rhode Island",
    "SC":"South Carolina","SD":"South Dakota","TN":"Tennessee","TX":"Texas",
    "UT":"Utah","VT":"Vermont","VA":"Virginia","WA":"Washington",
    "WV":"West Virginia","WI":"Wisconsin","WY":"Wyoming",
}
STATE_ABBRS = sorted(ABBR_TO_FULL.keys())

def safe_get(url:str, **kw):
    try:
        return requests.get(url, timeout=60, verify=certifi.where(), **kw)
    except SSLError:
        print("TLS verification failed – retrying without verification …")
        return requests.get(url, timeout=60, verify=False, **kw)

def extract_states(name:str) -> set[str]:
    # take tokens after last comma, collect all two‑letter codes
    tail = name.split(",")[-1]
    return {ABBR_TO_FULL[t] for t in re.findall(r"\b[A-Z]{2}\b", tail) if t in ABBR_TO_FULL}

def haversine(lat1,lon1,lat2,lon2):
    R=6371.0088
    φ1,φ2=math.radians(lat1),math.radians(lat2)
    dφ=φ2-φ1
    dλ=math.radians(lon2-lon1)
    a=math.sin(dφ/2)**2 + math.cos(φ1)*math.cos(φ2)*math.sin(dλ/2)**2
    return 2*R*math.asin(math.sqrt(a))

# --- download Gazetteer
print("Downloading CBSA Gazetteer …")
resp=safe_get(GAZ_URL); resp.raise_for_status()
with zipfile.ZipFile(io.BytesIO(resp.content)) as zf:
    txt=next(n for n in zf.namelist() if n.endswith(".txt"))
    with zf.open(txt) as fh:
        gaz=pd.read_csv(fh,sep="\t",dtype=str,engine="python",na_filter=False)

gaz.columns=gaz.columns.str.strip()
lat_col ="INTPTLAT"  if "INTPTLAT"  in gaz.columns else "LAT"
lon_col ="INTPTLONG" if "INTPTLONG" in gaz.columns else "LON"

gaz=gaz[["GEOID","NAME",lat_col,lon_col]].rename(
    columns={lat_col:"lat",lon_col:"lon"})
gaz["lat"]=gaz["lat"].astype(float)
gaz["lon"]=gaz["lon"].astype(float)

# derive full‑name state sets
gaz["states"]=gaz["NAME"].apply(extract_states)
gaz=gaz[gaz["states"].apply(bool)].reset_index(drop=True)
gaz.to_csv("data/cbsa_points.csv",index=False)
print("Wrote cbsa_points.csv  (rows:",len(gaz),")")

# coordinates dict
coords=gaz.set_index("GEOID")[["lat","lon"]].to_dict("index")
cbsa_ids=list(coords)

# cache distances with lexicographically sorted key
dist={}
import math
for i,c1 in enumerate(cbsa_ids):
    lat1,lon1=coords[c1]["lat"],coords[c1]["lon"]
    for c2 in cbsa_ids[i+1:]:
        lat2,lon2=coords[c2]["lat"],coords[c2]["lon"]
        key=(c1,c2) if c1<c2 else (c2,c1)
        dist[key]=haversine(lat1,lon1,lat2,lon2)

# nearest distances – produce symmetric rows with full names
records=[]
state_full_list=sorted(ABBR_TO_FULL.values())
# map full→abbrev for faster filtering
FULL_TO_ABBR={v:k for k,v in ABBR_TO_FULL.items()}

for s1_full in state_full_list:
    cbsas1 = gaz[gaz["states"].apply(lambda st: s1_full in st)]["GEOID"]

    for s2_full in state_full_list:
        # self–self pair → distance 0, add row and skip the rest
        if s1_full == s2_full:
            records.append(
                {"x": s1_full, "y": s2_full, "min_cbsa_dist": 0}
            )
            continue                       # go to next s2_full

        #----- only for different states -----
        best = math.inf
        cbsas2 = gaz[gaz["states"].apply(lambda st: s2_full in st)]["GEOID"]

        for a in cbsas1:
            for b in cbsas2:
                if a == b:
                    continue
                key = (a, b) if a < b else (b, a)
                d = dist.get(key)
                if d is not None and d < best:
                    best = d

        if math.isfinite(best):
            records.append(
                {"x": s1_full, "y": s2_full, "min_cbsa_dist": round(best, 2)}
            )
pd.DataFrame(records).to_csv("data/state_cbsa_dist.tsv",sep="\t",index=False)
print("Wrote state_cbsa_dist.tsv  (rows:",len(records),")")
