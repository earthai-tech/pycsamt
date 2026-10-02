"""Prior artifact: one station response figure saved as station.png at 100 dpi."""
import matplotlib

matplotlib.use("Agg")

from pycsamt.emtools import ensure_sites
from pycsamt.emtools.inspect import plot_station_response

sites = ensure_sites("data/3edis", verbose=0)
fig = plot_station_response(sites)
fig.savefig("station.png", dpi=100, bbox_inches="tight")
