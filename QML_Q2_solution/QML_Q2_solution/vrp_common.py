"""
ME44207 Quantitative Methods for Logistics - Group Assignment Q2 (2026)
Group number : <GROUP NUMBER>
Members      : <NAME 1>, <NAME 2>, <NAME 3>, <NAME 4>

vrp_common.py
-------------
Shared code used by all parts:
  * reading the data files AS IS (geo_data.txt, geo_data_large.txt, range_data.txt)
  * the Instance class (node sets, arcs, distances, time windows)
  * the three-index VRP(TW) building blocks in Gurobi (x, D, T variables + constraints)
  * route extraction, an independent schedule/range re-computation, tables and plots

Notation (same as in the report):
  node 0      = start depot (first row of the data file)
  nodes 1..n  = pickup locations (customers)
  node n+1    = copy of the depot used as end depot; it is created in Python only,
                the input file itself is not changed.
"""
import os
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import gurobipy as gp
from gurobipy import GRB

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(HERE, "output")
os.makedirs(OUT_DIR, exist_ok=True)

pd.set_option("display.width", 200)
pd.set_option("display.max_columns", 30)
pd.set_option("display.max_rows", 200)

STATUS_TXT = {GRB.OPTIMAL: "OPTIMAL", GRB.INFEASIBLE: "INFEASIBLE",
              GRB.TIME_LIMIT: "TIME_LIMIT", GRB.INF_OR_UNBD: "INF_OR_UNBD"}


# --------------------------------------------------------------------------------------
# Data
# --------------------------------------------------------------------------------------
def data_path(fname):
    """Data files are expected in the same folder as the .py files."""
    return os.path.join(HERE, fname)


def load_geo(fname="geo_data.txt"):
    """Read a location file as is: whitespace separated, rows may have trailing tabs and
    the file may end with an empty line. Columns: LOC_ID X Y READY DUE SERVICE."""
    rows = []
    with open(data_path(fname)) as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 6:
                rows.append([float(p) for p in parts[:6]])
    df = pd.DataFrame(rows, columns=["LOC_ID", "XCOORD", "YCOORD",
                                     "READYTIME", "DUETIME", "SERVICETIME"])
    df["LOC_ID"] = df["LOC_ID"].astype(int)
    return df


def load_range(fname="range_data.txt"):
    rows = []
    with open(data_path(fname)) as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 2:
                rows.append([parts[0], float(parts[1])])
    return pd.DataFrame(rows, columns=["MONTH", "RANGE"])


class Instance:
    """Holds sets and parameters of one problem instance.

    time_windows=False -> all time-window constraints are dropped (and no arcs are pruned).
    """

    def __init__(self, geo, time_windows=True):
        self.geo = geo
        self.tw = time_windows
        self.n = len(geo) - 1                        # number of pickup locations
        self.C = list(range(1, self.n + 1))          # customers
        self.end = self.n + 1                        # end-depot copy (node n+1)
        self.V = [0] + self.C + [self.end]           # all nodes incl. both depot copies

        xy = geo[["XCOORD", "YCOORD"]].to_numpy()
        xy = np.vstack([xy, xy[0]])                  # node n+1 has the depot coordinates
        self.xy = xy
        self.e = list(geo["READYTIME"]) + [geo["READYTIME"].iloc[0]]
        self.l = list(geo["DUETIME"]) + [geo["DUETIME"].iloc[0]]
        self.s = list(geo["SERVICETIME"]) + [0.0]

        # Euclidean distance = travel time (1 distance unit = 1 time unit)
        self.d = {(i, j): math.hypot(xy[i, 0] - xy[j, 0], xy[i, 1] - xy[j, 1])
                  for i in self.V for j in self.V}

        # Arc set A:  (0,j), (i,j) between customers, (i,n+1), plus (0,n+1) = "vehicle unused"
        arcs = [(0, j) for j in self.C] + [(0, self.end)]
        arcs += [(i, j) for i in self.C for j in self.C if i != j]
        arcs += [(i, self.end) for i in self.C]
        if self.tw:
            # pre-processing: an arc that can never be time-feasible is removed
            arcs = [(i, j) for (i, j) in arcs
                    if self.e[i] + self.s[i] + self.d[i, j] <= self.l[j] + 1e-9]
        self.arcs = arcs
        self.out_arcs = {i: [j for (a, j) in arcs if a == i] for i in self.V}
        self.in_arcs = {j: [i for (i, b) in arcs if b == j] for j in self.V}

    def route_length_upper_bound(self):
        """Valid upper bound on the length of any route (used as big-M when the range
        is not a hard constraint, i.e. in part 3)."""
        return sum(max(self.d[i, j] for j in self.out_arcs[i]) for i in [0] + self.C)


# --------------------------------------------------------------------------------------
# Model building blocks (three-index formulation)
# --------------------------------------------------------------------------------------
def build_routing_model(inst, K, name, Dmax, symmetry_breaking=True):
    """Creates the routing part shared by all models:
         x[i,j,k]  binary  vehicle k uses arc (i,j)
         D[i,k]    >= 0    cumulative distance of vehicle k on arrival at node i (<= Dmax)
         T[i,k]            start of service of vehicle k at node i (only with time windows)
    Constraints: visit once, depot out/in, flow conservation, distance propagation
    (MTZ-type, also eliminates subtours), time propagation + time windows.
    The objective and the range / recourse constraints are added by the separate parts."""
    m = gp.Model(name)
    Ks = range(K)
    A, C, end, d = inst.arcs, inst.C, inst.end, inst.d

    x = m.addVars([(i, j, k) for (i, j) in A for k in Ks], vtype=GRB.BINARY, name="x")
    # tightened bounds on D (valid inequalities): a vehicle that arrives at customer i has
    # driven at least d_0i and must still be able to drive back, so d_0i <= D_ik <= Dmax - d_i0
    Dlb = {i: (d[0, i] if i in C else 0.0) for i in inst.V}
    Dub = {i: (Dmax - d[i, end] if i in C else Dmax) for i in inst.V}
    D = m.addVars([(i, k) for i in inst.V for k in Ks],
                  lb={(i, k): Dlb[i] for i in inst.V for k in Ks},
                  ub={(i, k): Dub[i] for i in inst.V for k in Ks}, name="D")

    # (C1) every pickup location is visited exactly once
    for i in C:
        m.addConstr(gp.quicksum(x[i, j, k] for j in inst.out_arcs[i] for k in Ks) == 1,
                    name=f"visit[{i}]")
    for k in Ks:
        # (C2) every vehicle leaves the start depot once (arc (0,n+1) = vehicle not used)
        m.addConstr(gp.quicksum(x[0, j, k] for j in inst.out_arcs[0]) == 1, name=f"start[{k}]")
        # (C4) every vehicle arrives once at the end depot
        m.addConstr(gp.quicksum(x[i, end, k] for i in inst.in_arcs[end]) == 1, name=f"end[{k}]")
        # (C3) flow conservation at the pickup locations
        for h in C:
            m.addConstr(gp.quicksum(x[i, h, k] for i in inst.in_arcs[h]) ==
                        gp.quicksum(x[h, j, k] for j in inst.out_arcs[h]), name=f"flow[{h},{k}]")
        # (C5) vehicles start with an empty odometer
        m.addConstr(D[0, k] == 0, name=f"D0[{k}]")
        # (C6) cumulative distance propagation (also eliminates subtours since d_ij > 0)
        #      big-M = UB(D_i) + d_ij - LB(D_j)
        for (i, j) in A:
            M = max(0.0, Dub[i] + d[i, j] - Dlb[j])
            m.addConstr(D[j, k] >= D[i, k] + d[i, j] - M * (1 - x[i, j, k]),
                        name=f"dist[{i},{j},{k}]")

    # optional symmetry breaking (vehicles are identical): number the vehicles by
    # non-increasing route length, L_k >= L_k+1. This also implies that vehicle k+1 is
    # only used if vehicle k is used.
    if symmetry_breaking:
        for k in range(K - 1):
            m.addConstr(gp.quicksum(d[i, j] * x[i, j, k] for (i, j) in A) >=
                        gp.quicksum(d[i, j] * x[i, j, k + 1] for (i, j) in A), name=f"sym[{k}]")

    T = None
    if inst.tw:
        T = m.addVars([(i, k) for i in inst.V for k in Ks],
                      lb={(i, k): inst.e[i] for i in inst.V for k in Ks},
                      ub={(i, k): inst.l[i] for i in inst.V for k in Ks}, name="T")   # (C8)
        for k in Ks:
            for (i, j) in A:
                # (C7) start time propagation, big-M = max(0, l_i + s_i + t_ij - e_j)
                M = max(0.0, inst.l[i] + inst.s[i] + d[i, j] - inst.e[j])
                m.addConstr(T[j, k] >= T[i, k] + inst.s[i] + d[i, j] - M * (1 - x[i, j, k]),
                            name=f"time[{i},{j},{k}]")
    return m, x, D, T


def route_length_expr(inst, x, k):
    return gp.quicksum(inst.d[i, j] * x[i, j, k] for (i, j) in inst.arcs)


def total_distance_expr(inst, x, K):
    return gp.quicksum(inst.d[i, j] * x[i, j, k] for (i, j) in inst.arcs for k in range(K))


def set_params(m, time_limit=1800, verbose=False):
    m.Params.OutputFlag = 1 if verbose else 0
    m.Params.TimeLimit = time_limit          # 30 min limit as required by the assignment


# --------------------------------------------------------------------------------------
# Solution extraction and (independent) re-computation of times and ranges
# --------------------------------------------------------------------------------------
def extract_routes(inst, x, K):
    """Follows the x-arcs of every vehicle from node 0 to node n+1. Empty routes = []."""
    routes = []
    for k in range(K):
        succ = {i: j for (i, j) in inst.arcs if x[i, j, k].X > 0.5}
        route, cur = [], succ.get(0, inst.end)
        while cur != inst.end:
            route.append(cur)
            cur = succ[cur]
        routes.append(route)
    return routes


def route_distance(inst, route):
    path = [0] + route + [inst.end]
    return sum(inst.d[a, b] for a, b in zip(path[:-1], path[1:]))


def schedule(inst, route, use_tw=True):
    """Recomputes the timetable of a route from the routing decision only.
    The vehicle leaves the depot just in time for the first pickup (no waiting at the
    first stop); afterwards service starts at max(arrival, ready time).
    Without time windows service simply starts on arrival (vehicle leaves at t = 0)."""
    path = [0] + route + [inst.end]
    rows = []
    if use_tw and route:
        t_dep = max(inst.e[0], inst.e[route[0]] - inst.d[0, route[0]])
    else:
        t_dep = inst.e[0] if use_tw else 0.0
    t = t_dep
    for a, b in zip(path[:-1], path[1:]):
        arr = t + inst.d[a, b]
        start = max(arr, inst.e[b]) if use_tw else arr
        rows.append((b, arr, start, start - arr))
        t = start + inst.s[b]
    return t_dep, rows


def node_table(inst, routes, R=None, use_tw=True):
    """One row per visited node: arrival, service start, waiting time, time window,
    cumulative distance and remaining range (distance budget) on arrival."""
    rows = []
    for k, route in enumerate(routes):
        if not route:
            continue
        t_dep, sched = schedule(inst, route, use_tw)
        cum = 0.0
        rows.append(dict(vehicle=k + 1, stop=0, node="0 (depot)", arrival=np.nan,
                         start=t_dep, wait=0.0,
                         tw=f"[{inst.e[0]:.0f},{inst.l[0]:.0f}]" if use_tw else "-",
                         cum_dist=0.0, rem_range=R if R is not None else np.nan, tw_ok=True))
        prev = 0
        for pos, (node, arr, start, wait) in enumerate(sched, start=1):
            cum += inst.d[prev, node]
            prev = node
            label = f"{node}" if node != inst.end else "0 (depot)"
            ok = (start <= inst.l[node] + 1e-6) if use_tw else True
            rows.append(dict(vehicle=k + 1, stop=pos, node=label, arrival=arr, start=start,
                             wait=wait,
                             tw=f"[{inst.e[node]:.0f},{inst.l[node]:.0f}]" if use_tw else "-",
                             cum_dist=cum, rem_range=(R - cum) if R is not None else np.nan,
                             tw_ok=ok))
    df = pd.DataFrame(rows)
    return df.round(2)


def vehicle_table(inst, routes, R=None, use_tw=True):
    rows = []
    for k, route in enumerate(routes):
        if not route:
            rows.append(dict(vehicle=k + 1, used=False, route="0-0", n_stops=0,
                             distance=0.0, rem_range_end=R if R is not None else np.nan,
                             depart=np.nan, back=np.nan))
            continue
        t_dep, sched = schedule(inst, route, use_tw)
        dist = route_distance(inst, route)
        rows.append(dict(vehicle=k + 1, used=True,
                         route="0-" + "-".join(map(str, route)) + "-0",
                         n_stops=len(route), distance=dist,
                         rem_range_end=(R - dist) if R is not None else np.nan,
                         depart=t_dep, back=sched[-1][1]))
    return pd.DataFrame(rows).round(2)


def check_solution(inst, routes, R=None, use_tw=True):
    """Independent feasibility check of a solution (not using the solver's T/D values)."""
    visited = sorted(c for r in routes for c in r)
    assert visited == inst.C, f"not every location visited exactly once: {visited}"
    for r in routes:
        if not r:
            continue
        if R is not None:
            assert route_distance(inst, r) <= R + 1e-6, "range violated"
        if use_tw:
            _, sched = schedule(inst, r, True)
            for node, arr, start, _ in sched:
                assert start <= inst.l[node] + 1e-6, f"time window violated at {node}"
    return True


# --------------------------------------------------------------------------------------
# Plotting
# --------------------------------------------------------------------------------------
COLORS = ["#1f77b4", "#d62728", "#2ca02c", "#ff7f0e", "#9467bd", "#8c564b",
          "#e377c2", "#17becf", "#7f7f7f", "#bcbd22"]


def plot_routes(inst, routes, title, fname, show_tw=True):
    fig, ax = plt.subplots(figsize=(7.5, 7))
    xy = inst.xy
    ax.scatter(xy[1:inst.end, 0], xy[1:inst.end, 1], s=60, c="white", edgecolors="black",
               zorder=3)
    ax.scatter(xy[0, 0], xy[0, 1], s=180, marker="s", c="black", zorder=4, label="depot (0)")
    for i in inst.C:
        lbl = f"{i}"
        if show_tw and inst.tw:
            lbl += f"\n[{inst.e[i]:.0f},{inst.l[i]:.0f}]"
        ax.annotate(lbl, (xy[i, 0], xy[i, 1]), textcoords="offset points", xytext=(6, 4),
                    fontsize=8)
    for k, route in enumerate(routes):
        if not route:
            continue
        path = [0] + route + [inst.end]
        col = COLORS[k % len(COLORS)]
        for a, b in zip(path[:-1], path[1:]):
            ax.annotate("", xy=(xy[b, 0], xy[b, 1]), xytext=(xy[a, 0], xy[a, 1]),
                        arrowprops=dict(arrowstyle="-|>", color=col, lw=1.8,
                                        shrinkA=6, shrinkB=6))
        ax.plot([], [], color=col, lw=2,
                label=f"vehicle {k + 1}: {route_distance(inst, route):.1f}")
    ax.set_title(title, fontsize=11)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_aspect("equal")
    ax.grid(alpha=0.3)
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, fname), dpi=150)
    return fig


def header(txt):
    print("\n" + "=" * 100)
    print(txt)
    print("=" * 100)
