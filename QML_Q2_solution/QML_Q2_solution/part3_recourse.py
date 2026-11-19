"""
ME44207 Quantitative Methods for Logistics - Group Assignment Q2 (2026)
Group number : <GROUP NUMBER>
Members      : <NAME 1>, <NAME 2>, <NAME 3>, <NAME 4>

Part 3b / 3c - two-stage stochastic VRP(TW) with recourse "emergency charging".

First stage : the routes x[i,j,k] (one fixed plan for the whole year)
Second stage: for every month s (equally likely scenario, p_s = 1/12) the number of
              emergency recharges n[k,s] of vehicle k.
Because charging is possible anywhere, restores the range to R_s and costs a fixed penalty P,
a route of length L_k needs exactly ceil(L_k / R_s) - 1 recharges in month s, modelled as
              R_s * (1 + n[k,s]) >= L_k ,  n[k,s] integer >= 0
Objective   : min  sum_ij d_ij x_ijk  +  sum_s p_s * P * sum_k n[k,s]
"""
import math
import pandas as pd
import matplotlib.pyplot as plt
from gurobipy import GRB
import gurobipy as gp
import vrp_common as vc

TIME_LIMIT = 1800
SHOW_PLOTS = True
VERBOSE = False

# ---- experiment settings --------------------------------------------------------------
K_3B = 3                          # same fleet as in part 1
PENALTIES = [2, 60, 150]          # penalty per emergency recharge (distance units)
K_3C = 5                          # fleet for the large instance (see minimum_fleet())
RUN_3B, RUN_3C = True, True


def solve_recourse(inst, K, ranges, P, name):
    Ks, S = range(K), range(len(ranges))
    p = 1.0 / len(ranges)                         # equally likely monthly scenarios
    Dmax = inst.route_length_upper_bound()        # range is no longer a hard constraint
    m, x, D, T = vc.build_routing_model(inst, K, name, Dmax=Dmax)
    nmax = math.ceil(Dmax / min(ranges))
    n = m.addVars([(k, s) for k in Ks for s in S], vtype=GRB.INTEGER, lb=0, ub=nmax, name="n")
    for k in Ks:
        L_k = vc.route_length_expr(inst, x, k)
        used_k = 1 - x[0, inst.end, k]           # 1 if vehicle k leaves the depot
        for s in S:
            # (R1) recharges needed by vehicle k in scenario s:  R_s (u_k + n_ks) >= L_k
            #      (u_k instead of 1 is equivalent for integer x but gives a tighter LP)
            m.addConstr(ranges[s] * (used_k + n[k, s]) >= L_k, name=f"recharge[{k},{s}]")
    first_stage = vc.total_distance_expr(inst, x, K)
    second_stage = gp.quicksum(p * P * n[k, s] for k in Ks for s in S)
    m.setObjective(first_stage + second_stage, GRB.MINIMIZE)
    vc.set_params(m, TIME_LIMIT, VERBOSE)
    m.optimize()
    return m, x, n


def simulate_emergency(inst, route, R):
    """Drive the route with range R; recharge (to R) exactly when the battery is empty.
    Returns a list of (from_node, to_node, cumulative_km_at_recharge)."""
    path = [0] + route + [inst.end]
    battery, cum, events = R, 0.0, []
    for a, b in zip(path[:-1], path[1:]):
        dist = inst.d[a, b]
        while battery < dist - 1e-9:
            cum += battery
            dist -= battery
            events.append((a, 0 if b == inst.end else b, round(cum, 2)))
            battery = R
        battery -= dist
        cum += dist
    return events


def report_recourse(inst, m, x, n, K, ranges, months, P, label, fig_name):
    vc.header(f"{label} | penalty P = {P}")
    st = vc.STATUS_TXT.get(m.Status, m.Status)
    print(f"Status: {st} | runtime {m.Runtime:.2f} s | #vars {m.NumVars} "
          f"(#int {m.NumIntVars}) | #constraints {m.NumConstrs}")
    if m.SolCount == 0:
        print("No feasible solution.")
        return dict(case=label, P=P, status=st, runtime_s=round(m.Runtime, 2))
    if m.Status != GRB.OPTIMAL:
        print(f"Optimality gap: {100 * m.MIPGap:.2f} %")
    routes = vc.extract_routes(inst, x, K)
    vc.check_solution(inst, routes, None, inst.tw)
    dist = sum(vc.route_distance(inst, r) for r in routes)

    # second stage: recharges per vehicle and month (model value and simulation must agree)
    tab = pd.DataFrame(index=months)
    where = []
    for k, r in enumerate(routes):
        if not r:
            continue
        col = []
        for s, Rs in enumerate(ranges):
            ev = simulate_emergency(inst, r, Rs)
            if len(ev) != round(n[k, s].X):
                print(f"WARNING: model n[{k},{s}] = {n[k, s].X:.0f} but simulation needs "
                      f"{len(ev)} recharges (vehicle {k + 1}, {months[s]})")
            col.append(len(ev))
            for a, b, km in ev:
                where.append(dict(month=months[s], vehicle=k + 1, R_s=Rs,
                                  on_leg=f"{a}->{b}", at_km=km))
        tab[f"veh {k + 1} (L={vc.route_distance(inst, r):.1f})"] = col
    tab["R_s"] = ranges
    tab["total"] = tab.drop(columns="R_s").sum(axis=1)
    tab["penalty"] = tab["total"] * P
    total_ch = int(tab["total"].sum())
    exp_pen = P * total_ch / len(ranges)

    print(f"Total distance (first stage) : {dist:.2f}")
    print(f"Emergency recharges over the 12 monthly scenarios: {total_ch}  "
          f"-> expected per day: {total_ch / len(ranges):.3f}")
    print(f"Expected penalty per day     : {exp_pen:.2f}")
    print(f"Objective (expected cost)    : {m.ObjVal:.2f}")
    print("\nFirst-stage routes:")
    print(vc.vehicle_table(inst, routes, None, inst.tw)
          .drop(columns=["rem_range_end"]).to_string(index=False))
    print("\nEmergency recharges per vehicle and month:")
    print(tab.to_string())
    if where:
        print("\nWhere the recharges happen (leg and cumulative km):")
        print(pd.DataFrame(where).to_string(index=False))
    vc.plot_routes(inst, routes, f"{label} - P = {P}\ndistance {dist:.2f}, "
                   f"{total_ch} recharges/12 months, E[cost] {m.ObjVal:.2f}", fig_name)
    return dict(case=label, P=P, status=st, distance=round(dist, 2),
                vehicles_used=sum(1 for r in routes if r),
                recharges_12_months=total_ch, exp_recharges_per_day=round(total_ch / 12, 3),
                exp_penalty_per_day=round(exp_pen, 2), objective=round(m.ObjVal, 2),
                runtime_s=round(m.Runtime, 3), gap_pct=round(100 * (m.MIPGap or 0), 3),
                routes=" | ".join("0-" + "-".join(map(str, r)) + "-0" for r in routes if r))


def minimum_fleet(inst, kmax=12):
    """Smallest number of vehicles for which the time windows can be met (range ignored)."""
    for K in range(1, kmax + 1):
        m, x, D, T = vc.build_routing_model(inst, K, f"minfleet_{K}",
                                            Dmax=inst.route_length_upper_bound())
        m.setObjective(vc.total_distance_expr(inst, x, K), GRB.MINIMIZE)
        vc.set_params(m, 300, False)
        m.Params.SolutionLimit = 1                 # feasibility check only
        m.optimize()
        if m.SolCount > 0:
            return K
    return None


if __name__ == "__main__":
    rng = vc.load_range("range_data.txt")
    ranges, months = list(rng["RANGE"]), list(rng["MONTH"])
    summary = []

    if RUN_3B:
        geo = vc.load_geo("geo_data.txt")
        for tw in (True, False):
            inst = vc.Instance(geo, time_windows=tw)
            twl = "TW" if tw else "noTW"
            for P in PENALTIES:
                m, x, n = solve_recourse(inst, K_3B, ranges, P, f"rec_{twl}_{P}")
                label = f"PART 3b | {'with' if tw else 'without'} time windows | K = {K_3B}"
                summary.append(report_recourse(inst, m, x, n, K_3B, ranges, months, P, label,
                                               f"part3b_{twl}_P{P}.png"))

    if RUN_3C:
        geo_l = vc.load_geo("geo_data_large.txt")
        inst = vc.Instance(geo_l, time_windows=True)
        kmin = minimum_fleet(inst)
        vc.header(f"PART 3c | large instance: minimum fleet for time-window feasibility = {kmin}")
        for P in PENALTIES:
            m, x, n = solve_recourse(inst, K_3C, ranges, P, f"rec_large_{P}")
            label = f"PART 3c | large instance, with time windows | K = {K_3C}"
            summary.append(report_recourse(inst, m, x, n, K_3C, ranges, months, P, label,
                                           f"part3c_large_P{P}.png"))

    vc.header("PART 3b/3c | summary")
    print(pd.DataFrame(summary).to_string(index=False))
    if SHOW_PLOTS:
        plt.show()
