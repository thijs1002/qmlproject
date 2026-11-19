"""
ME44207 Quantitative Methods for Logistics - Group Assignment Q2 (2026)
Group number : <GROUP NUMBER>
Members      : <NAME 1>, <NAME 2>, <NAME 3>, <NAME 4>

Part 3d - new recourse strategy: PLANNED DEPOT RECHARGING DURING IDLE TIME
(with emergency charging kept as last-resort backup), case with time windows.

Idea: because of the time windows the vehicles have long idle (waiting) periods. In a month
with a low range a vehicle can use such a gap to make a detour via the depot between two
pickups i -> 0 -> j and fast-charge there to full range (charging time TAU). This costs the
extra distance d_i0 + d_0j - d_ij plus a small handling cost F, and is only possible if the
time windows are still met in that month. If no detour fits, emergency charging (penalty P,
as in part 3b) is still available -> complete recourse.

First stage : routes x[i,j,k]  (fixed for the whole year)
Second stage (per month s): depot detours z[i,j,s], emergency charges w[i,j,s],
              remaining range b[i,s] on arrival and service start times Ts[i,s].
"""
import pandas as pd
import matplotlib.pyplot as plt
import gurobipy as gp
from gurobipy import GRB
import vrp_common as vc

TIME_LIMIT = 1800
SHOW_PLOTS = True
VERBOSE = False

K = 3                       # same fleet as in part 3b
PENALTIES = [2, 60, 150]    # emergency-charging penalty (same values as part 3b)
TAU = 30                    # depot fast-charging time (time units) -> full range
F_DEPOT = 5                 # fixed handling cost per depot recharge (distance units)
GEO_FILE = "geo_data.txt"
RUN_TAU_SENSITIVITY = True
P_SENS, TAUS = 60, [0, 90, 180]   # extra runs: depot charging time 0 / 90 / 180


def solve_depot_recourse(inst, K, ranges, P, tau, F, name, allow_depot=True):
    Ks, S = range(K), range(len(ranges))
    p = 1.0 / len(ranges)
    A, C, end, d, e, l, sv = inst.arcs, inst.C, inst.end, inst.d, inst.e, inst.l, inst.s
    m, x, D, T = vc.build_routing_model(inst, K, name, Dmax=inst.route_length_upper_bound())
    y = {(i, j): gp.quicksum(x[i, j, k] for k in Ks) for (i, j) in A}   # arc used by any vehicle

    # arcs on which a depot detour is possible (time-window pre-processing)
    Z = [(i, j) for (i, j) in A if i in C and j in C
         and e[i] + sv[i] + d[i, 0] + tau + d[0, j] <= l[j] + 1e-9] if allow_depot else []
    W = [(i, j) for (i, j) in A if (i, j) != (0, end)]                # emergency possible

    z = m.addVars([(i, j, s) for (i, j) in Z for s in S], vtype=GRB.BINARY, name="z")
    w = m.addVars([(i, j, s) for (i, j) in W for s in S], vtype=GRB.BINARY, name="w")
    b = m.addVars([(i, s) for i in C for s in S], lb=0.0,
                  ub={(i, s): ranges[s] for i in C for s in S}, name="b")
    Ts = m.addVars([(i, s) for i in C for s in S],
                   lb={(i, s): e[i] for i in C for s in S},
                   ub={(i, s): l[i] for i in C for s in S}, name="Ts")

    for s in S:
        Rs = ranges[s]
        for (i, j) in A:
            zij = z[i, j, s] if (i, j) in Z else 0
            wij = w[i, j, s] if (i, j) in W else 0
            if (i, j) == (0, end):
                continue
            # (S1) recourse only on arcs that are used, at most one action per arc
            m.addConstr(zij + wij <= y[i, j], name=f"act[{i},{j},{s}]")
            if i == 0:                     # leaving the depot with a full battery
                m.addConstr(b[j, s] <= Rs - d[0, j] + Rs * wij + d[0, j] * (1 - y[i, j]),
                            name=f"bat0[{j},{s}]")
                m.addConstr(Ts[j, s] >= e[0] + d[0, j] - (e[0] + d[0, j]) * (1 - y[i, j]),
                            name=f"t0[{j},{s}]")
            elif j == end:                 # must be able to drive back to the depot
                m.addConstr(b[i, s] + Rs * wij >= d[i, end] * y[i, j], name=f"batE[{i},{s}]")
                Mt = max(0.0, l[i] + sv[i] + d[i, end] - l[end])
                m.addConstr(Ts[i, s] + sv[i] + d[i, end] <= l[end] + Mt * (1 - y[i, j]),
                            name=f"tE[{i},{s}]")
            else:                          # between two pickups
                # (S2) direct drive (relaxed if a depot detour is made)
                m.addConstr(b[j, s] <= b[i, s] - d[i, j] + Rs * wij
                            + (Rs + d[i, j]) * (1 - y[i, j] + zij), name=f"bat[{i},{j},{s}]")
                Mt = max(0.0, l[i] + sv[i] + d[i, j] - e[j])
                m.addConstr(Ts[j, s] >= Ts[i, s] + sv[i] + d[i, j] - Mt * (1 - y[i, j]),
                            name=f"t[{i},{j},{s}]")
                if (i, j) in Z:
                    # (S3) depot detour: reach the depot, recharge to R_s, drive to j
                    m.addConstr(b[i, s] >= d[i, 0] * zij, name=f"reach0[{i},{j},{s}]")
                    m.addConstr(b[j, s] <= Rs - d[0, j] + d[0, j] * (1 - zij),
                                name=f"batZ[{i},{j},{s}]")
                    detour_t = sv[i] + d[i, 0] + tau + d[0, j]
                    Mz = max(0.0, l[i] + detour_t - e[j])
                    m.addConstr(Ts[j, s] >= Ts[i, s] + detour_t - Mz * (1 - zij),
                                name=f"tZ[{i},{j},{s}]")

    first = vc.total_distance_expr(inst, x, K)
    second = gp.quicksum(p * (F + d[i, 0] + d[0, j] - d[i, j]) * z[i, j, s]
                         for (i, j) in Z for s in S) \
        + gp.quicksum(p * P * w[i, j, s] for (i, j) in W for s in S)
    m.setObjective(first + second, GRB.MINIMIZE)
    vc.set_params(m, TIME_LIMIT, VERBOSE)
    m.optimize()
    return m, x, z, w, Z, W


def simulate(inst, route, Rs, detours, emergencies, tau):
    """Independent check of one route in one month given the recourse decisions:
    returns (battery_ok, time_windows_ok, timetable rows)."""
    path = [0] + route + [inst.end]
    bat = Rs
    rows = []
    t = max(inst.e[0], inst.e[route[0]] - inst.d[0, route[0]])
    ok_b = ok_t = True
    for a, c in zip(path[:-1], path[1:]):
        if (a, c) in detours:
            bat -= inst.d[a, 0]
            ok_b &= bat >= -1e-6
            t += inst.d[a, 0] + tau
            bat = Rs - inst.d[0, c]
            t += inst.d[0, c]
            act = "depot recharge"
        elif (a, c) in emergencies:
            bat = min(Rs, bat - inst.d[a, c] + Rs)      # recharge to R_s when empty
            t += inst.d[a, c]
            act = "emergency"
        else:
            bat -= inst.d[a, c]
            t += inst.d[a, c]
            act = ""
        ok_b &= bat >= -1e-6
        start = max(t, inst.e[c])
        ok_t &= start <= inst.l[c] + 1e-6
        rows.append(dict(node=c if c != inst.end else 0, action_on_leg=act,
                         arrival=round(t, 2), start=round(start, 2), tw=f"[{inst.e[c]:.0f},"
                         f"{inst.l[c]:.0f}]", rem_range=round(bat, 2)))
        t = start + inst.s[c]
    return ok_b, ok_t, rows


def report(inst, m, x, z, w, Z, W, K, ranges, months, P, tau, F, label, fig_name):
    vc.header(f"{label} | P = {P}, TAU = {tau}, F = {F}")
    st = vc.STATUS_TXT.get(m.Status, m.Status)
    print(f"Status: {st} | runtime {m.Runtime:.2f} s | #vars {m.NumVars} "
          f"(#int {m.NumIntVars}) | #constraints {m.NumConstrs}")
    if m.SolCount == 0:
        return dict(P=P, status=st)
    if m.Status != GRB.OPTIMAL:
        print(f"Optimality gap: {100 * m.MIPGap:.2f} %")
    routes = vc.extract_routes(inst, x, K)
    vc.check_solution(inst, routes, None, True)
    dist = sum(vc.route_distance(inst, r) for r in routes)
    print(f"Total distance (first stage): {dist:.2f}")
    print(vc.vehicle_table(inst, routes, None, True)
          .drop(columns=["rem_range_end"]).to_string(index=False))

    rows, n_dep, n_em, km_det = [], 0, 0, 0.0
    example = None
    for s, Rs in enumerate(ranges):
        det = {(i, j) for (i, j) in Z if z[i, j, s].X > 0.5}
        em = {(i, j) for (i, j) in W if w[i, j, s].X > 0.5}
        for k, r in enumerate(routes):
            if not r:
                continue
            path = [0] + r + [inst.end]
            legs = set(zip(path[:-1], path[1:]))
            dk, ek = det & legs, em & legs
            ok_b, ok_t, tt = simulate(inst, r, Rs, dk, ek, tau)
            assert ok_b and ok_t, f"recourse infeasible in month {months[s]}, vehicle {k + 1}"
            extra = sum(inst.d[i, 0] + inst.d[0, j] - inst.d[i, j] for (i, j) in dk)
            n_dep += len(dk)
            n_em += len(ek)
            km_det += extra
            rows.append(dict(month=months[s], R_s=Rs, vehicle=k + 1,
                             depot_recharges=", ".join(f"{i}->0->{j}" for i, j in sorted(dk)),
                             emergency=", ".join(f"{i}->{j if j != inst.end else 0}"
                                                 for i, j in sorted(ek)),
                             extra_km=round(extra, 2)))
            if dk and example is None:
                example = (months[s], k + 1, tt)
    tab = pd.DataFrame(rows)
    tab = tab[(tab.depot_recharges != "") | (tab.emergency != "")]
    exp_cost = (F * n_dep + km_det + P * n_em) / len(ranges)
    print(f"\nDepot recharges over 12 months: {n_dep} (extra km {km_det:.2f}) | "
          f"emergency recharges: {n_em}")
    print(f"Expected recourse cost per day: {exp_cost:.2f} | objective {m.ObjVal:.2f}")
    if len(tab):
        print("\nRecourse actions per month and vehicle:")
        print(tab.to_string(index=False))
    if example:
        print(f"\nExample timetable with a depot recharge ({example[0]}, vehicle {example[1]}):")
        print(pd.DataFrame(example[2]).to_string(index=False))
    vc.plot_routes(inst, routes, f"{label} - P = {P}\ndistance {dist:.2f}, {n_dep} depot + "
                   f"{n_em} emergency recharges / 12 months", fig_name)
    return dict(P=P, tau=tau, F=F, status=st, distance=round(dist, 2),
                vehicles_used=sum(1 for r in routes if r), depot_recharges_12m=n_dep,
                detour_km_12m=round(km_det, 2), emergency_12m=n_em,
                exp_recourse_cost=round(exp_cost, 2), objective=round(m.ObjVal, 2),
                runtime_s=round(m.Runtime, 3), gap_pct=round(100 * (m.MIPGap or 0), 3),
                routes=" | ".join("0-" + "-".join(map(str, r)) + "-0" for r in routes if r))


if __name__ == "__main__":
    rng = vc.load_range("range_data.txt")
    ranges, months = list(rng["RANGE"]), list(rng["MONTH"])
    inst = vc.Instance(vc.load_geo(GEO_FILE), time_windows=True)
    summary = []
    for P in PENALTIES:
        m, x, z, w, Z, W = solve_depot_recourse(inst, K, ranges, P, TAU, F_DEPOT, f"dep_{P}")
        summary.append(report(inst, m, x, z, w, Z, W, K, ranges, months, P, TAU, F_DEPOT,
                              f"PART 3d | depot recharging, with time windows, K = {K}",
                              f"part3d_P{P}.png"))
    vc.header("PART 3d | summary")
    print(pd.DataFrame(summary).to_string(index=False))

    # sensitivity: how much does the depot charging time matter? (P fixed)
    if RUN_TAU_SENSITIVITY:
        sens = []
        for tau in TAUS:
            m, x, z, w, Z, W = solve_depot_recourse(inst, K, ranges, P_SENS, tau, F_DEPOT,
                                                    f"dep_tau{tau}")
            sens.append(report(inst, m, x, z, w, Z, W, K, ranges, months, P_SENS, tau,
                               F_DEPOT, "PART 3d | sensitivity charging time",
                               f"part3d_P{P_SENS}_tau{tau}.png"))
        vc.header(f"PART 3d | sensitivity to depot charging time TAU (P = {P_SENS})")
        print(pd.DataFrame(sens).to_string(index=False))
    if SHOW_PLOTS:
        plt.show()
