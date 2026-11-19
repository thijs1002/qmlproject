"""
ME44207 Quantitative Methods for Logistics - Group Assignment Q2 (2026)
Group number : <GROUP NUMBER>
Members      : <NAME 1>, <NAME 2>, <NAME 3>, <NAME 4>

Part 1b - VRP(TW) with a deterministic range R = mean of the 12 monthly ranges.
Solved with 3 vehicles, with and without time windows.
"""
import matplotlib.pyplot as plt
from gurobipy import GRB
import vrp_common as vc

K = 3                    # number of vehicles
TIME_LIMIT = 1800        # seconds (assignment: 30 min)
SHOW_PLOTS = True
VERBOSE = False          # True -> show the Gurobi log


def solve_deterministic(inst, K, R, label):
    m, x, D, T = vc.build_routing_model(inst, K, f"VRP_det_{label}", Dmax=R)
    # (C9) range constraint: D[n+1,k] <= R (also enforced by the upper bound of D)
    for k in range(K):
        m.addConstr(D[inst.end, k] <= R, name=f"range[{k}]")
    m.setObjective(vc.total_distance_expr(inst, x, K), GRB.MINIMIZE)
    vc.set_params(m, TIME_LIMIT, VERBOSE)
    m.optimize()
    return m, x


def report(inst, m, x, K, R, label, fig_name):
    vc.header(f"PART 1b | {label}")
    print(f"Status: {vc.STATUS_TXT.get(m.Status, m.Status)} | runtime {m.Runtime:.2f} s | "
          f"#vars {m.NumVars} (#binary {m.NumBinVars}) | #constraints {m.NumConstrs}")
    if m.SolCount == 0:
        print("No feasible solution found.")
        return None
    if m.Status != GRB.OPTIMAL:
        print(f"Optimality gap: {100 * m.MIPGap:.2f} %")
    routes = vc.extract_routes(inst, x, K)
    vc.check_solution(inst, routes, R, inst.tw)
    print(f"Total distance travelled: {m.ObjVal:.2f}   (range per vehicle R = {R:.2f})")
    print("\nRoutes per vehicle:")
    print(vc.vehicle_table(inst, routes, R, inst.tw).to_string(index=False))
    print("\nTime of visit and remaining range at each location:")
    print(vc.node_table(inst, routes, R, inst.tw).to_string(index=False))
    vc.plot_routes(inst, routes, f"Part 1b - {label}\ntotal distance {m.ObjVal:.2f}", fig_name)
    return dict(case=label, status=vc.STATUS_TXT.get(m.Status, m.Status),
                distance=round(m.ObjVal, 2), vehicles_used=sum(1 for r in routes if r),
                runtime_s=round(m.Runtime, 3), gap_pct=round(100 * (m.MIPGap or 0), 3))


if __name__ == "__main__":
    geo = vc.load_geo("geo_data.txt")
    rng = vc.load_range("range_data.txt")
    R = rng["RANGE"].mean()
    print(f"Deterministic range = mean of monthly ranges = {R:.4f}")

    summary = []
    for tw in (True, False):
        inst = vc.Instance(geo, time_windows=tw)
        label = "with time windows" if tw else "without time windows"
        m, x = solve_deterministic(inst, K, R, label)
        res = report(inst, m, x, K, R, label, f"part1b_{'tw' if tw else 'notw'}.png")
        if res:
            summary.append(res)

    vc.header("PART 1b | comparison")
    import pandas as pd
    print(pd.DataFrame(summary).to_string(index=False))
    if SHOW_PLOTS:
        plt.show()
