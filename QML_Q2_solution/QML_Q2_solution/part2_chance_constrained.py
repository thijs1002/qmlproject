"""
ME44207 Quantitative Methods for Logistics - Group Assignment Q2 (2026)
Group number : <GROUP NUMBER>
Members      : <NAME 1>, <NAME 2>, <NAME 3>, <NAME 4>

Part 2b - chance-constrained VRP(TW).
Range R ~ N(mu, sigma^2), mu and sigma estimated from the 12 monthly ranges.
Chance constraint per vehicle k:  P( L_k <= R ) >= alpha
Deterministic equivalent:          L_k <= mu + sigma * Phi^-1(1 - alpha) = mu - z_alpha * sigma
"""
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import norm
from gurobipy import GRB
import vrp_common as vc

K = 3
ALPHAS = [0.5, 0.75, 0.90, 0.99]
TIME_LIMIT = 1800
SHOW_PLOTS = True
VERBOSE = False
EXTRA_FLEET_IF_INFEASIBLE = True   # if infeasible with K vehicles, also find the minimum fleet


def range_quantile(mu, sigma, alpha):
    """Deterministic equivalent right-hand side R_alpha = mu - z_alpha * sigma."""
    return mu + sigma * norm.ppf(1 - alpha)


def solve_cc(inst, K, R_alpha, label):
    m, x, D, T = vc.build_routing_model(inst, K, f"VRP_cc_{label}", Dmax=R_alpha)
    for k in range(K):
        # (C9') deterministic equivalent of the chance constraint
        m.addConstr(D[inst.end, k] <= R_alpha, name=f"cc_range[{k}]")
    m.setObjective(vc.total_distance_expr(inst, x, K), GRB.MINIMIZE)
    vc.set_params(m, TIME_LIMIT, VERBOSE)
    m.optimize()
    return m, x


if __name__ == "__main__":
    geo = vc.load_geo("geo_data.txt")
    rng = vc.load_range("range_data.txt")
    mu = rng["RANGE"].mean()
    sigma = rng["RANGE"].std(ddof=1)       # unbiased sample estimate (n - 1)
    vc.header("PART 2 | estimated range distribution")
    print(f"mu = {mu:.4f}   sigma (sample, ddof=1) = {sigma:.4f}   "
          f"[population sigma (ddof=0) = {rng['RANGE'].std(ddof=0):.4f}]")
    tab = pd.DataFrame({"alpha": ALPHAS,
                        "z_alpha": [norm.ppf(a) for a in ALPHAS],
                        "R_alpha": [range_quantile(mu, sigma, a) for a in ALPHAS]})
    tab["months_with_range_below_R_alpha"] = [int((rng["RANGE"] < r).sum()) for r in tab.R_alpha]
    print(tab.round(3).to_string(index=False))

    summary = []
    for tw in (True, False):
        inst = vc.Instance(geo, time_windows=tw)
        twl = "TW" if tw else "noTW"
        for a in ALPHAS:
            Ra = range_quantile(mu, sigma, a)
            label = f"alpha={a} | {'with' if tw else 'without'} time windows"
            m, x = solve_cc(inst, K, Ra, f"{twl}_{a}")
            vc.header(f"PART 2b | {label} | R_alpha = {Ra:.2f}")
            st = vc.STATUS_TXT.get(m.Status, m.Status)
            print(f"Status: {st} | runtime {m.Runtime:.2f} s")
            row = dict(time_windows=tw, alpha=a, R_alpha=round(Ra, 2), status=st,
                       distance=None, vehicles_used=None, runtime_s=round(m.Runtime, 3),
                       gap_pct=None, routes="")
            if m.SolCount > 0:
                routes = vc.extract_routes(inst, x, K)
                vc.check_solution(inst, routes, Ra, tw)
                print(f"Total distance travelled: {m.ObjVal:.2f}")
                if m.Status != GRB.OPTIMAL:
                    print(f"Optimality gap: {100 * m.MIPGap:.2f} %")
                print("\nRoutes per vehicle (remaining budget w.r.t. R_alpha):")
                print(vc.vehicle_table(inst, routes, Ra, tw).to_string(index=False))
                print("\nTime of visit and remaining range at each location:")
                print(vc.node_table(inst, routes, Ra, tw).to_string(index=False))
                vc.plot_routes(inst, routes,
                               f"Part 2b - {label}\nR_alpha = {Ra:.2f}, distance {m.ObjVal:.2f}",
                               f"part2b_{twl}_alpha{int(a * 100)}.png")
                row.update(distance=round(m.ObjVal, 2),
                           vehicles_used=sum(1 for r in routes if r),
                           gap_pct=round(100 * (m.MIPGap or 0), 3),
                           routes=" | ".join("0-" + "-".join(map(str, r)) + "-0"
                                             for r in routes if r))
            else:
                print(f"No feasible solution with {K} vehicles and R_alpha = {Ra:.2f}.")
            summary.append(row)

            # extra experiment: smallest fleet for which this confidence level is attainable
            if m.Status == GRB.INFEASIBLE and EXTRA_FLEET_IF_INFEASIBLE:
                for K2 in range(K + 1, inst.n + 1):
                    m2, x2 = solve_cc(inst, K2, Ra, f"{twl}_{a}_K{K2}")
                    if m2.SolCount > 0:
                        routes = vc.extract_routes(inst, x2, K2)
                        vc.check_solution(inst, routes, Ra, tw)
                        print(f"\n-> extra: feasible with K = {K2} vehicles: distance "
                              f"{m2.ObjVal:.2f}, runtime {m2.Runtime:.2f} s, "
                              f"status {vc.STATUS_TXT.get(m2.Status, m2.Status)}")
                        print(vc.vehicle_table(inst, routes, Ra, tw).to_string(index=False))
                        print(vc.node_table(inst, routes, Ra, tw).to_string(index=False))
                        vc.plot_routes(inst, routes,
                                       f"Part 2b - {label} - K = {K2}\n"
                                       f"R_alpha = {Ra:.2f}, distance {m2.ObjVal:.2f}",
                                       f"part2b_{twl}_alpha{int(a * 100)}_K{K2}.png")
                        summary.append(dict(time_windows=tw, alpha=a, R_alpha=round(Ra, 2),
                                            status=f"{vc.STATUS_TXT.get(m2.Status)} (K={K2})",
                                            distance=round(m2.ObjVal, 2),
                                            vehicles_used=sum(1 for r in routes if r),
                                            runtime_s=round(m2.Runtime, 3),
                                            gap_pct=round(100 * (m2.MIPGap or 0), 3),
                                            routes=" | ".join("0-" + "-".join(map(str, r)) + "-0"
                                                              for r in routes if r)))
                        break
                    print(f"   K = {K2}: {vc.STATUS_TXT.get(m2.Status, m2.Status)}")

    vc.header("PART 2b | summary over all experiments")
    print(pd.DataFrame(summary).to_string(index=False))
    if SHOW_PLOTS:
        plt.show()
