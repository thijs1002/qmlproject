# ME44207 Group Assignment Q2 – solution notes

These notes go with the Python/Gurobi files. They contain the mathematical models (compact, as the
assignment asks), the justification of the modelling choices (kept separate from the models), the
results we obtained and the points worth discussing in the report.

> **Important – computation times.** Gurobi could not be installed where this solution was
> prepared. Every model was therefore run through a small test layer that sends the *same* gurobipy
> model to the open-source solver HiGHS. Objective values and routes are solver-independent (all
> small-instance optima were also confirmed by an exact brute-force dynamic program), but **all
> run times below are HiGHS times**. Re-run every script in VS Code with Gurobi and report *those*
> run times (Gurobi is normally much faster). When several optimal solutions exist, Gurobi may
> show a different but equally good route set.

## 0. Files and how to run

| File | Content |
|---|---|
| `vrp_common.py` | data reading (files used as is), `Instance` (sets, arcs, distances), the shared three-index routing model, route extraction, tables, plots |
| `part1_deterministic.py` | 1b – deterministic range, 3 vehicles, with/without time windows |
| `part2_chance_constrained.py` | 2b – chance constraint for α ∈ {0.5, 0.75, 0.90, 0.99}, with/without time windows (+ minimum fleet when infeasible) |
| `part3_recourse.py` | 3b (small instance, P ∈ {2, 60, 150}, with/without time windows) and 3c (large instance, K = 5) |
| `part3d_depot_recharging.py` | 3d – new recourse: depot recharging during idle time, emergency charging as backup |

Put `geo_data.txt`, `geo_data_large.txt` and `range_data.txt` in the same folder as the `.py` files,
open the folder in VS Code and run a part script (▶ or `python part1_deterministic.py`).
Tables are printed to the terminal, figures are saved to `output/` and shown at the end
(`SHOW_PLOTS = False` switches that off; `VERBOSE = True` shows the Gurobi log).
Only `gurobipy`, `numpy`, `pandas`, `math`, `scipy` and `matplotlib` are used.
Use the free Gurobi **academic license**: the size-limited license that comes with
`pip install gurobipy` stops at 2,000 variables/constraints, and the 3c and 3d models are larger.
Fill in group number and names in the header of every `.py` file.

---

## Part 1a – mathematical model (deterministic range)

**Sets**
- $N=\{0,1,\dots,n\}$: locations as in the data file; $0$ = depot, $C=N\setminus\{0\}$ = pickup locations
- $V=N\cup\{n+1\}$: node $n+1$ is a copy of the depot (route end), created in Python after reading the file
- $A=\{(0,j):j\in C\cup\{n+1\}\}\cup\{(i,j):i\in C,\ j\in C\cup\{n+1\},\ i\neq j\}$: arcs
- $K=\{1,\dots,m\}$: vehicles

**Parameters**
- $d_{ij}$: Euclidean distance = travel time of arc $(i,j)$
- $[e_i,l_i]$: time window, $s_i$: service time ($e_{n+1}=e_0$, $l_{n+1}=l_0$, $s_0=s_{n+1}=0$)
- $R=\tfrac{1}{12}\sum_{t=1}^{12}R_t = 119.28$: vehicle range
- $M_{ij}=\max\{0,\ l_i+s_i+d_{ij}-e_j\}$, $\ M^D_{ij}=R+d_{ij}$: big-M constants

**Decision variables**
- $x_{ijk}\in\{0,1\}$: 1 if vehicle $k$ drives arc $(i,j)$
- $T_{ik}\ge 0$: start of service of vehicle $k$ at node $i$
- $D_{ik}\ge 0$: distance driven by vehicle $k$ on arrival at node $i$ (remaining range $=R-D_{ik}$)

**Model**

$$\min \sum_{k\in K}\sum_{(i,j)\in A} d_{ij}\,x_{ijk}$$

$$
\begin{aligned}
&\textstyle\sum_{k\in K}\sum_{j:(i,j)\in A} x_{ijk}=1 && \forall i\in C &(1)\\
&\textstyle\sum_{j:(0,j)\in A} x_{0jk}=1 && \forall k\in K &(2)\\
&\textstyle\sum_{i:(i,h)\in A} x_{ihk}-\sum_{j:(h,j)\in A} x_{hjk}=0 && \forall h\in C,\ k\in K &(3)\\
&\textstyle\sum_{i:(i,n+1)\in A} x_{i,n+1,k}=1 && \forall k\in K &(4)\\
&D_{0k}=0 && \forall k\in K &(5)\\
&D_{jk}\ge D_{ik}+d_{ij}-M^D_{ij}(1-x_{ijk}) && \forall (i,j)\in A,\ k\in K &(6)\\
&D_{n+1,k}\le R && \forall k\in K &(7)\\
&T_{jk}\ge T_{ik}+s_i+d_{ij}-M_{ij}(1-x_{ijk}) && \forall (i,j)\in A,\ k\in K &(8)\\
&e_i\le T_{ik}\le l_i && \forall i\in V,\ k\in K &(9)\\
&x_{ijk}\in\{0,1\},\ D_{ik}\ge 0,\ T_{ik}\ge 0 &&&(10)
\end{aligned}
$$

Without time windows: drop (8)–(9).

### Explanation / modelling choices (separate section in the report)
- (1) every location visited exactly once; (2)–(4) every vehicle leaves the depot once, flow is conserved, and it returns once. Arc $(0,n+1)$ with $d=0$ means "vehicle not used", so the model works for any $m$ and $n$.
- (6) accumulates distance along the route. It replaces subtour elimination because $d_{ij}>0$. (8) does the same with time when time windows are used.
- (7) is the range constraint: the full route length of every vehicle stays within $R$.
- Code-level strengthening (valid, does not change the optimum): arcs that can never meet a time window ($e_i+s_i+d_{ij}>l_j$) are removed; tighter bounds $d_{0i}\le D_{ik}\le R-d_{i0}$ with matching big-M; vehicles are numbered by non-increasing route length ($L_k\ge L_{k+1}$) to remove symmetric copies. This made the no-time-window case about 10× faster.
- Reported times are recomputed from the routes: the vehicle leaves the depot just in time for its first pickup, and service starts at $\max(\text{arrival}, e_j)$. Every solution is independently checked for feasibility in the code.

## Part 1b – results (3 vehicles, R = 119.28)

| Case | Total distance | Vehicles used | Routes (length) | Run time (HiGHS) |
|---|---|---|---|---|
| With time windows | **212.53** | 2 | 0-3-4-2-9-5-0 (106.91) · 0-10-7-6-8-1-0 (105.62) | 0.2 s |
| Without time windows | **190.80** | 2 | 0-10-7-3-4-2-9-5-0 (114.93) · 0-1-6-8-0 (75.87) | 11 s |

Remaining range at the return: 12.38 / 13.66 (with time windows) and 4.35 / 43.41 (without).
The full per-node table (arrival, start, waiting, time window, cumulative distance, remaining range) is printed by the script.

**Discussion points**
- Time windows cost 21.73 distance (+11.4 %). They force the visiting order (3→4→2/10→7→6/8/1→9/5), so the routes zig-zag through the area instead of forming compact loops.
- In both cases only 2 of the 3 vehicles are used. The range (not the time) is the binding resource: one vehicle could serve all 10 locations within the time windows, but that route is 208.38 > R.
- There is a lot of idle time with time windows (e.g. 598 min waiting before node 9, 173 min before node 7). That matters again in 3d.
- Run time: the time-window case solves much faster even though it is "more constrained". The time windows remove many arcs (pre-processing) and make the LP relaxation much tighter. Without time windows only the weak distance-based MTZ constraints (6) eliminate subtours, and many more route combinations are possible.

---

## Part 2a – chance-constrained model

**Changes to the Part 1 model**
- New parameters: $\mu=119.28$, $\sigma=17.28$ (sample standard deviation of the 12 monthly ranges), confidence level $\alpha$, $z_\alpha=\Phi^{-1}(\alpha)$, and $R_\alpha=\mu-z_\alpha\sigma$.
- Sets, variables, objective and constraints (1)–(6), (8)–(10) are unchanged.
- Constraint (7) becomes the chance constraint
  $$\mathbb{P}\left(D_{n+1,k}\le \tilde R\right)\ge \alpha\quad\forall k\in K,\qquad \tilde R\sim N(\mu,\sigma^2)$$
  with deterministic equivalent
  $$D_{n+1,k}\le \mu-z_\alpha\,\sigma = R_\alpha\qquad \forall k\in K \qquad (7')$$
  and $M^D_{ij}=R_\alpha+d_{ij}$.

**Justification**
- $\mathbb{P}(\tilde R\ge D)=1-\Phi\!\left(\frac{D-\mu}{\sigma}\right)\ge\alpha \iff \frac{D-\mu}{\sigma}\le\Phi^{-1}(1-\alpha)=-z_\alpha \iff D\le \mu-z_\alpha\sigma$. The randomness is only in the right-hand side and its distribution is known, so this is an exact equivalent (not an approximation), and the model stays a MILP. $\mathbb{P}(\text{exceed})\le 1-\alpha$ is the same statement.
- Individual vs joint: all vehicles face the same range on a given day, so $\mathbb{P}(\max_k D_{n+1,k}\le\tilde R)\ge\alpha \iff \max_k D_{n+1,k}\le R_\alpha$. The joint chance constraint gives exactly the same model as one constraint per vehicle.
- $\sigma$ uses $n-1$ (unbiased sample estimate, 12 observations). Population $\sigma$ would be 16.54, giving slightly less conservative $R_\alpha$.
- $\alpha=0.5$ gives $z_\alpha=0$, so $R_\alpha=\mu$: identical to Part 1.
- Caveat for the discussion: the 12 monthly values are seasonal (low in winter, high in summer), not really normal. For $\alpha=0.99$ the normal model gives $R_\alpha=79.09$, below the lowest month ever observed (89.35), so the normal tail is very conservative here.

## Part 2b – results (3 vehicles)

| α | $R_\alpha$ | With time windows: distance / vehicles | Without time windows: distance / vehicles |
|---|---|---|---|
| 0.50 | 119.28 | 212.53 / 2 (= Part 1) | 190.80 / 2 (= Part 1) |
| 0.75 | 107.63 | 212.53 / 2 (same routes, 106.91 ≤ 107.63) | 195.03 / 2 – 0-7-3-1-6-8-0 (99.14) · 0-10-4-2-9-5-0 (95.89) |
| 0.90 | 97.14 | 246.41 / 3 – 0-3-4-7-1-0 (91.80) · 0-10-2-9-5-0 (81.11) · 0-8-6-0 (73.51) | 227.34 / 3 – 0-7-3-4-2-10-0 (91.90) · 0-1-6-8-0 (75.87) · 0-5-9-0 (59.57) |
| 0.99 | 79.09 | **infeasible with 3 vehicles**; K = 4: 293.84 | **infeasible with 3 vehicles**; K = 4: 255.80 |

HiGHS run times: with time windows ≤ 1 s each; without time windows 2–42 s.
Proving infeasibility at α = 0.99 is fast (0.05 s / 1.9 s).

**Discussion points**
- A higher α means a smaller distance budget $R_\alpha$. The cost rises in steps: nothing changes as long as the current routes still fit (α = 0.75 with time windows), then the solution jumps to more vehicles (α = 0.90 needs all 3).
- α = 0.99 cannot be met with 3 vehicles: the routes would have to be ≤ 79.09, and 4 vehicles are needed (+38 % distance with time windows).
- Check against the data (useful for the report): with the 12 real months, the α = 0.5 plan (= Part 1) exceeds the range in 3 of 12 months with time windows (Jan, Feb, Dec) and 5 of 12 without. The α = 0.90 plans fail only in January (1/12 = 8.3 %). The α = 0.99 plans never fail.
- Run times rise with α in the no-time-window case: a tighter budget makes the problem harder to prove optimal.

---

## Part 3a – two-stage stochastic model with recourse (emergency charging)

**New sets / parameters**: $S$ = months (scenarios), $p_s=1/|S|=1/12$, $R_s$ = range in month $s$, $P$ = penalty per emergency recharge.

**Variables**: first stage $x_{ijk}, T_{ik}, D_{ik}$ as in Part 1. Second stage $n_{ks}\in\mathbb{Z}_{\ge0}$ = number of emergency recharges of vehicle $k$ in scenario $s$.

$$\min\ \sum_{k\in K}\sum_{(i,j)\in A} d_{ij}x_{ijk}\;+\;\sum_{s\in S}p_s\sum_{k\in K}P\,n_{ks}$$

$$
\begin{aligned}
&(1)\text{–}(6),\ (8)\text{–}(10)\ \text{ with } M^D_{ij}=U+d_{ij}\ \text{ (range constraint (7) removed)}\\
&R_s\,(u_k+n_{ks})\ \ge\ \textstyle\sum_{(i,j)\in A} d_{ij}x_{ijk} && \forall k\in K,\ s\in S &(11)\\
&u_k = 1-x_{0,n+1,k} && \forall k\in K &(12)\\
&n_{ks}\in\mathbb{Z}_{\ge 0} && \forall k\in K,\ s\in S &(13)
\end{aligned}
$$

$U$ is an upper bound on any route length ($\sum_i\max_j d_{ij}$), and $u_k$ = 1 if vehicle $k$ is used.

**Justification**
- One plan for the whole year: $x$ (and the schedule) are first-stage and the same in every month. Only the recharges depend on the month.
- Charging is possible anywhere and restores the full range $R_s$. A route of length $L_k$ therefore needs exactly $\lceil L_k/R_s\rceil-1$ recharges: charge at the moment the battery is empty, and every charge adds at most $R_s$. (11) with integer $n_{ks}$ and a minimising objective gives exactly this value. The code checks it by simulating every route in every month.
- $u_k$ instead of 1 in (11) is equivalent for integer $x$ (unused vehicle: $L_k=0$, $n_{ks}=0$) but gives a tighter LP relaxation. Together with the route-length symmetry breaking, this cut the no-time-window run at P = 60 from 436 s to 19 s (HiGHS).
- Both stages are solved jointly (extensive form, deterministic equivalent). With equal probabilities, the expected penalty is $\frac{P}{12}\times$(total number of recharges over the 12 months). The objective trades distance against the total number of recharges, as the hint says.
- Units: the objective is expected cost per day in distance units, so $P$ is "how many km one recharge is worth".

## Part 3b – results (small instance, K = 3)

| P | Time windows | Distance | Vehicles | Recharges / 12 months | E[penalty]/day | Objective | Routes |
|---|---|---|---|---|---|---|---|
| 2 | with | 208.38 | 1 | 15 | 2.50 | **210.88** | 0-3-4-2-10-7-6-8-1-9-5-0 |
| 60 | with | 212.53 | 2 | 6 (2 each in Jan, Feb, Dec) | 30.00 | **242.53** | = Part 1 routes |
| 150 | with | 251.89 | 3 | 0 | 0 | **251.89** | 0-6-8-1-0 (85.81) · 0-3-4-7-0 (84.96) · 0-10-2-9-5-0 (81.11) |
| 2 | without | 158.29 | 1 | 12 | 2.00 | **160.29** | 0-5-9-2-10-4-3-7-1-6-8-0 |
| 60 | without | 193.71 | 2 | 4 (Jan, Feb, Nov, Dec) | 20.00 | **213.71** | 0-4-3-7-1-6-8-0 (112.60) · 0-10-2-9-5-0 (81.11) |
| 150 | without | 234.80 | 3 | 0 | 0 | **234.80** | 0-10-2-9-5-0 · 0-4-3-7-0 · 0-8-6-1-0 |

HiGHS run times: with time windows 0.2–1.7 s, without 11–21 s.
The script also prints, for each month and vehicle, how many recharges there are and on which leg / at which km they happen.

**Why these penalties**: a scan over P showed where the optimal plan changes.
- With time windows: 1 vehicle (P ≲ 4), a 2-vehicle unbalanced plan 81.1/130.0 (P ≈ 4–17), the Part 1 plan (P ≈ 17–79), then 3 vehicles with 0 recharges (P ≳ 79).
- Without time windows: 1 vehicle (P ≲ 53), 2 vehicles (P ≈ 53–123), 3 vehicles (P ≳ 123).
- P = 2, 60, 150 give 1, 2 and 3 vehicles in both cases.

**Discussion points**
- Low P: cheap recharging makes the range irrelevant. The model drives the single shortest feasible tour and accepts recharges every month.
- High P: the model buys "range robustness" with distance. All routes ≤ 89.35 (the worst month), so no recharge ever happens. This is more robust than α = 0.90 (fails in January) and needs fewer vehicles than α = 0.99.
- Compare with Part 1/2 through the *value of the stochastic solution*. Part 1 is exactly the expected-value problem (mean range). Its expected cost under the recourse model is:
  - With time windows: 212.53 + 6P/12 = 213.53 / 242.53 / 287.53 for P = 2/60/150, vs. the optimum 210.88 / 242.53 / 251.89. VSS = 2.65 / 0 / **35.64**.
  - Without time windows: 190.80 + 5P/12 = 191.63 / 215.80 / 253.30, vs. 160.29 / 213.71 / 234.80. VSS = 31.34 / 2.09 / 18.50.
  - Taking the distribution into account clearly pays off.
- Chance constraints only control *how often* the range is exceeded, not *how costly* that is. The recourse model prices exceedances and so picks the right plan for each penalty level.
- Computationally the recourse model is larger (12 scenarios, 36 integer $n_{ks}$) but still easy with time windows. Without time windows it is harder, because the range no longer bounds the route lengths.

## Part 3c – large instance (21 locations, time windows)

**Choosing the fleet**: time windows alone need only 2 vehicles; the deterministic mean range (Part 1 logic) needs 3. A plan with zero recharges in every month (all routes ≤ 89.35) needs **5 vehicles**; this was verified exactly by enumerating all routes ≤ 89.35. So **K = 5**: the same reasoning as K = 3 for the small instance, and it lets the penalty decide freely between "recharge" and "add a vehicle".

| P | Distance | Vehicles used | Recharges / 12 months | Objective | Run time (HiGHS) |
|---|---|---|---|---|---|
| 2 | 331.51 | 3 | 11 | 333.34 | 21 s |
| 60 | 337.86 | 4 | 9 | 382.86 | 245 s |
| 150 | 397.21 | 5 | 0 | 397.21 | 343 s |

Routes:
- P = 2: 0-3-4-2-10-19-15-9-5-0 (127.67) · 0-18-16-20-11-13-14-12-17-0 (108.22) · 0-21-7-6-8-1-0 (95.63)
- P = 60: 0-18-16-4-13-14-12-17-0 (108.47) · 0-20-11-2-10-19-15-9-5-0 (107.80) · 0-3-7-6-8-1-0 (107.45) · 0-21-0 (14.14)
- P = 150: 0-18-16-4-2-12-17-0 (88.37) · 0-20-3-7-1-0 (84.31) · 0-11-13-14-15-0 (83.63) · 0-21-6-8-0 (72.42) · 0-10-19-9-5-0 (68.48)

Remarkable at P = 60: the 4th vehicle serves only node 21 (next to the depot, 14.1 km round trip). The three other routes stay at about 108, so they need a recharge only in Jan, Feb and Dec (3 × 3 = 9). Adding node 21 to one of them would most likely push it above the November/March ranges (112.15/113.55), costing more recharges than the short extra trip.

**Discussion points vs 3b**
- Same qualitative pattern: more vehicles and more distance as P rises, and zero recharges at P = 150.
- Unlike the small instance, low P never leads to a single vehicle: one vehicle cannot meet all 21 time windows.
- Per location, recharges are lower than in 3b at P = 2 (11 vs 15 for twice as many locations), because the work is spread over 3 vehicles anyway.
- Run time grows strongly: about 2× the locations and 5 vehicles means roughly 1,350 binary arc variables vs 216. Report the Gurobi times and, if any run hits 1800 s, the gap.

---

## Part 3d – new recourse: planned depot recharging during idle time

**Idea and justification**
- Part 1 showed long idle periods caused by the time windows (up to ~600 min waiting).
- A realistic recourse for an e-fleet operator: in a month with low range, a vehicle uses such a gap to detour via the depot between two pickups ($i\to0\to j$) and fast-charge to full range (charging time $\tau$).
- Cost: the extra distance $d_{i0}+d_{0j}-d_{ij}$ plus a small handling cost $F$, usually far cheaper than an emergency call-out.
- Unlike 3b, the detour **takes time**. It is only possible if all later time windows can still be met *in that month*, so battery levels and service times become scenario-dependent (second stage).
- Emergency charging (penalty $P$, as in 3b) is kept as a backup. That gives *complete recourse*: every first-stage plan has a feasible second stage.
- With detours switched off, the model reproduces the 3b optimum exactly (checked for all P).

**New parameters**: $\tau=30$ (charging time), $F=5$; sets $A^Z=\{(i,j)\in A: i,j\in C,\ e_i+s_i+d_{i0}+\tau+d_{0j}\le l_j\}$ (arcs where a detour can fit), $A^W=A\setminus\{(0,n+1)\}$.

**Second-stage variables (per $s$)**: $z_{ijs}\in\{0,1\}$ depot detour on arc $(i,j)$; $w_{ijs}\in\{0,1\}$ emergency recharge on arc $(i,j)$; $b_{is}\in[0,R_s]$ remaining range on arrival at $i\in C$; $T^s_i\in[e_i,l_i]$ service start at $i\in C$ in month $s$. Write $y_{ij}=\sum_k x_{ijk}$ (each arc between pickups is used by at most one vehicle, so the second-stage variables need no vehicle index).

$$\min\ \sum_{k}\sum_{(i,j)\in A} d_{ij}x_{ijk}+\sum_{s}p_s\Big[\sum_{(i,j)\in A^Z}\big(F+d_{i0}+d_{0j}-d_{ij}\big)z_{ijs}+\sum_{(i,j)\in A^W}P\,w_{ijs}\Big]$$

subject to (1)–(6), (8)–(10) and, for all $s\in S$:

$$
\begin{aligned}
& z_{ijs}+w_{ijs}\le y_{ij} && (i,j)\in A^W &(14)\\
& b_{js}\le R_s-d_{0j}+R_s w_{0js}+d_{0j}(1-y_{0j}) && j\in C &(15)\\
& b_{js}\le b_{is}-d_{ij}+R_s w_{ijs}+(R_s+d_{ij})(1-y_{ij}+z_{ijs}) && i,j\in C &(16)\\
& b_{is}\ge d_{i0}\,z_{ijs},\qquad b_{js}\le R_s-d_{0j}+d_{0j}(1-z_{ijs}) && (i,j)\in A^Z &(17)\\
& b_{is}+R_s w_{i,n+1,s}\ge d_{i0}\,y_{i,n+1} && i\in C &(18)\\
& T^s_j\ge T^s_i+s_i+d_{ij}-M_{ij}(1-y_{ij}) && i,j\in C &(19)\\
& T^s_j\ge T^s_i+s_i+d_{i0}+\tau+d_{0j}-M^Z_{ij}(1-z_{ijs}) && (i,j)\in A^Z &(20)\\
& T^s_j\ge (e_0+d_{0j})\,y_{0j},\quad T^s_i+s_i+d_{i0}\le l_0+M^E_i(1-y_{i,n+1}) && i,j\in C &(21)\\
& e_i\le T^s_i\le l_i,\ 0\le b_{is}\le R_s,\ z_{ijs},w_{ijs}\in\{0,1\} &&&(22)
\end{aligned}
$$

What the constraints do:
- (15)–(18) track the battery: start full; drive directly; or go to the depot (must be reachable) and restart full there; emergency adds $R_s$, capped at $R_s$; the vehicle must reach the end depot.
- (19)–(21) are the month-specific schedule, including detour and charging time.
- Emergency charging has no time impact, as in 3b.

**Results (small instance, time windows, K = 3, τ = 30, F = 5)**

| P | Plan | Depot recharges / 12 m | Emergency / 12 m | E[recourse]/day | Objective | 3b objective |
|---|---|---|---|---|---|---|
| 2 | same as 3b (1 vehicle) | 0 | 15 | 2.50 | 210.88 | 210.88 |
| 60 | same routes as 3b/Part 1 | 4 (Jan, Feb ×2, Dec) | 2 | 25.55 | **238.08** | 242.53 |
| 150 | **new plan**: 0-3-4-7-1-0 · 0-10-2-9-5-0 · 0-8-6-0 (246.41) | 1 (Jan, 7→0→1) | 0 | 3.84 | **250.26** | 251.89 |

Sensitivity at P = 60: τ = 0 or 90 gives the same result. With τ = 180 only 2 detours still fit the time windows (4 emergencies remain), and the cost becomes 240.62.

**Discussion points**
- With cheap emergencies (P = 2) the new option is never used: a detour costs at least F plus extra km.
- At P = 60 the routes stay the same, but 4 of 6 emergencies become depot recharges in waiting time. The remaining 2 are on vehicle 2 late in its route (6→8, 8→1), where there is no slack: node 1 starts at 969 with window [819, 979].
- At P = 150 the first-stage plan changes: the model accepts a route of 91.80 (> January range) because January can be fixed with one depot detour. It saves 5.48 km every day of the year. This is the α = 0.90 plan from Part 2, but now with a priced and feasible recourse.
- Interaction with time windows: depot recharging depends on idle time. Longer charging times remove options (the τ sensitivity), so the timetable is now a real part of the recourse decision.
- Model size: ~2,000 variables (1,644 integer) vs ~330 in 3b, and HiGHS times of 2–43 s. Fine for the small instance; for the large one expect long runs.

---

## Report checklist
- Group number and names in every `.py` file and in the report; contribution statement.
- The assignment asks for a **generative-AI statement per question**: name the tool, the prompts, and how the answers were used. This chat counts, so keep the prompts.
- Re-run all scripts with Gurobi and replace the run times. If a run hits the 1800 s limit, report the gap (printed automatically).
- Figures: `output/*.png` (route maps for every experiment).
