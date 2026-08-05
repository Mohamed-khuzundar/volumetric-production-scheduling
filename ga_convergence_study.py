# ==========================================================================
# GA CONVERGENCE STUDY  -  standalone, single-file script
#
# WHAT IT DOES
#   Runs the genetic algorithm (population 100) on the extreme mixed-disruption
#   scenario (3 MEP stoppages + 3 material shortages + 3 reworks) for five random
#   seeds, recording the best floor tardiness at every generation. It reports the
#   generation at which the search converges and saves the per-generation curves.
#
# HOW TO RUN
#   1.  Install the one dependency:   pip install matplotlib
#   2.  Run the file:                 python GA_convergence_study.py
#       (takes ~1 minute; no arguments needed)
#
# OUTPUTS (written into the current folder)
#   convergence_curves.csv    best-so-far tardiness per generation, per seed
#   fig_ga_convergence.png    the convergence figure
#
# The simulation engine is embedded in this file, so nothing else is required.
# ==========================================================================
import random, json, sys, time, os
from collections import defaultdict, deque
MIN=480
PROC={"A":[67,178,155,44],"B":[225,599,524,150],"C":[353,941,823,235]}
NF=6; FLOORP=["A"]*8+["B"]*2+["C"]*7
def baseline_items():
    items=[]; i=0
    for f in range(1,NF+1):
        for t in FLOORP: items.append(dict(f=f,t=t,kind="N",rid=None,oid=i)); i+=1
    return items

def run_combined(seq, windows):
    stoppages=sorted([w for w in windows if w["kind"]=="stoppage"],key=lambda w:w["start"])
    delays=[w for w in windows if w["kind"]=="delay"]
    n=len(seq); free=[0.0]*4; C=[{} for _ in range(n)]; F_insp={}; mep={}; insp={}
    for j in range(n):
        it=seq[j]; route=[1,3] if it["kind"]=="R" else [0,1,2,3]; prev=None
        for s in route:
            arr=prev if prev is not None else 0.0
            if s==1 and it["kind"]=="R": arr=max(arr,F_insp.get(it["rid"],0.0))
            start=arr if arr>free[s] else free[s]
            if s==1:
                moved=True
                while moved:
                    moved=False
                    for w in stoppages:
                        if w["start"]<=start<w["end"]: start=w["end"]; moved=True; break
                    if moved: continue
                    for w in delays:
                        if w["type"]==it["t"] and w["start"]<=start<w["end"]: start=w["end"]; moved=True; break
                mep[j]=start
                rem=PROC[it["t"]][1]; tcur=start; comp=None
                for w in stoppages:
                    if w["end"]<=tcur: continue
                    if tcur<w["start"]:
                        avail=w["start"]-tcur
                        if avail>=rem: comp=tcur+rem; break
                        rem-=avail; tcur=w["end"]
                    else: tcur=w["end"]
                if comp is None: comp=tcur+rem
                c=comp
            else: c=start+PROC[it["t"]][s]
            free[s]=c; C[j][s]=c; prev=c
        insp[j]=C[j][3]
        if it["kind"]=="F": F_insp[it["rid"]]=C[j][3]
    ff={}
    for j in range(n):
        it=seq[j]
        if it["kind"] in ("N","R"): ff[it["f"]]=max(ff.get(it["f"],0.0),C[j][3])
    return [ff[f]/MIN for f in range(1,NF+1)], mep, insp

base=baseline_items()
TARGET,_,base_insp=run_combined(base,[])
def total_tardiness(fl): return sum(max(0.0,f-t) for f,t in zip(fl,TARGET))

def repair(items):
    out=items[:]
    while True:
        moved=False
        for i,it in enumerate(out):
            if it["kind"]=="R":
                fp=next((k for k,x in enumerate(out) if x["kind"]=="F" and x["rid"]==it["rid"]),None)
                if fp is not None and i<fp:
                    r=out.pop(i); fp=next(k for k,x in enumerate(out) if x["kind"]=="F" and x["rid"]==it["rid"])
                    out.insert(fp+1,r); moved=True; break
        if not moved: return out

def edd_items(free): return repair(sorted(free, key=lambda it:(it["f"], {"F":0,"N":1,"R":2}[it["kind"]])))

def ga_order(locked, free, windows, gens=55, pop=40, seed=1, extra_seeds=None):
    pool=free; m=len(pool)
    if m<=1: return pool[:]
    def fit(o): return total_tardiness(run_combined(locked+repair([pool[i] for i in o]), windows)[0])
    byfloor=sorted(range(m),key=lambda i:(pool[i]["f"],{"F":0,"N":1,"R":2}[pool[i]["kind"]]))
    def seeds():
        S=[list(range(m)), byfloor, byfloor[::-1]]
        if extra_seeds:
            for es in extra_seeds:
                if sorted(es)==list(range(m)): S.append(list(es))
        for _ in range(5): r=list(range(m)); random.shuffle(r); S.append(r)
        return S
    def ox(p1,p2):
        a,b=sorted(random.sample(range(m),2)); seg=set(p1[a:b+1])
        c=[None]*m; c[a:b+1]=p1[a:b+1]; fill=[x for x in p2 if x not in seg]; k=0
        for i in range(m):
            if c[i] is None: c[i]=fill[k]; k+=1
        return c
    def mut(o):
        o=o[:]
        if random.random()<0.5: i=random.randrange(m); x=o.pop(i); o.insert(random.randrange(m),x)
        else: i,j=random.sample(range(m),2); o[i],o[j]=o[j],o[i]
        return o
    random.seed(seed); P=seeds()
    while len(P)<pop: P.append(random.sample(range(m),m))
    best=None; bt=float("inf")
    for g in range(gens):
        sc=sorted([(fit(o),o) for o in P],key=lambda x:x[0])
        if sc[0][0]<bt: bt,best=sc[0][0],sc[0][1][:]
        nP=[o for _,o in sc[:8]]
        while len(nP)<pop:
            t1=min(random.sample(sc,5),key=lambda x:x[0])[1]
            t2=min(random.sample(sc,5),key=lambda x:x[0])[1]
            c=ox(t1,t2)
            if random.random()<0.6: c=mut(c)
            nP.append(c)
        P=nP
    return repair([pool[i] for i in best])

def determine_reworks(rework_rules):
    obi=sorted(range(len(base)), key=lambda j: base_insp[j])
    rms=[]; used=set(); rid=0
    for (typ,day) in rework_rules:
        for j in obi:
            if j not in used and base[j]["t"]==typ and base_insp[j]>=day*MIN-1e-9:
                used.add(j); rms.append(dict(oid=j,t=typ,f=base[j]["f"],rid=rid,day=day,T=base_insp[j])); rid+=1; break
    return rms

def pos_by_oid(order): return {order[p]["oid"]:p for p in range(len(order)) if order[p].get("oid") is not None}

def reactive(windows_all, rework_modules, policy):
    order=[dict(base[j]) for j in range(len(base))]
    rid_of={rm["oid"]:rm["rid"] for rm in rework_modules}
    f_of={rm["oid"]:rm["f"] for rm in rework_modules}; t_of={rm["oid"]:rm["t"] for rm in rework_modules}
    known=[]; unprocessed=sorted(windows_all,key=lambda w:w["start"]); unhandled=set(rm["oid"] for rm in rework_modules)
    while unprocessed or unhandled:
        _,_,insp=run_combined(order, known)
        pos=pos_by_oid(order)
        next_w=unprocessed[0] if unprocessed else None
        next_r_oid=None; next_r_t=float("inf")
        for oid in unhandled:
            t=insp[pos[oid]]
            if t<next_r_t: next_r_t=t; next_r_oid=oid
        use_window = next_w is not None and (next_r_oid is None or next_w["start"]<=next_r_t)
        if use_window:
            et=next_w["start"]; known=known+[next_w]; unprocessed=unprocessed[1:]; addR=None
        else:
            et=next_r_t; oid=next_r_oid; p=pos[oid]; rid=rid_of[oid]
            order[p]=dict(order[p],kind="F",rid=rid)
            addR=dict(f=f_of[oid],t=t_of[oid],kind="R",rid=rid,oid=None); unhandled.discard(oid)
        _,mep,_=run_combined(order, known)
        locked=[]; free=[]
        for q in range(len(order)):
            (locked if mep[q]<et-1e-9 else free).append(order[q])
        if addR is not None: free=free+[addR]
        free=repair(free)
        free = ga_order(locked,free,known) if policy=="ga" else edd_items(free)   # honest: only KNOWN windows
        order=locked+free
    return order

def passive(windows_all, rework_modules):
    rid_of={rm["oid"]:rm["rid"] for rm in rework_modules}
    items=[]
    for j in range(len(base)):
        if j in rid_of: items.append(dict(f=base[j]["f"],t=base[j]["t"],kind="F",rid=rid_of[j],oid=j))
        else: items.append(dict(base[j]))
    R=[dict(f=rm["f"],t=rm["t"],kind="R",rid=rm["rid"],oid=None) for rm in rework_modules]
    return repair(items+R)

def predictive(windows_all, rework_modules, seed_seq):
    order=[dict(base[j]) for j in range(len(base))]
    rid_of={rm["oid"]:rm["rid"] for rm in rework_modules}
    for p in range(len(order)):
        if order[p]["oid"] in rid_of: order[p]=dict(order[p],kind="F",rid=rid_of[order[p]["oid"]])
    ev=[w["start"] for w in windows_all]+[rm["T"] for rm in rework_modules]
    T1=min(ev) if ev else 0.0
    _,mep,_=run_combined(order, [])
    locked=[]; free=[]
    for q in range(len(order)):
        (locked if mep[q]<T1-1e-9 else free).append(order[q])
    R=[dict(f=rm["f"],t=rm["t"],kind="R",rid=rm["rid"],oid=None) for rm in rework_modules]
    free=repair(free+R)
    extra=None
    if seed_seq is not None:
        def key(it): return ("R",it["rid"]) if it["kind"]=="R" else (it["kind"],it["oid"])
        fk={key(free[i]):i for i in range(len(free))}
        es=[fk[key(it)] for it in seed_seq if key(it) in fk]
        if sorted(es)==list(range(len(free))): extra=[es]
    free=ga_order(locked,free,windows_all,extra_seeds=extra)
    return locked+free

def evaluate(windows, rms):
    react_seq=reactive(windows,rms,"ga")
    m1_seq=reactive(windows,rms,"edd")
    pass_seq=passive(windows,rms)
    pred_seq=predictive(windows,rms,react_seq)
    return dict(passive=total_tardiness(run_combined(pass_seq,windows)[0]),
                m1=total_tardiness(run_combined(m1_seq,windows)[0]),
                react=total_tardiness(run_combined(react_seq,windows)[0]),
                pred=total_tardiness(run_combined(pred_seq,windows)[0]))

def random_scenario(rng):
    for _ in range(40):
        k=rng.randint(1,5); windows=[]; rules=[]
        for _ in range(k):
            kind=rng.choice(["delay","stoppage","rework"])
            if kind=="delay":
                st=round(rng.uniform(3,95),1); du=round(rng.uniform(2,5),1)
                windows.append(dict(kind="delay",type=rng.choice("ABC"),start=st*MIN,end=(st+du)*MIN))
            elif kind=="stoppage":
                st=round(rng.uniform(3,95),1); du=round(rng.uniform(2,5),1)
                windows.append(dict(kind="stoppage",start=st*MIN,end=(st+du)*MIN))
            else:
                rules.append((rng.choice("ABC"), round(rng.uniform(3,95),1)))
        rms=determine_reworks(rules)
        if total_tardiness(run_combined(passive(windows,rms),windows)[0])>1.0:
            return windows, rules, rms
    return windows, rules, rms

def describe(windows, rules):
    parts=[]
    for w in sorted(windows,key=lambda w:w["start"]):
        d=(w["end"]-w["start"])/MIN
        parts.append((f"D{w['type']}" if w["kind"]=="delay" else "S")+f"d{w['start']/MIN:.0f}+{d:.0f}")
    for (ty,dy) in rules: parts.append(f"RW{ty}d{dy:.0f}")
    return "; ".join(parts)


# ===================== STUDY CODE =====================
"""
GA CONVERGENCE STUDY (generations only, population fixed at 100).
Runs the validated GA on the HARDEST rescheduling problem:
 3 MEP stoppages + 3 material delays + 3 Type-C reworks, timed EARLY,
 with the entire set of 102 modules (+5 reworks) free to reorder
 (locked = [] -> largest possible free pool = hardest search).
Records best-so-far tardiness Z at every generation, over several seeds.
Reuses the uploaded, validated engine in ga_sweep_combined.py.
"""
import random, time, json
import sys as _sys; g = _sys.modules[__name__]  # let g.* refer to the embedded engine above
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

MIN = g.MIN

def build_hardest(scenario_seed=42):
    # Extreme severity (aligned with the mixed-study design): 3 stoppages
    # + 3 delays + 3 reworks = 9 disruptions in one run.
    # Each duration is drawn ONCE from Triangular(0.25, 0.5, 3.0) days (mean 1.25)
    # using a fixed scenario seed, so the scenario is held constant across GA seeds.
    srng = random.Random(scenario_seed)
    windows = []
    dur = []
    for st in [4, 10, 16]:
        du = srng.triangular(0.25, 3.0, 0.5)   # (low, high, mode)
        dur.append(round(du, 2))
        windows.append(dict(kind="stoppage", start=st*MIN, end=(st+du)*MIN))
    for st, ty in zip([7, 13, 19], ["C", "C", "B"]):
        du = srng.triangular(0.25, 3.0, 0.5)
        dur.append(round(du, 2))
        windows.append(dict(kind="delay", type=ty, start=st*MIN, end=(st+du)*MIN))
    print("sampled durations (days):", dur, " mean=%.2f" % (sum(dur)/len(dur)))
    rules = [("C", d) for d in [5, 9, 13]]
    rms = g.determine_reworks(rules)
    order = [dict(g.base[j]) for j in range(len(g.base))]
    rid_of = {rm["oid"]: rm["rid"] for rm in rms}
    for p in range(len(order)):
        if order[p]["oid"] in rid_of:
            order[p] = dict(order[p], kind="F", rid=rid_of[order[p]["oid"]])
    R = [dict(f=rm["f"], t=rm["t"], kind="R", rid=rm["rid"], oid=None) for rm in rms]
    free = g.repair(order + R)
    locked = []
    return locked, free, windows

def ga_trace(locked, free, windows, gens, pop, seed):
    """Same operators as the validated GA, but records best-Z per generation."""
    pool = free; m = len(pool)
    def fit(o):
        return g.total_tardiness(g.run_combined(locked + g.repair([pool[i] for i in o]), windows)[0])
    byfloor = sorted(range(m), key=lambda i: (pool[i]["f"], {"F": 0, "N": 1, "R": 2}[pool[i]["kind"]]))
    def make_seeds():
        S = [list(range(m)), byfloor, byfloor[::-1]]
        for _ in range(5):
            r = list(range(m)); random.shuffle(r); S.append(r)
        return S
    def ox(p1, p2):
        a, b = sorted(random.sample(range(m), 2)); seg = set(p1[a:b+1])
        c = [None]*m; c[a:b+1] = p1[a:b+1]; fill = [x for x in p2 if x not in seg]; k = 0
        for i in range(m):
            if c[i] is None: c[i] = fill[k]; k += 1
        return c
    def mut(o):
        o = o[:]
        if random.random() < 0.5:
            i = random.randrange(m); x = o.pop(i); o.insert(random.randrange(m), x)
        else:
            i, j = random.sample(range(m), 2); o[i], o[j] = o[j], o[i]
        return o
    random.seed(seed); P = make_seeds()
    while len(P) < pop: P.append(random.sample(range(m), m))
    bt = float("inf"); curve = []
    for gen in range(gens):
        sc = sorted([(fit(o), o) for o in P], key=lambda x: x[0])
        if sc[0][0] < bt: bt = sc[0][0]
        curve.append(bt)
        nP = [o for _, o in sc[:8]]
        while len(nP) < pop:
            t1 = min(random.sample(sc, 5), key=lambda x: x[0])[1]
            t2 = min(random.sample(sc, 5), key=lambda x: x[0])[1]
            c = ox(t1, t2)
            if random.random() < 0.6: c = mut(c)
            nP.append(c)
        P = nP
    return curve

GENS = 500; POP = 100; SEEDS = [1, 2, 3, 4, 5]
locked, free, windows = build_hardest()
print("free pool size:", len(free))
curves = {}; t0 = time.time()
for s in SEEDS:
    c = ga_trace(locked, free, windows, GENS, POP, s)
    curves[s] = c
    print("seed %d: final Z=%.2f  (%.0fs elapsed)" % (s, c[-1], time.time()-t0))

# Save the per-generation curves as CSV (one row per generation, one column per seed).
# This is the results file that accompanies the script in the repository.
import csv
with open("convergence_curves.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["generation"] + [f"seed_{s}" for s in SEEDS])
    for gi in range(GENS):
        w.writerow([gi + 1] + [round(curves[s][gi], 4) for s in SEEDS])
print("saved convergence_curves.csv")

# ---- convergence analysis ----
def gen_within(curve, frac):
    final = curve[-1]
    thresh = final * (1 + frac)
    for gi, v in enumerate(curve, start=1):
        if v <= thresh: return gi
    return len(curve)

print("\n--- convergence points (generation to reach within X% of each seed's final) ---")
for pct, frac in [("5%", 0.05), ("2%", 0.02), ("1%", 0.01)]:
    gens_needed = [gen_within(curves[s], frac) for s in SEEDS]
    print("within %s: per-seed %s  -> worst-case gen %d" % (pct, gens_needed, max(gens_needed)))
finals = [curves[s][-1] for s in SEEDS]
print("\nfinal Z across seeds: %s  (mean %.2f, spread %.2f)" %
      ([round(x,2) for x in finals], sum(finals)/len(finals), max(finals)-min(finals)))

# ---- plot (clean: adopted line at 250, no zoom inset) ----
import matplotlib.pyplot as plt
colors = ["#1D9E75", "#534AB7", "#D85A30", "#BA7517", "#3A7CA5"]
markers = ['o', 's', '^', 'D', 'v']
fig, ax = plt.subplots(figsize=(9.2, 5.0))
for s, col, mk in zip(SEEDS, colors, markers):
    ax.plot(range(1, GENS+1), curves[s], lw=1.5, color=col, label=f"seed {s}",
            marker=mk, markevery=45, markersize=5,
            markerfacecolor='white', markeredgecolor=col, markeredgewidth=1.2)

# convergence marker (worst-case generation to reach within 1% of the final)
conv_gen = max(gen_within(curves[s], 0.01) for s in SEEDS)

# adaptive y-range from the actual data, so the curves always fill the frame
all_vals = [v for s in SEEDS for v in curves[s]]
z_lo = min(all_vals); z_hi = max(all_vals)
pad = max(1.0, (z_hi - z_lo) * 0.08)
y_bottom = z_lo - pad
y_top = z_hi + pad
y_conv = z_lo + (y_top - z_lo) * 0.45     # arrow head near the converged floor
y_conv_txt = z_lo + (y_top - z_lo) * 0.72
y_adopt_txt = z_lo + (y_top - z_lo) * 0.30

ax.axvline(conv_gen, color="#444", ls="--", lw=1.6)
ax.annotate(f"convergence reached\nby ~{conv_gen} generations", xy=(conv_gen, y_conv),
            xytext=(conv_gen+60, y_conv_txt),
            fontsize=10.5, color="#333", ha="left",
            arrowprops=dict(arrowstyle="->", color="#555", lw=1.2))

# adopted number of generations = 250 (safety margin beyond convergence)
ADOPTED = 250
ax.axvline(ADOPTED, color="#999", ls=":", lw=1.8)
ax.annotate(f"{ADOPTED} generations\nadopted (safety margin)", xy=(ADOPTED, y_adopt_txt),
            xytext=(ADOPTED+8, y_adopt_txt),
            fontsize=10.5, color="#666", ha="left")

ax.set_xlabel("Number of generations", fontsize=12)
ax.set_ylabel("Best floor tardiness Z so far (days)", fontsize=12)
ax.set_title("Genetic Algorithm Convergence at Population = 100\n"
             "Extreme scenario: 3 MEP stoppages + 3 material shortages + 3 reworks, over five random seeds",
             fontsize=11)
ax.set_xlim(0, GENS); ax.set_ylim(y_bottom, y_top)
ax.grid(alpha=0.25)
ax.legend(title="GA runs", fontsize=9.5, title_fontsize=10, frameon=False, loc="upper right")
plt.tight_layout()
plt.savefig("fig_ga_convergence.png", dpi=185, bbox_inches="tight", facecolor="white")
print("\nsaved fig_ga_convergence.png  (adopted line at 250, no zoom inset)")
