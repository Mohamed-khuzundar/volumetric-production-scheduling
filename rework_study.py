# ==========================================================================
# REWORK STUDY  -  standalone, single-file simulation sweep
#
# WHAT IT DOES
#   Reproduces the full dataset behind the rework results in the dissertation.
#   A rework event is a module type + a day: the first module of that type
#   inspected on or after that day fails inspection and must repeat its MEP and
#   inspection stages, adding work at the bottleneck. The line is simulated at
#   five frequency levels (1..5 reworks per run) and three timing bands
#   (early / middle / late), 200 runs per level = 1000 runs. For each run three
#   scheduling strategies are scored on floor tardiness and makespan:
#       passive  - redo appended after the whole planned line
#       model1   - rule-based reactive (redo re-inserted by earliest due date)
#       model2   - genetic-algorithm reactive
#
#   The design is 67 early / 67 middle / 66 late per level, produced in two
#   parts, exactly as in the study, both handled automatically here:
#     * a main sweep (seed 7): 67 early + 67 middle + 20 late per level;
#     * a late-band top-up (seed 20260714): 46 further late runs per level.
#   The two are concatenated per level to give the 200-run design.
#
# HOW TO RUN
#   Needs only Python 3 (standard library only - nothing to install).
#       python rework_study.py
#   Runs the full 1000-run sweep (about 75-90 min) and writes
#   rework_results.csv next to this file. To run in shorter chunks:
#       python rework_study.py 0 200
#       python rework_study.py 200 400   ... up to 1000
#
# OUTPUT  (rework_results.csv, one row per run; tardiness & makespan in days)
#   run, k_reworks, realised_failures, timing_band, first_start_day, rework_types,
#   passive_tardiness, model1_tardiness, model2_tardiness,
#   passive_makespan, model1_makespan, model2_makespan
#
#   realised_failures may be below k when late reworks leave too little of a
#   type still in production to fail; this is recorded, not hidden.
#
# Running this script regenerates rework_results.csv exactly; that file is the
# dataset analysed in the dissertation's rework study.
# ==========================================================================
import random, json, sys, time, os
MIN=480
PROC={"A":[67,178,155,44],"B":[225,599,524,150],"C":[353,941,823,235]}
NF=6; FLOORP=["A"]*8+["B"]*2+["C"]*7
def baseline_items():
    items=[]; i=0
    for f in range(1,NF+1):
        for t in FLOORP: items.append(dict(f=f,t=t,kind="N",rid=None,oid=i)); i+=1
    return items

def run_rework(seq):
    n=len(seq); free=[0.0]*4; C=[{} for _ in range(n)]; F_insp={}; mep={}; insp={}
    for j in range(n):
        it=seq[j]; route=[1,3] if it["kind"]=="R" else [0,1,2,3]; prev=None
        for s in route:
            arr=prev if prev is not None else 0.0
            if s==1 and it["kind"]=="R": arr=max(arr,F_insp.get(it["rid"],0.0))
            start=arr if arr>free[s] else free[s]
            c=start+PROC[it["t"]][s]; free[s]=c; C[j][s]=c; prev=c
            if s==1: mep[j]=start
        insp[j]=C[j][3]
        if it["kind"]=="F": F_insp[it["rid"]]=C[j][3]
    ff={}
    for j in range(n):
        it=seq[j]
        if it["kind"] in ("N","R"): ff[it["f"]]=max(ff.get(it["f"],0.0),C[j][3])
    return [ff[f]/MIN for f in range(1,NF+1)], mep, insp

base=baseline_items()
TARGET,_,_=run_rework(base)
def total_tardiness(fl): return sum(max(0.0,f-t) for f,t in zip(fl,TARGET))
def makespan(fl): return max(fl)

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

# ---- GA: population 100, generations 250, elitism 8, tournament 5 ----
def ga_order(locked, free, gens=250, pop=100, seed=1, extra_seeds=None):
    pool=free; m=len(pool)
    if m<=1: return repair(pool[:])
    def fit(o): return total_tardiness(run_rework(locked+repair([pool[i] for i in o]))[0])
    byfloor=sorted(range(m),key=lambda i:(pool[i]["f"],{"F":0,"N":1,"R":2}[pool[i]["kind"]]))
    def seeds():
        S=[list(range(m)),byfloor]
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
        nP=[o for _,o in sc[:8]]                      # elitism = 8
        while len(nP)<pop:
            t1=min(random.sample(sc,5),key=lambda x:x[0])[1]   # tournament 5
            t2=min(random.sample(sc,5),key=lambda x:x[0])[1]
            c=ox(t1,t2)
            if random.random()<0.6: c=mut(c)                    # per-individual mutation 0.6
            nP.append(c)
        P=nP
    return repair([pool[i] for i in best])


# ---- Option B: identify the failed modules ONCE on the baseline schedule ----

# ---- Option X: which module types are still available to fail at a given day ----
_, _, _BASE_INSP = run_rework(base)
def types_available(day):
    """Types with at least one baseline module inspected on/after 'day'."""
    t = day * MIN
    avail = set()
    for j in range(len(base)):
        if _BASE_INSP[j] >= t - 1e-9:
            avail.add(base[j]["t"])
    return avail

def failures_on_baseline(rules):
    """Return list of dicts {oid,f,t,rid,T} for the modules that fail, decided on the
       baseline order so every model faces the identical disruption."""
    order=[dict(base[j]) for j in range(len(base))]
    _,_,insp=run_rework(order)
    used=set(); remaining=list(range(len(rules))); rid=0; out=[]
    while remaining:
        best=None
        for ri in remaining:
            typ,day=rules[ri]
            cands=[(insp[p],p) for p in range(len(order))
                   if order[p]["oid"] not in used and order[p]["t"]==typ and insp[p]>=day*MIN-1e-9]
            if cands:
                Tf,p=min(cands)
                if best is None or Tf<best[0]: best=(Tf,ri,p)
        if best is None: break
        Tf,ri,p=best; used.add(order[p]["oid"])
        out.append(dict(oid=order[p]["oid"],f=order[p]["f"],t=order[p]["t"],rid=rid,T=Tf))
        remaining.remove(ri); rid+=1
    return out

def reactive_tp(fails, policy):
    order=[dict(base[j]) for j in range(len(base))]
    fails=sorted(fails,key=lambda d:d["T"])
    for d in fails:
        # locate the fixed failed module by oid, mark it failed now
        p=next(q for q in range(len(order)) if order[q]["oid"]==d["oid"])
        Tf_row=run_rework(order)[2][p]           # its inspection time in THIS schedule
        order[p]=dict(order[p],kind="F",rid=d["rid"])
        _,mep,_=run_rework(order)
        locked=[]; free=[]
        for q in range(len(order)):
            (locked if mep[q]<Tf_row-1e-9 else free).append(order[q])
        free=repair(free+[dict(f=d["f"],t=d["t"],kind="R",rid=d["rid"],oid=None)])
        if policy=="ga": free=ga_order(locked,free)
        else: free=repair(sorted(free,key=lambda it:(it["f"],{"F":0,"N":1,"R":2}[it["kind"]])))
        order=locked+free
    return order, fails

def passive_tp(fails):
    order=[dict(base[j]) for j in range(len(base))]
    rid_of={d["oid"]:d["rid"] for d in fails}
    items=[]
    for j in range(len(base)):
        if base[j]["oid"] in rid_of:
            items.append(dict(f=base[j]["f"],t=base[j]["t"],kind="F",rid=rid_of[base[j]["oid"]],oid=base[j]["oid"]))
            items.append(dict(f=base[j]["f"],t=base[j]["t"],kind="R",rid=rid_of[base[j]["oid"]],oid=None))
        else:
            items.append(dict(base[j]))
    return repair(items), fails

def evaluate(rules):
    rules=sorted(rules,key=lambda r:r[1])
    fails=failures_on_baseline(rules)
    react_seq, _ = reactive_tp(fails,"ga")
    m1_seq,    _ = reactive_tp(fails,"edd")
    pass_seq,  _ = passive_tp(fails)
    flP=run_rework(pass_seq)[0]; fl1=run_rework(m1_seq)[0]; flR=run_rework(react_seq)[0]
    return dict(passive=total_tardiness(flP), m1=total_tardiness(fl1), react=total_tardiness(flR),
                ms_passive=round(makespan(flP),3), ms_m1=round(makespan(fl1),3), ms_react=round(makespan(flR),3),
                nfail=len(fails))

# ---- stratified scenario: first rework in [lo,hi]; the rest at/after it ----
def scenario_stratified(rng, k, lo, hi):
    first_day = rng.uniform(lo, hi)
    days = sorted([round(first_day,1)] + [round(rng.uniform(first_day,117),1) for _ in range(k-1)])
    rules=[]
    for d in days:
        avail = sorted(types_available(d))          # only types still in production at day d
        if not avail: avail = ["C"]                  # safety (very late): C is last to finish
        rules.append((rng.choice(avail), d))
    return rules

def describe(rules): return "; ".join(f"{t} d{d:.0f}" for t,d in sorted(rules,key=lambda r:r[1]))

LEVELS=[1,2,3,4,5]
BANDS=[("early",100,(1,30)),("middle",80,(31,90)),("late",20,(91,117))]
SEED=7


# ==========================================================================
# STUDY DRIVER  -  builds the 67/67/66 design from two seeds and writes CSV
# ==========================================================================
MAIN_SEED  = SEED
TOPUP_SEED = 20260714
NEW_LATE_PER_LEVEL = 46


def build_design():
    # main stream (seed 7): 100 early / 80 middle / 20 late per level
    rng_main = random.Random(MAIN_SEED)
    main = {}
    for k in LEVELS:
        for bname, count, (lo, hi) in BANDS:
            lst = main.setdefault((k, bname), [])
            for _ in range(count):
                lst.append(scenario_stratified(rng_main, k, lo, hi))

    # top-up stream (seed 20260714): 46 new late per level.
    # (The rework top-up uses no duplicate rejection, matching the study.)
    rng_top = random.Random(TOPUP_SEED)
    topup = {k: [] for k in LEVELS}
    for k in LEVELS:
        for _ in range(NEW_LATE_PER_LEVEL):
            topup[k].append(scenario_stratified(rng_top, k, 91, 117))

    design = []
    for k in LEVELS:
        for rules in main[(k, "early")][:67]:
            design.append((k, "early", rules))
        for rules in main[(k, "middle")][:67]:
            design.append((k, "middle", rules))
        for rules in main[(k, "late")][:20]:
            design.append((k, "late", rules))
        for rules in topup[k][:NEW_LATE_PER_LEVEL]:
            design.append((k, "late", rules))
    return design


if __name__ == "__main__":
    import csv
    design = build_design()
    N = len(design)

    if len(sys.argv) >= 3:
        start_i, end_i = int(sys.argv[1]), int(sys.argv[2])
    else:
        start_i, end_i = 0, N

    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "rework_results.csv")
    header_row = ["run", "k_reworks", "realised_failures", "timing_band",
                  "first_start_day", "rework_types",
                  "passive_tardiness", "model1_tardiness", "model2_tardiness",
                  "passive_makespan", "model1_makespan", "model2_makespan"]
    mode = "w" if start_i == 0 else "a"

    t0 = time.time()
    with open(path, mode, newline="") as f:
        w = csv.writer(f)
        if start_i == 0:
            w.writerow(header_row)
        for i in range(start_i, end_i):
            k, band, rules = design[i]
            r = evaluate(rules)
            first_start = min(d for _, d in rules)
            types_str = "|".join(t for t, _ in sorted(rules, key=lambda x: x[1]))
            w.writerow([i + 1, k, r["nfail"], band, round(first_start, 1), types_str,
                        round(r["passive"], 3), round(r["m1"], 3), round(r["react"], 3),
                        round(r["ms_passive"], 3), round(r["ms_m1"], 3), round(r["ms_react"], 3)])
            f.flush()
    print("runs %d-%d done in %.0fs  (total: %d, output: rework_results.csv)"
          % (start_i, end_i, time.time() - t0, N))
