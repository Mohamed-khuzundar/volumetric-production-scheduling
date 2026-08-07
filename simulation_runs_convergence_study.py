# ==========================================================================
# STANDALONE FILE - no external imports needed. The validated simulation
# engine (run_combined, reactive, ga_order, etc.) is embedded directly below,
# followed by the study code. Just run this file on its own.
# ==========================================================================

"""
============================================================================
 COMBINED SWEEP  -  Monte Carlo over random MIXED-disruption scenarios
 Each scenario is 1-5 disruptions drawn from {material delay, MEP stoppage,
 rework} in any mix. For each it scores Passive, Model 1 (EDD), GA Reactive
 (honest: re-plans at each event knowing only what has happened) and GA
 Predictive (knows the whole scenario up front), and records the gaps.
 Uses the validated combined recursion. Writes sweep_combined_results.jsonl.
 Run with the play button, or "py ga_sweep_combined.py 0 40" for a chunk.
============================================================================
"""
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
RUN-COUNT CONVERGENCE STUDY.
Fix severity = EXTREME (3 stoppages + 3 material shortages + 3 reworks).
Each run draws fresh random timings, triangular durations, delay types and
rework modules. Run Model 2 (GA reactive) on each and record its tardiness Z.
The running mean of Z is then plotted against the number of runs to see where
it stabilises. Lighter GA settings are used (the running-mean SHAPE depends on
scenario variance, not GA quality) to keep the study tractable.
Chunked: python3 runcount_study.py START END
"""
import random, json, sys, time, os
import sys as _sys; g = _sys.modules[__name__]  # let g.* refer to the embedded engine above
MIN = g.MIN

# lighter GA for tractability (shape of running mean is unaffected)
_orig = g.ga_order
def light_ga(locked, free, windows, gens=25, pop=30, seed=1, extra_seeds=None):
    return _orig(locked, free, windows, gens=gens, pop=pop, seed=seed, extra_seeds=extra_seeds)
g.ga_order = light_ga

def extreme_scenario(rng):
    windows = []
    for _ in range(3):
        st = round(rng.uniform(3, 105), 1); du = round(rng.triangular(0.25, 3.0, 0.5), 2)
        windows.append(dict(kind="stoppage", start=st*MIN, end=(st+du)*MIN))
    for _ in range(3):
        st = round(rng.uniform(3, 105), 1); du = round(rng.triangular(0.25, 3.0, 0.5), 2)
        windows.append(dict(kind="delay", type=rng.choice("ABC"), start=st*MIN, end=(st+du)*MIN))
    rules = [(rng.choice("ABC"), round(rng.uniform(3, 105), 1)) for _ in range(3)]
    rms = g.determine_reworks(rules)
    return windows, rms

TOTAL, MASTER = 500, 20
path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "runcount_results.jsonl")

def make_plot():
    """Deduplicate results (in case chunks overlapped) and plot the running mean."""
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    seen = {}
    for line in open(path):
        r = json.loads(line)
        if r["i"] not in seen: seen[r["i"]] = r["Z"]
    idx = sorted(seen); Z = [seen[i] for i in idx]; n = len(Z)
    run = [sum(Z[:i+1])/(i+1) for i in range(n)]; final = run[-1]
    def settle(band):
        for i in range(n):
            if all(abs(run[j]-final) <= band for j in range(i, n)): return i+1
        return n
    s2 = settle(0.02*final); ADOPT = 200
    print("unique runs=%d, final mean=%.2f, within +/-2%% by run %d" % (n, final, s2))
    fig, ax = plt.subplots(figsize=(9.8, 5.4))
    ax.plot(range(1, n+1), Z, color="#C7CEE8", lw=0.6, alpha=0.8, label="tardiness of each individual run")
    ax.plot(range(1, n+1), run, color="#3B34A0", lw=2.3, label="cumulative running mean")
    ax.axhline(final, color="#999", ls=":", lw=1.1)
    ax.axvline(ADOPT, color="#333", ls="--", lw=1.5)
    ax.annotate("running mean stable to within \u00b12%%\nby ~%d runs; %d runs adopted" % (s2, ADOPT),
                xy=(ADOPT, final+2.4), xytext=(ADOPT+45, final+12), fontsize=11,
                arrowprops=dict(arrowstyle="->", lw=1.2))
    ax.set_xlabel("Number of simulation runs"); ax.set_ylabel("Mean floor tardiness Z of Model 2 (days)")
    ax.set_title("Stabilisation of the Mean Result with Number of Simulation Runs\n"
                 "Model 2 (GA reactive) at the extreme severity level "
                 "(3 stoppages + 3 material shortages + 3 reworks)", fontsize=11)
    ax.grid(alpha=0.25); ax.set_ylim(final-18, final+22); ax.legend(frameon=False, loc="upper right")
    plt.tight_layout(); plt.savefig("fig_runcount_convergence.png", dpi=180, bbox_inches="tight")
    print("saved fig_runcount_convergence.png")

if __name__ == "__main__":
    if len(sys.argv) >= 2 and sys.argv[1] == "plot":
        make_plot(); sys.exit(0)
    start_i, end_i = (int(sys.argv[1]), int(sys.argv[2])) if len(sys.argv) >= 3 else (0, TOTAL)
    if start_i == 0: open(path, "w").close()
    # regenerate the SAME scenario stream every time (fixed master seed) so runs are reproducible
    rng = random.Random(MASTER)
    scenarios = [extreme_scenario(rng) for _ in range(TOTAL)]
    t0 = time.time()
    with open(path, "a") as f:
        for i in range(start_i, end_i):
            w, rms = scenarios[i]
            seq = g.reactive(w, rms, "ga")
            Z = g.total_tardiness(g.run_combined(seq, w)[0])
            f.write(json.dumps({"i": i, "Z": round(Z, 3)}) + "\n"); f.flush()
    print("chunk %d-%d done in %.0fs" % (start_i, end_i, time.time()-t0))
