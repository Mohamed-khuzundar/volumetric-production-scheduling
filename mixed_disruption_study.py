# ==========================================================================
# MIXED-DISRUPTION STUDY  -  standalone, single-file simulation sweep
#
# WHAT IT DOES
#   Reproduces the full dataset behind the mixed-disruption results in the
#   dissertation. All three disruption types strike in one run together:
#   MEP stoppages, material shortages, and reworks, in equal number. Three
#   severity levels are tested:
#       mild     = 1 of each type  (3 disruptions per run)
#       moderate = 2 of each type  (6 disruptions per run)
#       extreme  = 3 of each type  (9 disruptions per run)
#   Each level has 200 runs across three timing bands (early / middle / late),
#   67 / 67 / 66 per band = 600 runs total. Each run scores three scheduling
#   strategies on floor tardiness and makespan:
#       passive  - no rescheduling
#       model1   - rule-based reactive (earliest due date)
#       model2   - genetic-algorithm reactive, re-planning at each disruption
#
# HOW TO RUN
#   Needs only Python 3 (standard library only - nothing to install).
#       python mixed_study.py
#   Runs the full 600-run sweep (heavy: Model 2 re-plans up to 9 times per run)
#   and writes mixed_results.csv next to this file. To run in shorter chunks:
#       python mixed_study.py 0 100
#       python mixed_study.py 100 200   ... up to 600
#
# OUTPUT  (mixed_results.csv, one row per run; tardiness & makespan in days)
#   run, severity_level, disruptions_each_type, timing_band,
#   n_stoppages, n_shortages, n_reworks, realised_failures,
#   passive_tardiness, model1_tardiness, model2_tardiness,
#   passive_makespan, model1_makespan, model2_makespan
#
# Running this script regenerates mixed_results.csv exactly; that file is the
# dataset analysed in the dissertation's mixed-disruption study.
# ==========================================================================
import json, sys, time, os
import random
MIN=480
PROC={"A":[67,178,155,44],"B":[225,599,524,150],"C":[353,941,823,235]}
NF=6; FLOORP=["A"]*8+["B"]*2+["C"]*7

def baseline_items():
    out=[]; oid=0
    for f in range(1,NF+1):
        for t in FLOORP:
            out.append(dict(f=f,t=t,kind="N",rid=None,oid=oid)); oid+=1
    return out
base=baseline_items()

def _mep_push(start, typ, stop_ws, delay_ws):
    """push a desired MEP start later past any active stoppage or type-block window."""
    moved=True
    while moved:
        moved=False
        for w in stop_ws:                       # stoppage: blocks ALL types
            if w["start"]<=start<w["end"]: start=w["end"]; moved=True; break
        if moved: continue
        for w in delay_ws:                      # delay: blocks its type only
            if w["type"]==typ and w["start"]<=start<w["end"]: start=w["end"]; moved=True; break
    return start

def _mep_finish_with_stoppage(start, proc, stop_ws):
    """MEP processing that PAUSES during stoppage windows (preemptive), like the
       stoppage engine. Returns (effective_start, completion)."""
    t=start; moved=True
    while moved:
        moved=False
        for w in stop_ws:
            if w["start"]<=t<w["end"]: t=w["end"]; moved=True; break
    es=t; rem=proc
    for w in sorted(stop_ws,key=lambda x:x["start"]):
        if w["end"]<=t: continue
        if t<w["start"]:
            avail=w["start"]-t
            if avail>=rem: return es, t+rem
            rem-=avail; t=w["end"]
        else: t=w["end"]
    return es, t+rem

def run_mixed(seq, stop_ws=None, delay_ws=None):
    """seq: list of dict items (may include kind 'F' failed and 'R' repair).
       stop_ws: stoppage windows [{start,end}]; delay_ws: [{type,start,end}]."""
    stop_ws=stop_ws or []; delay_ws=delay_ws or []
    n=len(seq); C=[{} for _ in range(n)]; free=[0.0]*4; F_insp={}; mep={}
    for j in range(n):
        it=seq[j]; route=[1,3] if it["kind"]=="R" else [0,1,2,3]; prev=None
        for s in route:
            arr=prev if prev is not None else 0.0
            base_start=arr if arr>free[s] else free[s]
            if s==1:
                if it["kind"]=="R":                       # repair readiness
                    base_start=max(base_start, F_insp.get(it["rid"],0.0))
                # push past stoppage + type-block, then process WITH stoppage pause
                pushed=_mep_push(base_start, it["t"], stop_ws, delay_ws)
                es,c=_mep_finish_with_stoppage(pushed, PROC[it["t"]][s], stop_ws)
                mep[j]=es; free[s]=c; C[j][s]=c; prev=c
            else:
                c=base_start+PROC[it["t"]][s]; free[s]=c; C[j][s]=c; prev=c
        if it["kind"]=="F": F_insp[it["rid"]]=C[j][3]
    ff={}
    for j in range(n):
        it=seq[j]
        if it["kind"] in ("N","R"): ff[it["f"]]=max(ff.get(it["f"],0.0),C[j][3])
    return [ff[f]/MIN for f in range(1,NF+1)], mep

TARGET,_=run_mixed(base)
def total_tardiness(fl): return sum(max(0.0,f-t) for f,t in zip(fl,TARGET))
def makespan(fl): return max(fl)



# ---- baseline inspection times (to know when each module would fail) ----
def _baseline_insp():
    n=len(base); free=[0.0]*4; C=[{} for _ in range(n)]
    for j in range(n):
        it=base[j]; prev=None
        for s in range(4):
            arr=prev if prev is not None else 0.0
            start=arr if arr>free[s] else free[s]
            c=start+PROC[it["t"]][s]; free[s]=c; C[j][s]=c; prev=c
        pass
    return [C[j][3] for j in range(n)]
BASE_INSP=_baseline_insp()

def types_available(day):
    t=day*MIN; return sorted({base[j]["t"] for j in range(len(base)) if BASE_INSP[j]>=t-1e-9})

def _next_type_mep_day(typ, day):
    """delay anchor: day the next module of 'typ' reaches MEP in baseline, on/after day."""
    n=len(base); free=[0.0]*4; C=[{} for _ in range(n)]
    for j in range(n):
        it=base[j]; prev=None
        for s in range(4):
            arr=prev if prev is not None else 0.0
            start=arr if arr>free[s] else free[s]
            if s==1 and it["t"]==typ and start>=day*MIN-1e-9:
                return start/MIN
            c=start+PROC[it["t"]][s]; free[s]=c; C[j][s]=c; prev=c
    return None

# ---------------------------------------------------------------- #
# Scenario: k of each type, first disruption in the chosen band
# ---------------------------------------------------------------- #
def mixed_scenario(rng, k, lo, hi):
    first_day = rng.uniform(lo, hi)
    # all 3k disruption days: first is the band anchor, rest on/after it to 117
    n_total = 3*k
    days = sorted([round(first_day,1)] + [round(rng.uniform(first_day,117),1) for _ in range(n_total-1)])
    # assign the 3k days to types so that EACH type gets exactly k of them,
    # keeping chronological interleave: round-robin by sorted day index is biased,
    # so shuffle type-labels of a balanced multiset.
    labels = ["stop"]*k + ["delay"]*k + ["rework"]*k
    rng.shuffle(labels)
    stop_ws=[]; delay_ws=[]; rework_days=[]
    for d,lab in zip(days, labels):
        if lab=="stop":
            dur=round(rng.triangular(0.25,3.0,0.5),2)
            stop_ws.append(dict(start=d*MIN, end=(d+dur)*MIN))
        elif lab=="delay":
            avail=[t for t in types_available(d)] or ["C"]
            typ=rng.choice(avail)
            anchor=_next_type_mep_day(typ,d) or d
            dur=round(rng.triangular(0.25,3.0,0.5),2)
            delay_ws.append(dict(type=typ, start=anchor*MIN, end=(anchor+dur)*MIN))
        else:  # rework: choose a type available at d, fail its next module
            avail=[t for t in types_available(d)] or ["C"]
            typ=rng.choice(avail)
            rework_days.append((typ,d))
    return stop_ws, delay_ws, rework_days

# ---- resolve which modules fail, on the baseline (common to all models) ----
def failures_on_baseline(rework_days):
    used=set(); out=[]; rid=0
    remaining=list(rework_days)
    while remaining:
        best=None
        for idx,(typ,day) in enumerate(remaining):
            cands=[(BASE_INSP[p],p) for p in range(len(base))
                   if base[p]["oid"] not in used and base[p]["t"]==typ and BASE_INSP[p]>=day*MIN-1e-9]
            if cands:
                Tf,p=min(cands)
                if best is None or Tf<best[0]: best=(Tf,idx,p)
        if best is None: break
        Tf,idx,p=best; used.add(base[p]["oid"])
        out.append(dict(oid=base[p]["oid"],f=base[p]["f"],t=base[p]["t"],rid=rid,T=Tf))
        remaining.pop(idx); rid+=1
    return out

def edd_key(it): return (it["f"], {"F":0,"N":1,"R":2}[it["kind"]])

def _repair_after_F(items):
    """ensure each R sits immediately after its matching F (passive rule)."""
    out=items[:]
    while True:
        moved=False
        for i,it in enumerate(out):
            if it["kind"]=="R":
                fp=next((k for k,x in enumerate(out) if x["kind"]=="F" and x["rid"]==it["rid"]),None)
                if fp is not None and i<fp:
                    r=out.pop(i)
                    fp=next(k for k,x in enumerate(out) if x["kind"]=="F" and x["rid"]==it["rid"])
                    out.insert(fp+1,r); moved=True; break
        if not moved: return out

# ---------------- Passive ----------------
def passive_seq(fails):
    rid_of={d["oid"]:d["rid"] for d in fails}
    items=[]
    for j in range(len(base)):
        b=base[j]
        if b["oid"] in rid_of:
            items.append(dict(f=b["f"],t=b["t"],kind="F",rid=rid_of[b["oid"]],oid=b["oid"]))
            items.append(dict(f=b["f"],t=b["t"],kind="R",rid=rid_of[b["oid"]],oid=None))
        else:
            items.append(dict(b))
    return _repair_after_F(items)

# ---------------- GA over a free pool (full settings) ----------------
def ga_pool(locked, free, stop_ws, delay_ws, gens=250, pop=100, seed=1):
    m=len(free)
    if m<=1: return free[:]
    def fit(o): return total_tardiness(run_mixed(locked+[free[i] for i in o], stop_ws, delay_ws)[0])
    byfloor=sorted(range(m), key=lambda i: edd_key(free[i]))
    def seeds():
        S=[list(range(m)), byfloor]
        for _ in range(5):
            r=list(range(m)); random.shuffle(r); S.append(r)
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
        sc=sorted([(fit(o),o) for o in P], key=lambda x:x[0])
        if sc[0][0]<bt: bt,best=sc[0][0],sc[0][1][:]
        nP=[o for _,o in sc[:8]]
        while len(nP)<pop:
            t1=min(random.sample(sc,5),key=lambda x:x[0])[1]
            t2=min(random.sample(sc,5),key=lambda x:x[0])[1]
            c=ox(t1,t2)
            if random.random()<0.6: c=mut(c)
            nP.append(c)
        P=nP
    return [free[i] for i in best]

# ---------------- Reactive (EDD or GA), re-planning at each onset ----------------
def reactive_seq(fails, stop_ws, delay_ws, policy):
    """Re-plan at each disruption ONSET (stoppage start, delay start, or rework
       inspection time), revealing only disruptions up to that time."""
    # order of onsets
    onsets=[]
    for w in stop_ws: onsets.append(("stop", w["start"], w))
    for w in delay_ws: onsets.append(("delay", w["start"], w))
    for d in fails: onsets.append(("rework", d["T"], d))
    onsets.sort(key=lambda x:x[1])

    # start from planned items (no failures marked yet)
    order=[dict(base[j]) for j in range(len(base))]
    known_stop=[]; known_delay=[]; known_fail_oids=set()

    def current_seq_with_known_failures():
        rid_of={d["oid"]:d["rid"] for d in fails if d["oid"] in known_fail_oids}
        items=[]
        for it in order:
            if it["kind"] in ("N","F") and it["oid"] in rid_of:
                items.append(dict(f=it["f"],t=it["t"],kind="F",rid=rid_of[it["oid"]],oid=it["oid"]))
                items.append(dict(f=it["f"],t=it["t"],kind="R",rid=rid_of[it["oid"]],oid=None))
            elif it["kind"] in ("N","F"):
                items.append(dict(f=it["f"],t=it["t"],kind="N",rid=None,oid=it["oid"]))
        return _repair_after_F(items)

    seq=[dict(base[j]) for j in range(len(base))]
    for kind,tt,obj in onsets:
        if kind=="stop": known_stop.append(obj)
        elif kind=="delay": known_delay.append(obj)
        else: known_fail_oids.add(obj["oid"])
        # rebuild items reflecting known failures
        rid_of={d["oid"]:d["rid"] for d in fails if d["oid"] in known_fail_oids}
        items=[]
        for j in range(len(base)):
            b=base[j]
            if b["oid"] in rid_of:
                items.append(dict(f=b["f"],t=b["t"],kind="F",rid=rid_of[b["oid"]],oid=b["oid"]))
                items.append(dict(f=b["f"],t=b["t"],kind="R",rid=rid_of[b["oid"]],oid=None))
            else:
                items.append(dict(b))
        items=_repair_after_F(items)
        # split: lock everything already started at MEP before this onset time tt
        _,mep=run_mixed(items, known_stop, known_delay)
        locked=[]; free=[]
        for j,it in enumerate(items):
            (locked if mep[j]<tt-1e-9 else free).append(it)
        if policy=="ga":
            free=ga_pool(locked, free, known_stop, known_delay)
        else:  # edd
            free=_repair_after_F(sorted(free, key=edd_key))
        seq=locked+free
        order=seq  # carry forward
    # final run with ALL disruptions
    return seq

def evaluate_mixed(stop_ws, delay_ws, rework_days):
    fails=failures_on_baseline(rework_days)
    pS=passive_seq(fails)
    flP=run_mixed(pS, stop_ws, delay_ws)[0]
    s1=reactive_seq(fails, stop_ws, delay_ws, "edd")
    fl1=run_mixed(s1, stop_ws, delay_ws)[0]
    sR=reactive_seq(fails, stop_ws, delay_ws, "ga")
    flR=run_mixed(sR, stop_ws, delay_ws)[0]
    return dict(passive=total_tardiness(flP), m1=total_tardiness(fl1), react=total_tardiness(flR),
                ms_passive=round(max(flP),3), ms_m1=round(max(fl1),3), ms_react=round(max(flR),3),
                nfail=len(fails))

LEVELS=[("mild",1),("moderate",2),("extreme",3)]
BANDS=[("early",67,(1,30)),("middle",67,(31,90)),("late",66,(91,117))]
SEED=7

if __name__=="__main__":
    import csv
    rng=random.Random(SEED)
    scenarios=[]
    for lname,k in LEVELS:
        for bname,count,(lo,hi) in BANDS:
            for _ in range(count):
                scenarios.append((lname,k,bname,mixed_scenario(rng,k,lo,hi)))
    N=len(scenarios)
    a,b=(int(sys.argv[1]),int(sys.argv[2])) if len(sys.argv)>=3 else (0,N)

    path=os.path.join(os.path.dirname(os.path.abspath(__file__)),"mixed_results.csv")
    header_row=["run","severity_level","disruptions_each_type","timing_band",
                "n_stoppages","n_shortages","n_reworks","realised_failures",
                "passive_tardiness","model1_tardiness","model2_tardiness",
                "passive_makespan","model1_makespan","model2_makespan"]
    mode="w" if a==0 else "a"

    t0=time.time()
    with open(path,mode,newline="") as f:
        w=csv.writer(f)
        if a==0: w.writerow(header_row)
        for i in range(a,min(b,N)):
            lname,k,bname,(stop_ws,delay_ws,rework_days)=scenarios[i]
            r=evaluate_mixed(stop_ws,delay_ws,rework_days)
            w.writerow([i+1, lname, k, bname,
                        len(stop_ws), len(delay_ws), len(rework_days), r["nfail"],
                        round(r["passive"],3), round(r["m1"],3), round(r["react"],3),
                        round(r["ms_passive"],3), round(r["ms_m1"],3), round(r["ms_react"],3)])
            f.flush()
    print("runs %d-%d done in %.0fs  (total: %d, output: mixed_results.csv)"
          % (a, min(b,N), time.time()-t0, N))
