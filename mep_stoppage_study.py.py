# ==========================================================================
# MEP STOPPAGE STUDY  -  standalone, single-file simulation sweep
#
# WHAT IT DOES
#   Reproduces the full dataset behind the MEP-stoppage results in the
#   dissertation. The line is simulated under MEP stoppages at five frequency
#   levels (1..5 stoppages per run) and three timing bands (early / middle /
#   late), 200 runs per level = 1000 runs. For each run three scheduling
#   strategies are scored on floor tardiness:
#       passive  - no rescheduling
#       model1   - rule-based reactive (earliest due date)
#       model2   - genetic-algorithm reactive
#
#   The design is 67 early / 67 middle / 66 late per level. It is produced in
#   two parts, exactly as in the study, both handled automatically here:
#     * a main sweep (seed 7): 67 early + 67 middle + 20 late per level;
#     * a late-band top-up (seed 20260714): 46 further late runs per level,
#       with any scenario identical to a main-sweep run rejected and redrawn.
#   The two are concatenated per level to give the 200-run design.
#
# HOW TO RUN
#   Needs only Python 3 (standard library only - nothing to install).
#       python stoppage_study.py
#   Runs the full 1000-run sweep (about 75-90 min) and writes
#   stoppage_results.csv next to this file. To run in shorter chunks:
#       python stoppage_study.py 0 200
#       python stoppage_study.py 200 400   ... up to 1000
#
# OUTPUT  (stoppage_results.csv, one row per run)
#   run, k_stoppages, timing_band, first_start_day,
#   passive_tardiness, model1_tardiness, model2_tardiness
#
# Running this script regenerates stoppage_results.csv exactly; that file is
# the dataset analysed in the dissertation's stoppage study.
# ==========================================================================
import random, json, sys, time, os
MIN = 480
PROC = {"A":[67,178,155,44], "B":[225,599,524,150], "C":[353,941,823,235]}
NF = 6; FLOORP = ["A"]*8 + ["B"]*2 + ["C"]*7
def baseline_sequence(): return [(f,t) for f in range(1,NF+1) for t in FLOORP]

def run(sequence, windows=None):
    """preemptive: MEP fully off during a window; in-progress job pauses & resumes."""
    windows = windows or []; ws = sorted(windows, key=lambda w: w["start"])
    n = len(sequence); C = [[0.0]*4 for _ in range(n)]; mep_start = {}
    def mep_finish(start, proc):
        t=start; moved=True
        while moved:
            moved=False
            for w in ws:
                if w["start"]<=t<w["end"]: t=w["end"]; moved=True; break
        es=t; rem=proc
        for w in ws:
            if w["end"]<=t: continue
            if t<w["start"]:
                avail=w["start"]-t
                if avail>=rem: return es, t+rem
                rem-=avail; t=w["end"]
            else: t=w["end"]
        return es, t+rem
    for j in range(n):
        f,typ=sequence[j]
        for s in range(4):
            ps=C[j-1][s] if j>0 else 0.0; pt=C[j][s-1] if s>0 else 0.0
            start=ps if ps>pt else pt
            if s==1: mep_start[j],C[j][s]=mep_finish(start,PROC[typ][s])
            else: C[j][s]=start+PROC[typ][s]
    ff={}
    for j in range(n): ff[sequence[j][0]]=max(ff.get(sequence[j][0],0.0),C[j][3])
    return [ff[f]/MIN for f in range(1,NF+1)], mep_start

TARGET,_ = run(baseline_sequence())
def total_tardiness(fl): return sum(max(0.0,f-t) for f,t in zip(fl,TARGET))

def split_at(sequence, windows, day):
    _, mep = run(sequence, windows)
    lo=[sequence[j] for j in range(len(sequence)) if mep[j]<day*MIN-1e-9]
    fr=[sequence[j] for j in range(len(sequence)) if mep[j]>=day*MIN-1e-9]
    return lo, fr

def edd_order(fp): return sorted(range(len(fp)), key=lambda i:(fp[i][0], i))

# ---- GA: population 100, generations 250, elitism 8, tournament 5 ----
def ga(fp, locked, windows, gens=250, pop=100, seed=1):
    m=len(fp)
    if m<=1:                                    # 0 or 1 modules: only one possible order
        order=list(range(m))
        return order, total_tardiness(run(locked+[fp[i] for i in order], windows)[0])
    def fit(o): return total_tardiness(run(locked+[fp[i] for i in o], windows)[0])
    def seeds():
        ed=edd_order(fp); S=[list(range(m)), ed, ed[::-1]]
        for _ in range(6): r=list(range(m)); random.shuffle(r); S.append(r)
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
        nP=[o for _,o in sc[:8]]                       # elitism = 8
        while len(nP)<pop:
            t1=min(random.sample(sc,5),key=lambda x:x[0])[1]   # tournament 5
            t2=min(random.sample(sc,5),key=lambda x:x[0])[1]
            c=ox(t1,t2)
            if random.random()<0.6: c=mut(c)                    # per-individual mutation 0.6
            nP.append(c)
        P=nP
    return best, bt

base = baseline_sequence()

def evaluate(windows):
    windows=sorted(windows,key=lambda w:w["start"])
    lock_day=windows[0]["start"]/MIN
    locked_m, free_m = split_at(base, [], lock_day)
    t_pass=total_tardiness(run(base, windows)[0])
    seq_m1=locked_m+[free_m[i] for i in edd_order(free_m)]
    t_m1=total_tardiness(run(seq_m1, windows)[0])
    # Reactive: GA re-plans at each stoppage (knows only stoppages seen so far)
    cur_locked, cur_free = locked_m, free_m
    known=[windows[0]]
    b1,_=ga(cur_free, cur_locked, known)
    cur_seq=cur_locked+[cur_free[i] for i in b1]
    for w in windows[1:]:
        known=known+[w]
        cur_locked, cur_free = split_at(cur_seq, known[:-1], w["start"]/MIN)
        bn,_=ga(cur_free, cur_locked, known)
        cur_seq=cur_locked+[cur_free[i] for i in bn]
    t_react=total_tardiness(run(cur_seq, windows)[0])
    return dict(passive=t_pass, m1=t_m1, react=t_react)

# ---- stratified scenario: first stoppage in [lo,hi]; the rest at/after it ----
def scenario_stratified(rng, k, lo, hi):
    first_day = rng.uniform(lo, hi)
    days = [first_day] + [rng.uniform(first_day, 117) for _ in range(k-1)]
    days.sort()
    ws=[]
    for d in days:
        d=round(d,1); dur=round(rng.triangular(0.25,3.0,0.5),2)
        ws.append(dict(station=1, start=d*MIN, end=(d+dur)*MIN))
    return ws

LEVELS = [1,2,3,4,5]
# exact counts per level: (band name, count, day range)
BANDS = [("early",100,(1,30)), ("middle",80,(31,90)), ("late",20,(91,117))]
SEED = 7


# ==========================================================================
# STUDY DRIVER  -  builds the 67/67/66 design from two seeds and writes CSV
# ==========================================================================
MAIN_SEED   = SEED            # 7, the main-sweep stream
TOPUP_SEED  = 20260714        # independent stream for the extra late runs
NEW_LATE_PER_LEVEL = 46       # 20 main-late + 46 top-up = 66 late per level


def _fingerprint(ws):
    return tuple(sorted((round(w["start"], 4), round(w["end"], 4)) for w in ws))


def build_design():
    """Return the full ordered list of (k, band, windows), 1000 scenarios,
       identical to the dataset analysed in the study."""
    # --- main stream (seed 7): 100 early / 80 middle / 20 late per level ---
    rng_main = random.Random(MAIN_SEED)
    main = {}
    existing_late = {k: set() for k in LEVELS}
    for k in LEVELS:
        for bname, count, (lo, hi) in BANDS:
            lst = main.setdefault((k, bname), [])
            for _ in range(count):
                ws = scenario_stratified(rng_main, k, lo, hi)
                lst.append(ws)
                if bname == "late":
                    existing_late[k].add(_fingerprint(ws))

    # --- top-up stream (seed 20260714): 46 new late per level, no duplicates ---
    rng_top = random.Random(TOPUP_SEED)
    topup = {k: [] for k in LEVELS}
    for k in LEVELS:
        seen = set(existing_late[k])
        made = 0
        while made < NEW_LATE_PER_LEVEL:
            ws = scenario_stratified(rng_top, k, 91, 117)
            fp = _fingerprint(ws)
            if fp in seen:
                continue
            seen.add(fp)
            topup[k].append(ws)
            made += 1

    # --- assemble the 67/67/66 design, per level ---
    design = []
    for k in LEVELS:
        for ws in main[(k, "early")][:67]:
            design.append((k, "early", ws))
        for ws in main[(k, "middle")][:67]:
            design.append((k, "middle", ws))
        for ws in main[(k, "late")][:20]:
            design.append((k, "late", ws))
        for ws in topup[k][:NEW_LATE_PER_LEVEL]:
            design.append((k, "late", ws))
    return design


if __name__ == "__main__":
    import csv
    design = build_design()
    N = len(design)                                  # 1000

    if len(sys.argv) >= 3:
        start_i, end_i = int(sys.argv[1]), int(sys.argv[2])
    else:
        start_i, end_i = 0, N

    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "stoppage_results.csv")
    header_row = ["run", "k_stoppages", "timing_band", "first_start_day",
                  "passive_tardiness", "model1_tardiness", "model2_tardiness"]
    mode = "w" if start_i == 0 else "a"

    t0 = time.time()
    with open(path, mode, newline="") as f:
        w = csv.writer(f)
        if start_i == 0:
            w.writerow(header_row)
        for i in range(start_i, end_i):
            k, band, ws = design[i]
            r = evaluate(ws)
            first_start = min(x["start"] for x in ws) / MIN
            w.writerow([i + 1, k, band, round(first_start, 1),
                        round(r["passive"], 3), round(r["m1"], 3), round(r["react"], 3)])
            f.flush()
    print("runs %d-%d done in %.0fs  (total: %d, output: stoppage_results.csv)"
          % (start_i, end_i, time.time() - t0, N))
