# ==========================================================================
# MATERIAL SHORTAGE STUDY  -  standalone, single-file simulation sweep
#
# WHAT IT DOES
#   Reproduces the full dataset behind the material-shortage results in the
#   dissertation. A shortage blocks a whole module type from entering the MEP
#   bottleneck for the duration of a window (the material for that type is
#   unavailable). The line is simulated at five frequency levels (1..5
#   shortages per run) and three timing bands (early / middle / late),
#   200 runs per level = 1000 runs. For each run three scheduling strategies
#   are scored on floor tardiness and makespan:
#       passive  - no rescheduling
#       model1   - rule-based reactive (earliest due date)
#       model2   - genetic-algorithm reactive
#
#   The design is 67 early / 67 middle / 66 late per level, produced in two
#   parts, exactly as in the study, both handled automatically here:
#     * a main sweep (seed 7): 67 early + 67 middle + 20 late per level;
#     * a late-band top-up (seed 20260714): 46 further late runs per level,
#       any scenario identical to a main-sweep run rejected and redrawn.
#   The two are concatenated per level to give the 200-run design.
#
# HOW TO RUN
#   Needs only Python 3 (standard library only - nothing to install).
#       python shortage_study.py
#   Runs the full 1000-run sweep (about 75-90 min) and writes
#   shortage_results.csv next to this file. To run in shorter chunks:
#       python shortage_study.py 0 200
#       python shortage_study.py 200 400   ... up to 1000
#
# OUTPUT  (shortage_results.csv, one row per run; tardiness & makespan in days)
#   run, k_shortages, timing_band, first_start_day, shortage_types,
#   passive_tardiness, model1_tardiness, model2_tardiness,
#   passive_makespan, model1_makespan, model2_makespan
#
# Running this script regenerates shortage_results.csv exactly; that file is
# the dataset analysed in the dissertation's material-shortage study.
# ==========================================================================
import random, json, sys, time, os
MIN = 480
PROC = {"A":[67,178,155,44], "B":[225,599,524,150], "C":[353,941,823,235]}
NF = 6; FLOORP = ["A"]*8 + ["B"]*2 + ["C"]*7
def baseline_sequence(): return [(f,t) for f in range(1,NF+1) for t in FLOORP]

def run(sequence, windows=None):
    """Material delay: a blocked-type module cannot START at MEP during its window;
       it is pushed to the window end. Other types proceed. No overtaking within
       the given sequence (permutation flow shop)."""
    windows = windows or []
    n = len(sequence); C = [[0.0]*4 for _ in range(n)]; mep_start = {}
    for j in range(n):
        floor, typ = sequence[j]
        for s in range(4):
            prev_seq = C[j-1][s] if j>0 else 0.0
            prev_stn = C[j][s-1] if s>0 else 0.0
            start = prev_seq if prev_seq>prev_stn else prev_stn
            if s == 1:
                moved = True
                while moved:
                    moved = False
                    for w in windows:
                        if w["type"]==typ and w["start"] <= start < w["end"]:
                            start = w["end"]; moved = True; break
                mep_start[j] = start
            C[j][s] = start + PROC[typ][s]
    ff = {}
    for j in range(n):
        f = sequence[j][0]; ff[f] = max(ff.get(f,0.0), C[j][3])
    return [ff[f]/MIN for f in range(1,NF+1)], mep_start

TARGET, _ = run(baseline_sequence())
def total_tardiness(fl): return sum(max(0.0, f-t) for f,t in zip(fl, TARGET))

def split_at(sequence, windows, day):
    _, mep = run(sequence, windows)
    lo = [sequence[j] for j in range(len(sequence)) if mep[j] < day*MIN - 1e-9]
    fr = [sequence[j] for j in range(len(sequence)) if mep[j] >= day*MIN - 1e-9]
    return lo, fr

# ---- Model 1: greedy EDD dispatch that skips currently-blocked modules ----
def greedy_edd_skip(fp, windows, t0):
    t=t0; rem=list(range(len(fp))); order=[]
    blk=lambda i,tn: any(fp[i][1]==w["type"] and w["start"]<=tn<w["end"] for w in windows)
    while rem:
        av=[i for i in rem if not blk(i,t)]
        if not av:
            fe=[w["end"] for w in windows if w["end"]>t]
            if fe: t=min(fe); continue
            av=rem[:]
        p=min(av,key=lambda i:(fp[i][0],i)); order.append(p); rem.remove(p); t+=PROC[fp[p][1]][1]
    return order

# ---- GA: population 100, generations 250, elitism 8, tournament 5 ----
def ga(fp, locked, windows, t0, gens=250, pop=100, seed=1):
    m=len(fp)
    if m<=1:
        order=list(range(m))
        return order, total_tardiness(run(locked+[fp[i] for i in order], windows)[0])
    def fit(o): return total_tardiness(run(locked+[fp[i] for i in o], windows)[0])
    def seeds():
        fs=sorted(range(m),key=lambda i:(fp[i][0],i))
        S=[list(range(m)),fs,greedy_edd_skip(fp,windows,t0)]
        nonC=[i for i in fs if fp[i][1] not in [w["type"] for w in windows]]
        for K in (4,10,18):
            h=nonC[:K] if K<=len(nonC) else nonC
            S.append(h+[i for i in fs if i not in set(h)])
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

# ---- baseline MEP occupancy: which type is at the bottleneck at a given day ----
_, _BASE_MEP = run(base, [])
_BASE_INTERVALS = []
for _j in range(len(base)):
    _s = _BASE_MEP[_j]; _e = _s + PROC[base[_j][1]][1]
    _BASE_INTERVALS.append((_s, _e, base[_j][1]))
def type_at_mep(day):
    """Type the BASELINE schedule has entering/at MEP at 'day' (just-in-time delivery target).
       If the moment falls in a tiny gap, use the nearest following module's type."""
    t = day * MIN
    nxt = None
    for (s, e, ty) in _BASE_INTERVALS:
        if s <= t < e:
            return ty
        if s >= t and (nxt is None or s < nxt[0]):
            nxt = (s, ty)
    return nxt[1] if nxt else _BASE_INTERVALS[-1][2]

def evaluate(windows):
    windows = sorted(windows, key=lambda w: w["start"])
    lock_day = windows[0]["start"]/MIN
    locked_m, free_m = split_at(base, [], lock_day)
    flP = run(base, windows)[0]; t_pass = total_tardiness(flP)
    m1 = greedy_edd_skip(free_m, windows, windows[0]["start"])
    fl1 = run(locked_m+[free_m[i] for i in m1], windows)[0]; t_m1 = total_tardiness(fl1)
    # Reactive: GA re-plans at each delay (knows only delays seen so far)
    cur_locked, cur_free = locked_m, free_m
    known = [windows[0]]
    b1, _ = ga(cur_free, cur_locked, known, windows[0]["start"])
    cur_seq = cur_locked + [cur_free[i] for i in b1]
    for w in windows[1:]:
        known = known + [w]
        cur_locked, cur_free = split_at(cur_seq, known[:-1], w["start"]/MIN)
        bn, _ = ga(cur_free, cur_locked, known, w["start"])
        cur_seq = cur_locked + [cur_free[i] for i in bn]
    flR = run(cur_seq, windows)[0]; t_react = total_tardiness(flR)
    return dict(passive=t_pass, m1=t_m1, react=t_react,
                ms_passive=round(max(flP),3), ms_m1=round(max(fl1),3), ms_react=round(max(flR),3))

# ---- stratified scenario: first delay in [lo,hi]; the rest at/after it ----
def _next_type_mep_day(typ, day):
    """Baseline MEP day of the next module of 'typ' reaching MEP on/after 'day'."""
    t=day*MIN
    cands=[_BASE_MEP[j] for j in range(len(base)) if base[j][1]==typ and _BASE_MEP[j]>=t-1e-9]
    return (min(cands)/MIN) if cands else None

def _types_avail(day):
    t=day*MIN
    return sorted({base[j][1] for j in range(len(base)) if _BASE_MEP[j]>=t-1e-9})

def scenario_stratified(rng, k, lo, hi):
    first_day = rng.uniform(lo, hi)
    days = sorted([round(first_day,1)] + [round(rng.uniform(first_day,117),1) for _ in range(k-1)])
    ws=[]
    for d in days:
        avail=_types_avail(d) or ["C"]
        typ=rng.choice(avail)
        anchor=_next_type_mep_day(typ,d)      # window starts when next such module hits MEP
        if anchor is None: anchor=d
        dur=round(rng.triangular(0.25,3.0,0.5),2)
        ws.append(dict(type=typ, station=1, start=anchor*MIN, end=(anchor+dur)*MIN))
    return ws

LEVELS = [1,2,3,4,5]
BANDS = [("early",100,(1,30)), ("middle",80,(31,90)), ("late",20,(91,117))]
SEED = 7


# ==========================================================================
# STUDY DRIVER  -  builds the 67/67/66 design from two seeds and writes CSV
# ==========================================================================
MAIN_SEED  = SEED
TOPUP_SEED = 20260714
NEW_LATE_PER_LEVEL = 46


def _fingerprint(ws):
    # a shortage scenario's identity: each window's type, start and end
    return tuple(sorted((w["type"], round(w["start"], 4), round(w["end"], 4)) for w in ws))


def build_design():
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
    N = len(design)

    if len(sys.argv) >= 3:
        start_i, end_i = int(sys.argv[1]), int(sys.argv[2])
    else:
        start_i, end_i = 0, N

    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "shortage_results.csv")
    header_row = ["run", "k_shortages", "timing_band", "first_start_day", "shortage_types",
                  "passive_tardiness", "model1_tardiness", "model2_tardiness",
                  "passive_makespan", "model1_makespan", "model2_makespan"]
    mode = "w" if start_i == 0 else "a"

    t0 = time.time()
    with open(path, mode, newline="") as f:
        w = csv.writer(f)
        if start_i == 0:
            w.writerow(header_row)
        for i in range(start_i, end_i):
            k, band, ws = design[i]
            r = evaluate(ws)
            ws_sorted = sorted(ws, key=lambda x: x["start"])
            first_start = ws_sorted[0]["start"] / MIN
            types_str = "|".join(x["type"] for x in ws_sorted)
            w.writerow([i + 1, k, band, round(first_start, 1), types_str,
                        round(r["passive"], 3), round(r["m1"], 3), round(r["react"], 3),
                        round(r["ms_passive"], 3), round(r["ms_m1"], 3), round(r["ms_react"], 3)])
            f.flush()
    print("runs %d-%d done in %.0fs  (total: %d, output: shortage_results.csv)"
          % (start_i, end_i, time.time() - t0, N))
