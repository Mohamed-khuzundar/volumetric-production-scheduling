# ==========================================================================
# RESULTS ANALYSIS  -  one script for every study's results file
#
# WHAT IT DOES
#   Reads the results CSVs produced by the study scripts (stoppage, material
#   shortage, rework, mixed) and reproduces the analysis reported in the
#   dissertation:
#     * mean tardiness of Passive, Model 1 and Model 2 per condition
#     * percentage improvement (vs Passive, and Model 2 vs Model 1)
#     * win / tie / loss counts of Model 2 against Model 1 (paired runs)
#     * mean makespan per condition, when the file contains it
#   The study's two dissertation figures are saved next to the CSV
#   (<name>_fig_frequency.png and <name>_fig_timing.png), and the full
#   printed analysis is also saved to analysis_report.txt.
#
# HOW TO RUN  (the simple way)
#   Put this script in the same folder as the results CSVs and run it with
#   no arguments - it finds and analyses every results file automatically:
#       python analyze_results.py
#
#   Or analyse one specific file:
#       python analyze_results.py rework_results.csv
#
#   Requires Python 3. matplotlib is optional: without it the tables and the
#   report are still produced, only the figures are skipped.
# ==========================================================================
import csv, sys, os
import statistics as st

REPORT = []
def echo(s=""):
    print(s)
    REPORT.append(str(s))

def analyze(path):
    # ---------------- read the file and detect the study ----------------------
    rows = list(csv.DictReader(open(path)))
    if not rows:
        echo("no data rows found in " + path); return

    DETECT = [("k_stoppages", "MEP stoppage",      "stoppages"),
              ("k_shortages", "material shortage", "shortages"),
              ("k_reworks",   "rework",            "reworks"),
              ("severity_level", "mixed disruption", "severity")]
    for col, study, unit in DETECT:
        if col in rows[0]:
            GROUP, STUDY, UNIT = col, study, unit
            break
    else:
        echo(f"skipping {os.path.basename(path)}: not a recognised results file"); return

    HAS_MS = "passive_makespan" in rows[0]

    # condition order: numeric levels sorted, or mild/moderate/extreme for mixed
    levels = sorted({r[GROUP] for r in rows},
                    key=(lambda v: {"mild": 1, "moderate": 2, "extreme": 3}.get(v, 99))
                    if GROUP == "severity_level" else (lambda v: int(v)))

    f = lambda r, c: float(r[c])
    EPS = 1e-9

    # ---------------- per-condition table -------------------------------------
    echo("=" * 78)
    echo(f" {STUDY.upper()} STUDY  -  {len(rows)} runs, grouped by {GROUP}")
    echo("=" * 78)
    hdr = (f"{'cond':>8} {'n':>4} | {'Passive':>8} {'Model1':>8} {'Model2':>8} |"
           f" {'M1 vs P':>8} {'M2 vs P':>8} {'M2 vs M1':>8} | {'W':>4} {'T':>4} {'L':>4}")
    echo(hdr); echo("-" * len(hdr))

    tot_w = tot_t = tot_l = 0
    for lv in levels:
        g  = [r for r in rows if r[GROUP] == lv]
        P  = st.mean(f(r, "passive_tardiness") for r in g)
        M1 = st.mean(f(r, "model1_tardiness")  for r in g)
        M2 = st.mean(f(r, "model2_tardiness")  for r in g)
        impM1 = (P - M1) / P * 100 if P else 0.0
        impM2 = (P - M2) / P * 100 if P else 0.0
        imp21 = (M1 - M2) / M1 * 100 if M1 else 0.0
        w = sum(1 for r in g if f(r, "model2_tardiness") < f(r, "model1_tardiness") - EPS)
        l = sum(1 for r in g if f(r, "model2_tardiness") > f(r, "model1_tardiness") + EPS)
        t = len(g) - w - l
        tot_w += w; tot_t += t; tot_l += l
        echo(f"{str(lv):>8} {len(g):>4} | {P:8.2f} {M1:8.2f} {M2:8.2f} |"
              f" {impM1:7.1f}% {impM2:7.1f}% {imp21:7.1f}% | {w:>4} {t:>4} {l:>4}")

    echo("-" * len(hdr))
    echo(f"{'all':>8} {len(rows):>4} | {'':8} {'':8} {'':8} |"
          f" {'':8} {'':8} {'':8} | {tot_w:>4} {tot_t:>4} {tot_l:>4}")
    echo("\n(tardiness in days; improvement = reduction vs the named model;"
          "\n W/T/L = Model 2 wins / ties / losses against Model 1, paired per run)")

    # ---------------- makespan table (when present) ----------------------------
    if HAS_MS:
        echo("\nMean makespan (days):")
        echo(f"{'cond':>8} | {'Passive':>8} {'Model1':>8} {'Model2':>8}")
        for lv in levels:
            g = [r for r in rows if r[GROUP] == lv]
            echo(f"{str(lv):>8} | "
                  f"{st.mean(f(r,'passive_makespan') for r in g):8.2f} "
                  f"{st.mean(f(r,'model1_makespan')  for r in g):8.2f} "
                  f"{st.mean(f(r,'model2_makespan')  for r in g):8.2f}")
    else:
        echo("\n(no makespan columns in this file)")

    # ---------------- figures (the two dissertation figures) ------------------
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        echo("\n(matplotlib not installed - figures skipped; tables above are complete)")
        return

    C_P, C_M1, C_M2 = "#B3432B", "#E09C41", "#1D9E75"
    base = os.path.splitext(path)[0]
    mean = lambda c, g: st.mean(f(r, c) for r in g)
    groups = [[r for r in rows if r[GROUP] == lv] for lv in levels]
    mP  = [mean("passive_tardiness", g) for g in groups]
    m1  = [mean("model1_tardiness",  g) for g in groups]
    m2  = [mean("model2_tardiness",  g) for g in groups]
    imp = [(a - b) / a * 100 if a else 0.0 for a, b in zip(m1, m2)]

    if GROUP == "severity_level":
        xlab = [f"{str(lv).capitalize()}\n({int(groups[i][0]['disruptions_each_type'])*3} total)"
                for i, lv in enumerate(levels)]
        xtitle = "severity level"
    else:
        xlab = [str(lv) for lv in levels]
        xtitle = f"number of {UNIT} per run (k)"

    # --- Figure 1: frequency/severity lines + improvement bars ---
    x = range(len(levels))
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.2, 4.6), dpi=200)
    a1.fill_between(x, mP, m1, color="#F2E2C4", alpha=0.85,
                    label="recovered by the rule (Model 1)")
    a1.fill_between(x, m1, m2, color="#BFE6D8", alpha=0.9,
                    label="added by optimisation (Model 2)")
    a1.plot(x, mP, "-o", color=C_P,  lw=2.2, ms=7, label="Passive")
    a1.plot(x, m1, "-s", color=C_M1, lw=2.2, ms=7, label="Model 1 (EDD)")
    a1.plot(x, m2, "-^", color=C_M2, lw=2.2, ms=7, label="Model 2 (GA)")
    a1.set_xticks(list(x)); a1.set_xticklabels(xlab)
    a1.set_xlabel(xtitle); a1.set_ylabel("Mean floor tardiness (days)")
    a1.set_title("(a) Mean tardiness by " +
                 ("severity" if GROUP == "severity_level" else "disruption frequency"))
    a1.set_ylim(0, max(mP) * 1.15); a1.grid(alpha=0.3)
    a1.legend(fontsize=8.5, frameon=True)

    a2.bar(list(x), imp, 0.5, color=C_M2)
    for i, v in zip(x, imp):
        a2.annotate(f"{v:.1f}%", (i, v), ha="center", va="bottom",
                    fontsize=11, fontweight="bold", color="#1a6b50")
    a2.set_xticks(list(x)); a2.set_xticklabels(xlab)
    a2.set_xlabel(xtitle)
    a2.set_ylabel("Model 2 improvement over Model 1 (%)")
    a2.set_title("(b) Optimisation gain over the rule")
    a2.set_ylim(0, max(imp) * 1.25 if max(imp) > 0 else 1)
    a2.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    f1 = base + "_fig_frequency.png"
    fig.savefig(f1, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    echo(f"\nfigure saved: {f1}")

    # --- Figure 2: mean tardiness by timing band ---
    BANDS = [("early", "Early (days 1\u201330)"),
             ("middle", "Middle (days 31\u201390)"),
             ("late", "Late (days 91\u2013117)")]
    bg = [[r for r in rows if r["timing_band"] == b] for b, _ in BANDS]
    if all(len(g) > 0 for g in bg):
        bP  = [mean("passive_tardiness", g) for g in bg]
        b1  = [mean("model1_tardiness",  g) for g in bg]
        b2  = [mean("model2_tardiness",  g) for g in bg]
        fig, ax = plt.subplots(figsize=(7.6, 4.6), dpi=200)
        xb = range(3); w = 0.26
        for off, vals, col, lab in [(-w, bP, C_P, "Passive"),
                                    (0, b1, C_M1, "Model 1 (EDD)"),
                                    (w, b2, C_M2, "Model 2 (GA)")]:
            ax.bar([i + off for i in xb], vals, w, color=col, label=lab)
            for i, v in zip(xb, vals):
                ax.annotate(f"{v:.1f}", (i + off, v), ha="center", va="bottom", fontsize=9)
        ax.set_xticks(list(xb))
        ax.set_xticklabels([f"{lbl}\nn={len(g)}" for (b, lbl), g in zip(BANDS, bg)])
        ax.set_ylabel("Mean floor tardiness (days)")
        ax.set_title("Mean tardiness by timing of the first disruption")
        ax.legend(); ax.grid(axis="y", alpha=0.3)
        fig.tight_layout()
        f2 = base + "_fig_timing.png"
        fig.savefig(f2, dpi=200, bbox_inches="tight", facecolor="white")
        plt.close(fig)
        echo(f"figure saved: {f2}")

if __name__ == "__main__":
    here = os.path.dirname(os.path.abspath(__file__))
    if len(sys.argv) >= 2:
        targets = [sys.argv[1]]
    else:
        names = ["stoppage_results.csv", "shortage_results.csv",
                 "rework_results.csv", "mixed_results.csv"]
        targets = [os.path.join(here, n) for n in names
                   if os.path.exists(os.path.join(here, n))]
        if not targets:   # fall back: any *_results.csv beside the script
            targets = sorted(os.path.join(here, f) for f in os.listdir(here)
                             if f.endswith("_results.csv"))
        if not targets:
            sys.exit("no results CSV files found next to this script")
        echo(f"found {len(targets)} results file(s) to analyse\n")
    for t in targets:
        analyze(t)
        echo("")
    rp = os.path.join(here, "analysis_report.txt")
    open(rp, "w").write("\n".join(REPORT) + "\n")
    print("full report saved:", rp)
    # keep the window open when the script was double-clicked
    if len(sys.argv) < 2 and sys.stdin is not None and sys.stdin.isatty():
        try:
            input("\nPress Enter to close...")
        except EOFError:
            pass
