"""
============================================================================
 BASELINE SIMULATION  —  Volumetric Modular Construction Production Line
============================================================================
 Discrete-event simulation (SimPy) of a four-station flow shop.

 PURPOSE (model verification step):
   Reproduce, in SimPy, the deterministic makespan and floor-completion
   times that were independently calculated in the Excel flow-shop model.
   Matching that known result confirms the simulation logic is correct,
   before any variability or disruptions are added.

 CONFIGURATION:
   - 102 modules  =  6 floors  x  17 modules per floor (8 A, 2 B, 7 C)
   - Fixed initial sequence per floor: 8A - 2B - 7C (ascending size)
   - 4 stations in fixed order: Framing -> MEP -> Finishing -> Inspection
   - One module processed at a time per station (single position)
   - Processing times are FIXED per module type (Option B): the learning
     curve justifies the block sequence but is NOT applied to in-run times,
     so reordering under disruption never invalidates the times.
   - No variability (CV off) and no disruptions in this baseline.

 VERIFICATION TARGET (from the no-learning deterministic Excel model):
   Makespan       ~ 117.5 working days
   Floor finishes ~ [21.5, 40.7, 59.9, 79.1, 98.3, 117.5] days
============================================================================
"""

import simpy

# ---------------------------------------------------------------------------
# 1. FIXED INPUTS
# ---------------------------------------------------------------------------
MINUTES_PER_DAY = 480                      # one 8-hour working day

# Processing time (minutes) for each module type at each station.
# Station order:        Framing,  MEP,  Finishing,  Inspection
PROCESSING_TIME = {
    "A": [67,  178, 155, 44],              # Type A  (small,  4.3 m^2)
    "B": [225, 599, 524, 150],             # Type B  (medium, 7.9 m^2)
    "C": [353, 941, 823, 235],             # Type C  (large,  9.9 m^2)
}
STATION_NAMES = ["Framing", "MEP", "Finishing", "Inspection"]
N_STATIONS = len(STATION_NAMES)

# One floor = 8 Type A, then 2 Type B, then 7 Type C  (block order, ascending size)
MODULES_PER_FLOOR = ["A"] * 8 + ["B"] * 2 + ["C"] * 7      # 17 modules
N_FLOORS = 6                                               # 6 x 17 = 102


# ---------------------------------------------------------------------------
# 2. BUILD THE PRODUCTION SEQUENCE  (the 102 modules, in order)
# ---------------------------------------------------------------------------
def build_sequence():
    """Return the ordered list of 102 modules, each tagged with type and floor."""
    sequence = []
    for floor in range(1, N_FLOORS + 1):
        for module_type in MODULES_PER_FLOOR:
            sequence.append({"type": module_type, "floor": floor})
    return sequence


# ---------------------------------------------------------------------------
# 3. THE SIMULATION
# ---------------------------------------------------------------------------
def run_baseline():
    """Run the deterministic baseline and return completion results."""
    env = simpy.Environment()

    # Each station is a resource that can hold ONE module at a time.
    stations = [simpy.Resource(env, capacity=1) for _ in range(N_STATIONS)]

    sequence = build_sequence()

    results = {
        "module_finish": {},    # module index -> time it left Inspection
        "floor_finish": {},     # floor number -> time its LAST module finished
    }

    def module_flow(env, index, info):
        """A single module's journey: Framing -> MEP -> Finishing -> Inspection."""
        module_type = info["type"]
        for s in range(N_STATIONS):
            # Wait for this station to be free, then occupy it...
            with stations[s].request() as request:
                yield request
                # ...and spend this station's processing time.
                yield env.timeout(PROCESSING_TIME[module_type][s])

        # Module has cleared Inspection — record its finish time.
        finish_time = env.now
        results["module_finish"][index] = finish_time

        # A floor is "complete" when its slowest (last-finishing) module is done.
        floor = info["floor"]
        previous = results["floor_finish"].get(floor, 0)
        results["floor_finish"][floor] = max(previous, finish_time)

    # Launch all 102 module processes, in sequence order.
    # Because each station holds one module and serves first-come-first-served,
    # the modules keep their sequence order through the whole line (no overtaking).
    for index, info in enumerate(sequence):
        env.process(module_flow(env, index, info))

    env.run()    # run until every module has finished
    return results


# ---------------------------------------------------------------------------
# 4. REPORT THE RESULTS
# ---------------------------------------------------------------------------
def total_tardiness(floor_finish, target):
    """Total floor tardiness: sum over floors of (finish - target), late only.
       For the baseline the finishes equal the target, so this is exactly zero.
       The same measure is used to score every disruption scenario against it."""
    return sum(max(0.0, floor_finish[f] - target[f]) for f in range(1, N_FLOORS + 1))


def main():
    results = run_baseline()
    module_finish = results["module_finish"]
    floor_finish = results["floor_finish"]

    makespan_min = max(module_finish.values())

    # The baseline floor finishes ARE the on-time targets that every
    # disruption scenario is later measured against.
    target = dict(floor_finish)
    tardiness_min = total_tardiness(floor_finish, target)

    print("=" * 60)
    print(" BASELINE RESULT  (deterministic, no variability, no disruption)")
    print("=" * 60)
    print(f" Modules processed : {len(module_finish)}")
    print(f" Makespan          : {makespan_min / MINUTES_PER_DAY:>8.2f} working days")
    print(f" Total tardiness   : {tardiness_min / MINUTES_PER_DAY:>8.2f} days"
          f"   (baseline reference: on time by definition)")
    print()
    print(" Floor completion times")
    print(f"   {'Floor':<7}{'days':>9}")
    print("   " + "-" * 16)
    for floor in range(1, N_FLOORS + 1):
        m = floor_finish[floor]
        print(f"   {floor:<7}{m / MINUTES_PER_DAY:>9.2f}")
    print("=" * 60)


if __name__ == "__main__":
    main()
