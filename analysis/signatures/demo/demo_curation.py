#!/usr/bin/env python3
"""Replay of the signature-curation workflow with no model behind it.

The rules are the ones in ../curate_signatures.py: overlap coefficient on
uppercased symbols, review of pairs above 0.50, candidate additions rejected
above 0.20 (exactly 0.20 passes), every change gated by a person's decision.
The model's replies and the person's decisions are read from toy_run.json, so
the run is deterministic, offline and needs no key. Standard library only.

    python demo_curation.py            # scripted replay
    python demo_curation.py --fast     # no pauses
    python demo_curation.py --interactive   # you make the decisions
"""
import argparse, copy, json, os, time

HERE = os.path.dirname(os.path.abspath(__file__))
DEDUP_THRESHOLD = 0.50
NEW_SIG_THRESHOLD = 0.20


def overlap_coefficient(a, b):
    if not a or not b:
        return 0.0
    return len(a & b) / min(len(a), len(b))


def up(genes):
    return set(str(g).upper() for g in genes)


def say(text="", pause=0.0):
    print(text, flush=True)
    time.sleep(pause)


def candidate_pairs(catalogue):
    """Computed once, before any review, as in the original script."""
    pairs = []
    for cell, sigs in catalogue.items():
        names = list(sigs)
        for i, s1 in enumerate(names):
            for s2 in names[i + 1:]:
                g1, g2 = up(sigs[s1]), up(sigs[s2])
                ov = overlap_coefficient(g1, g2)
                if ov > DEDUP_THRESHOLD:
                    pairs.append(dict(cell=cell, s1=s1, s2=s2, g1=g1, g2=g2, ov=ov))
    return pairs


def breakdown(p):
    shared, u1, u2 = sorted(p["g1"] & p["g2"]), sorted(p["g1"] - p["g2"]), sorted(p["g2"] - p["g1"])
    say(f"   shared ({len(shared)}): {', '.join(shared)}")
    say(f"   only in {p['s1']} ({len(u1)}): {', '.join(u1) or '(none) - subset'}")
    say(f"   only in {p['s2']} ({len(u2)}): {', '.join(u2) or '(none) - subset'}")


MENU = ("1 keep first only", "2 keep second only", "3 keep both", "4 merge into first", "5 merge into second")


def ask_person(prompt, valid):
    while True:
        c = input(prompt).strip().upper()
        if c in valid:
            return c


def phase1(catalogue, run, log, interactive, pace):
    say("=" * 72 + "\nPHASE 1: REDUNDANCY REVIEW\n" + "=" * 72, pace)
    pairs = candidate_pairs(catalogue)
    replies = {(r["cell_type"], r["sig1"], r["sig2"]): r for r in run["dedup_replies"]}
    for n, p in enumerate(pairs, 1):
        c, s1, s2 = p["cell"], p["s1"], p["s2"]
        if s1 not in catalogue[c] or s2 not in catalogue[c]:
            continue
        say(f"\n[{n}/{len(pairs)}] {c}: '{s1}' vs '{s2}'   overlap coefficient {p['ov']:.1%}", pace)
        breakdown(p)
        r = replies.get((c, s1, s2), {"choice": 3, "reasoning": "(no scripted reply)"})
        say(f"   scripted reply (no model called): choice {r['choice']} - {r['reasoning']}", pace)
        for m in MENU:
            say("     " + m)
        if interactive:
            choice = ask_person("   your decision (1-5, S to skip): ", set("12345S"))
            choice = None if choice == "S" else int(choice)
        else:
            choice = run["dedup_decisions"].get(f"{c}|{s1}|{s2}")
            say(f"   scripted person's decision: {choice}", pace)
        if choice == 1:
            del catalogue[c][s2]; log.append(f"DEDUP {c}: kept {s1}, deleted {s2}")
        elif choice == 2:
            del catalogue[c][s1]; log.append(f"DEDUP {c}: kept {s2}, deleted {s1}")
        elif choice == 3:
            log.append(f"DEDUP {c}: kept both")
        elif choice in (4, 5):
            keep, drop = (s1, s2) if choice == 4 else (s2, s1)
            catalogue[c][keep] = sorted(p["g1"] | p["g2"])
            del catalogue[c][drop]; log.append(f"DEDUP {c}: merged into {keep}")
        say(f"   -> {log[-1] if choice else 'skipped'}", pace)


def phase2(catalogue, run, log, interactive, pace):
    say("\n" + "=" * 72 + "\nPHASE 2: MISSING-SIGNATURE DISCOVERY\n" + "=" * 72, pace)
    for cell in list(catalogue):
        say(f"\n{cell}: existing signatures {', '.join(catalogue[cell])}", pace)
        for s in run["suggestions"].get(cell, []):
            name, genes = s["name"], sorted(s["genes"])
            if name in catalogue[cell]:
                continue
            new = up(genes)
            worst = max(overlap_coefficient(new, up(g)) for g in catalogue[cell].values())
            say(f"\n   suggestion '{name}' ({len(genes)} genes, scripted text): {', '.join(genes)}", pace)
            say(f"   highest overlap with an existing signature: {worst:.2f}  (limit 0.20, equal passes)", pace)
            if worst > NEW_SIG_THRESHOLD:
                say("   -> rejected automatically, never shown to the person", pace)
                continue
            if not (8 <= len(genes) <= 12):
                say(f"   note: {len(genes)} genes. The 8-12 limit is a prompt instruction; the script does not check it.", pace)
            if interactive:
                d = ask_person("   add this signature? [y/n]: ", {"Y", "N"}).lower()
            else:
                d = run["suggestion_decisions"].get(f"{cell}|{name}", "n")
                say(f"   scripted person's decision: {d}", pace)
            if d == "y":
                catalogue[cell][name] = genes; log.append(f"ADDED {cell}: {name}")
                say(f"   -> added '{name}'", pace)
            else:
                say("   -> not added", pace)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fast", action="store_true", help="no pauses")
    ap.add_argument("--interactive", action="store_true", help="make the decisions yourself")
    ap.add_argument("--out", help="write the refined catalogue to this JSON file")
    a = ap.parse_args()
    run = json.load(open(os.path.join(HERE, "toy_run.json")))
    catalogue, log = copy.deepcopy(run["catalogue"]), []
    pace = 0.0 if a.fast else 0.35
    say("Signature curation replay. Toy catalogue, scripted replies, no model, no network.\n", pace)
    before = sum(len(v) for v in catalogue.values())
    phase1(catalogue, run, log, a.interactive, pace)
    phase2(catalogue, run, log, a.interactive, pace)
    say("\n" + "=" * 72)
    say(f"signatures: {before} before, {sum(len(v) for v in catalogue.values())} after")
    for line in log:
        say("  " + line)
    if a.out:
        json.dump(catalogue, open(a.out, "w"), indent=2)
        say(f"saved {a.out}")


if __name__ == "__main__":
    main()
