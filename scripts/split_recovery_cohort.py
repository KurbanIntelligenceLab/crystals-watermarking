"""Split a reserved physical cohort by source order and verify disjoint producer classes."""

import argparse
import hashlib
from pathlib import Path

from basinmark.io import JsonStore


def split(source):
    original = JsonStore(source).read_object()
    rows = original["rows"]
    if len(rows) != 768 or len({r["id"] for r in rows}) != len(rows):
        raise ValueError("expected 768 distinct reserved source IDs")
    first, last = rows[:384], rows[384:]

    def fingerprint(row):
        return row.get("producer", {}).get("no_cell", row["id"])

    left = {fingerprint(r) for r in first}
    cross = [r for r in last if fingerprint(r) in left]
    last = [r for r in last if fingerprint(r) not in left]
    # Keep shared classes together by exchanging singleton first-half classes for crossovers.
    if cross:
        from collections import Counter

        counts = Counter(fingerprint(r) for r in rows)
        exchange = [r for r in reversed(first) if counts[fingerprint(r)] == 1][: len(cross)]
        if len(exchange) != len(cross):
            raise ValueError("cannot produce equal class-disjoint groups")
        move = {r["id"] for r in exchange}
        first = [r for r in first if r["id"] not in move] + cross
        last += exchange
    assert len(first) == len(last) == 384
    assert not ({fingerprint(r) for r in first} & {fingerprint(r) for r in last})
    return [
        {
            "metadata": {
                **original["metadata"],
                "sample_size": 384,
                "half": i,
                "crossovers": len(cross),
                "parent_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
            },
            "rows": part,
        }
        for i, part in enumerate((first, last), 1)
    ]


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", type=Path, required=True)
    p.add_argument("--out", type=Path, nargs=2, required=True)
    a = p.parse_args()
    if a.out[0] == a.out[1] or any(path.exists() for path in a.out):
        p.error("choose two distinct new outputs")
    for path, part in zip(a.out, split(a.source), strict=True):
        JsonStore(path).write(part)
