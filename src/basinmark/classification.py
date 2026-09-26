"""Classify fingerprint changes using shared relaxation and validated checkpoints."""

from __future__ import annotations

import hashlib
import json
import math
import time
from dataclasses import dataclass
from pathlib import Path

from .core import fingerprint_structure
from .io import JsonStore
from .relaxation import CalculatorRegistry, RelaxationEngine

CHECKPOINT_EVERY = 25


@dataclass(frozen=True)
class ClassificationConfig:
    engine: str
    potential_full: Path
    pairs_meta: Path
    cifs_dir: Path
    output: Path
    symprec: float = 0.1
    cell_bin: float = 0.05
    fmax: float = 0.02

    def __post_init__(self):
        if any(
            not math.isfinite(value) or value <= 0
            for value in (self.symprec, self.cell_bin, self.fmax)
        ):
            raise ValueError(
                "Symmetry tolerance, cell-bin width, and force threshold "
                "must be positive and finite"
            )


@dataclass(frozen=True)
class FlipComparison:
    """Compare discrete and full fingerprints without invoking a potential."""

    disc_mlip: str
    disc_dft: str
    full_mlip: str
    full_dft: str

    def as_record(self):
        discrete = int(self.disc_mlip != self.disc_dft)
        full = int(self.full_mlip != self.full_dft)
        if discrete and not full:
            raise ValueError("A discrete fingerprint change must also change the full fingerprint")
        return dict(
            recomputed_full_flip=full,
            discrete_flip=discrete,
            continuous_only_flip=int(not discrete and full),
            disc_mlip=self.disc_mlip,
            disc_dft=self.disc_dft,
            full_mlip=self.full_mlip,
            full_dft=self.full_dft,
        )


class FlipClassifier:
    """Own input validation, relaxation, classification, and resumable output."""

    def __init__(self, config: ClassificationConfig):
        self.config = config
        self.store = JsonStore(config.output)
        self.relaxer = RelaxationEngine(CalculatorRegistry(), config.engine, config.fmax)

    def _load_inputs(self):
        config = self.config
        if not config.cifs_dir.is_dir():
            raise ValueError("--cifs-dir must name an existing directory")
        with config.pairs_meta.open(encoding="utf-8") as handle:
            metadata = json.load(handle)
        if not isinstance(metadata, list):
            raise ValueError("Pair metadata must be a JSON list")
        meta = {}
        for record in metadata:
            if not isinstance(record, dict) or not isinstance(record.get("pair_id"), str):
                raise ValueError("Each metadata record must have a string pair_id")
            pair_id = record["pair_id"]
            if pair_id in meta:
                raise ValueError(f"Duplicate metadata pair_id: {pair_id}")
            meta[pair_id] = record
        full = JsonStore(config.potential_full).read_object()
        if config.engine not in full or not isinstance(full[config.engine], dict):
            raise ValueError("Potential input does not contain the requested engine")
        per_pair = full[config.engine].get("per_pair", [])
        if not isinstance(per_pair, list) or not per_pair:
            raise ValueError(
                "Potential input must contain nonempty per_pair records, not an aggregate summary"
            )
        ids = set()
        flipped = []
        structure_paths = set()
        for record in per_pair:
            if not isinstance(record, dict) or not isinstance(record.get("pair_id"), str):
                raise ValueError("Each potential record must have a string pair_id")
            pair_id = record["pair_id"]
            if pair_id in ids or record.get("flipped") not in (0, 1):
                raise ValueError(f"Duplicate pair_id or invalid flipped flag: {pair_id}")
            ids.add(pair_id)
            if not record["flipped"]:
                continue
            if pair_id not in meta:
                raise ValueError(f"No structure metadata for flipped pair: {pair_id}")
            for field in ("raw", "dft"):
                if not isinstance(meta[pair_id].get(field), str):
                    raise ValueError(f"Missing {field} structure path for pair: {pair_id}")
                path = config.cifs_dir / Path(meta[pair_id][field]).name
                if not path.is_file():
                    raise FileNotFoundError(f"Missing structure: {path}")
                structure_paths.add(path)
            flipped.append(record)
        digest = hashlib.sha256()
        for path in sorted(structure_paths):
            digest.update(path.name.encode() + b"\0")
            digest.update(hashlib.sha256(path.read_bytes()).digest())
        signature = dict(
            engine=config.engine,
            symprec=config.symprec,
            cell_bin=config.cell_bin,
            fmax=config.fmax,
            steps=self.relaxer.steps,
            potential_sha256=hashlib.sha256(config.potential_full.read_bytes()).hexdigest(),
            metadata_sha256=hashlib.sha256(config.pairs_meta.read_bytes()).hexdigest(),
            structures_sha256=digest.hexdigest(),
        )
        return meta, per_pair, flipped, signature

    def _load_checkpoint(self, signature, valid_ids):
        if not self.config.output.exists():
            return []
        previous = self.store.read_object()
        if previous.get("_configuration") != signature:
            raise ValueError(
                "Checkpoint inputs or settings differ, or metadata is absent; choose a new output"
            )
        records = previous.get("classified", [])
        if not isinstance(records, list):
            raise ValueError("Checkpoint classified records must be a list")
        ids = set()
        required = {
            "pair_id",
            "recomputed_full_flip",
            "discrete_flip",
            "continuous_only_flip",
            "disc_mlip",
            "disc_dft",
            "full_mlip",
            "full_dft",
        }
        for record in records:
            if not isinstance(record, dict) or not required.issubset(record):
                raise ValueError("Checkpoint contains an incomplete classification record")
            pair_id = record["pair_id"]
            if pair_id not in valid_ids or pair_id in ids:
                raise ValueError("Checkpoint contains duplicate or unrelated pair IDs")
            comparison = FlipComparison(
                *(record[key] for key in ("disc_mlip", "disc_dft", "full_mlip", "full_dft"))
            )
            expected = comparison.as_record()
            if any(record[key] != value for key, value in expected.items()):
                raise ValueError("Checkpoint classification is inconsistent with its fingerprints")
            ids.add(pair_id)
        return records

    def classify_pair(self, record, metadata):
        from pymatgen.core import Structure

        config = self.config
        raw = Structure.from_file(config.cifs_dir / Path(metadata["raw"]).name)
        dft = Structure.from_file(config.cifs_dir / Path(metadata["dft"]).name)
        relaxed = self.relaxer.run(raw).structure
        comparison = FlipComparison(
            fingerprint_structure(
                relaxed,
                symprec=config.symprec,
                cell_bin=config.cell_bin,
                use_cell_shape=False,
            ),
            fingerprint_structure(
                dft,
                symprec=config.symprec,
                cell_bin=config.cell_bin,
                use_cell_shape=False,
            ),
            fingerprint_structure(
                relaxed,
                symprec=config.symprec,
                cell_bin=config.cell_bin,
                use_cell_shape=True,
            ),
            fingerprint_structure(
                dft,
                symprec=config.symprec,
                cell_bin=config.cell_bin,
                use_cell_shape=True,
            ),
        )
        fields = comparison.as_record()
        if not fields["recomputed_full_flip"]:
            print(
                f"  WARNING: {record['pair_id']} no longer flips; retaining the mismatch",
                flush=True,
            )
        return dict(
            pair_id=record["pair_id"],
            generator=record.get("generator"),
            crystal_system=record.get("crystal_system"),
            potential_comparison_flipped=record["flipped"],
            **fields,
        )

    def run(self):
        meta, per_pair, flipped, signature = self._load_inputs()
        classified = self._load_checkpoint(signature, {record["pair_id"] for record in flipped})
        done_ids = {record["pair_id"] for record in classified}
        print(
            f"{self.config.engine}: {len(per_pair)} pairs, {len(flipped)} flips, "
            f"{len(done_ids)} already classified",
            flush=True,
        )
        started = time.time()
        completed = 0

        def write_checkpoint():
            output = dict(
                classified=classified,
                engine=self.config.engine,
                n_total_pairs=len(per_pair),
                n_flipped_per_potential_comparison=len(flipped),
                done=len(classified),
                complete=len(classified) == len(flipped),
                _configuration=signature,
            )
            self.store.write(output)
            return output

        for record in flipped:
            pair_id = record["pair_id"]
            if pair_id in done_ids:
                continue
            classified.append(self.classify_pair(record, meta[pair_id]))
            done_ids.add(pair_id)
            completed += 1
            if completed % CHECKPOINT_EVERY == 0:
                print(
                    f"  {self.config.engine}: {len(classified)}/{len(flipped)} classified, "
                    f"{time.time() - started:.0f}s elapsed",
                    flush=True,
                )
                write_checkpoint()
        return write_checkpoint()
