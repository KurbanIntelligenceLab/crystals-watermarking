"""Validated input readers and atomic JSON output."""

import json
import os
from dataclasses import dataclass
from pathlib import Path
from tempfile import NamedTemporaryFile


@dataclass(frozen=True)
class JsonStore:
    path: Path

    def read_object(self):
        with self.path.open(encoding="utf-8") as handle:
            result = json.load(handle)
        if not isinstance(result, dict):
            raise ValueError(f"{self.path}: expected a JSON object")
        return result

    def write(self, value):
        self.path.parent.mkdir(parents=True, exist_ok=True)
        with NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=self.path.parent,
            prefix=".checkpoint-",
            suffix=".tmp",
            delete=False,
        ) as handle:
            temporary = Path(handle.name)
            try:
                json.dump(value, handle, indent=2, default=float, allow_nan=False)
            except BaseException:
                temporary.unlink(missing_ok=True)
                raise
        try:
            os.replace(temporary, self.path)
        finally:
            temporary.unlink(missing_ok=True)
        print(f"  wrote {self.path}", flush=True)


@dataclass(frozen=True)
class PairManifest:
    path: Path

    def load(self, fields: tuple[str, ...]):
        """Validate a nonempty pair manifest and resolve its structure paths."""
        path = self.path
        manifest = Path(path)
        with manifest.open(encoding="utf-8") as handle:
            records = json.load(handle)
        if not isinstance(records, list) or not records:
            raise ValueError(f"{path}: expected a nonempty JSON list of pair records")
        for index, record in enumerate(records):
            if not isinstance(record, dict):
                raise ValueError(f"{path}: record {index} must be an object")
            if "relaxed" in fields and "relaxed" not in record and "dft" in record:
                record["relaxed"] = record["dft"]
            for field in fields:
                if field not in record or not isinstance(record[field], str):
                    raise ValueError(f"{path}: record {index} requires a string '{field}' path")
                candidate = Path(record[field])
                if not candidate.is_file() and not candidate.is_absolute():
                    candidate = manifest.parent / candidate
                if not candidate.is_file():
                    raise FileNotFoundError(f"{path}: record {index}, '{field}': {record[field]}")
                record[field] = str(candidate)
        return records


@dataclass(frozen=True)
class FingerprintCorpus:
    path: Path

    def load(self) -> list[str]:
        records = self.path.read_text(encoding="utf-8").splitlines()
        if not records or any(not line.strip() for line in records):
            raise ValueError(f"{self.path}: expected one nonempty fingerprint per line")
        return [line.strip() for line in records]
