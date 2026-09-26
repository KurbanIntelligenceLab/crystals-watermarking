"""Build release distributions and remove local identity from source-archive headers."""

import gzip
import os
import subprocess
import sys
import tarfile
from pathlib import Path
from tempfile import NamedTemporaryFile


def anonymize_source_archive(path: Path) -> None:
    """Preserve member contents while clearing owner and timestamp metadata."""
    with NamedTemporaryFile(dir=path.parent, suffix=".tmp", delete=False) as handle:
        temporary = Path(handle.name)
    try:
        with tarfile.open(path, "r:gz") as source, temporary.open("wb") as raw:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
                with tarfile.open(
                    fileobj=compressed, mode="w", format=tarfile.PAX_FORMAT
                ) as target:
                    for member in source:
                        member.uid = member.gid = 0
                        member.uname = member.gname = ""
                        member.mtime = 0
                        member.pax_headers = {}
                        contents = source.extractfile(member) if member.isfile() else None
                        try:
                            target.addfile(member, contents)
                        finally:
                            if contents is not None:
                                contents.close()
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    subprocess.run([sys.executable, "-m", "build"], cwd=root, check=True)
    for archive in (root / "dist").glob("*.tar.gz"):
        anonymize_source_archive(archive)
        print(f"Cleared local archive metadata: {archive.name}")


if __name__ == "__main__":
    main()
