"""Ensure release archives do not disclose the builder's local identity."""

import io
import runpy
import tarfile
from pathlib import Path


def test_source_archive_clears_ownership_without_changing_contents(tmp_path):
    archive = tmp_path / "source.tar.gz"
    payload = b"release contents\n"
    with tarfile.open(archive, "w:gz") as source:
        member = tarfile.TarInfo("package/README.md")
        member.size = len(payload)
        member.uid, member.gid = 501, 20
        member.uname, member.gname = "example-builder", "example-group"
        member.mtime = 1234567890
        member.pax_headers = {"mtime": "1234567890.25"}
        source.addfile(member, io.BytesIO(payload))
    script = Path(__file__).resolve().parents[1] / "scripts/build_release.py"
    runpy.run_path(str(script))["anonymize_source_archive"](archive)
    with tarfile.open(archive) as source:
        member = source.getmember("package/README.md")
        assert source.extractfile(member).read() == payload
        assert (member.uid, member.gid, member.uname, member.gname) == (0, 0, "", "")
        assert member.mtime == 0
        assert not member.pax_headers
    header = archive.read_bytes()[:10]
    assert header[3] == 0  # No original filename or other optional gzip header fields.
    assert header[4:8] == b"\0" * 4
