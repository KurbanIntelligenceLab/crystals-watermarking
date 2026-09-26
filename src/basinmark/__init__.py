"""Reusable tools for crystal watermark construction and auditing."""

from .core import (
    AuditResult,
    audit,
    audit_size,
    critical_value,
    fingerprint_structure,
    is_green,
    n_forge,
    power,
    prf_uniform,
    q_flipped,
    q_green,
    select_green,
)

__all__ = [
    "AuditResult",
    "audit",
    "audit_size",
    "critical_value",
    "fingerprint_structure",
    "is_green",
    "n_forge",
    "power",
    "prf_uniform",
    "q_flipped",
    "q_green",
    "select_green",
]

__version__ = "0.2.0"
