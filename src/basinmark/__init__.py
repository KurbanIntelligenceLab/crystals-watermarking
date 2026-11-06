"""Reusable tools for crystal watermark construction and auditing."""

from .allocation import select_collection
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
from .symmetry_audit import audit_deposited, paired_audit

__all__ = [
    "AuditResult",
    "audit",
    "audit_deposited",
    "audit_size",
    "critical_value",
    "fingerprint_structure",
    "is_green",
    "n_forge",
    "paired_audit",
    "power",
    "prf_uniform",
    "q_flipped",
    "q_green",
    "select_collection",
    "select_green",
]

__version__ = "0.5.0"
