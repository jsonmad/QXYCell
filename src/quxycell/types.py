"""Shared lightweight data structures."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class Message:
    """A check message emitted during project preflight."""

    level: str
    code: str
    message: str
    path: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class MeasurementFile:
    """Discovered QuPath measurement table."""

    path: Path
    delimiter: str
    n_columns: int
    columns: tuple[str, ...]
    n_rows: int | None = None

    def to_dict(self) -> dict[str, Any]:
        data = asdict(self)
        data["path"] = str(self.path)
        data["delimiter_name"] = "tab" if self.delimiter == "\t" else "comma"
        return data


@dataclass(frozen=True)
class ClassifierDefinition:
    """A simple single-measurement QuPath threshold classifier."""

    path: Path
    name: str
    measurement_column: str | None
    threshold: float | None
    is_simple: bool
    reason: str = ""

    def to_dict(self) -> dict[str, Any]:
        data = asdict(self)
        data["path"] = str(self.path)
        return data


@dataclass(frozen=True)
class GeoJsonFile:
    """Discovered GeoJSON file summary."""

    path: Path
    n_features: int | None
    object_type_counts: dict[str, int] = field(default_factory=dict)
    class_counts: dict[str, int] = field(default_factory=dict)
    name_counts: dict[str, int] = field(default_factory=dict)
    readable: bool = True
    error: str = ""

    def to_dict(self) -> dict[str, Any]:
        data = asdict(self)
        data["path"] = str(self.path)
        return data

