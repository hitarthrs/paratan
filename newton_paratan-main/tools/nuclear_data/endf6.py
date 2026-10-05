"""Minimal ENDF 6 reader for the DD and DT fusion angular evaluations"""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import hashlib
import re
from typing import Sequence

_ENDF_FLOAT_RE = re.compile(r"^([+-]?(?:\d+(?:\.\d*)?|\.\d+))([+-]\d+)$")

@dataclass(frozen=True)
class EndfControlRecord:
    c1: float
    c2: float
    l1: int
    l2: int
    n1: int
    n2: int
    mat: int
    mf: int
    mt: int
    ns: int

@dataclass(frozen=True)
class EndfInterpolation:
    breakpoints: tuple[int, ...]
    laws: tuple[int, ...]

    def __post_init__(self) -> None:
        if len(self.breakpoints) != len(self.laws):
            raise ValueError("interpolation breakpoints and laws must have equal length")
        if not self.breakpoints:
            raise ValueError("at least one interpolation region is required")
        if any(value <= 0 for value in self.breakpoints):
            raise ValueError("interpolation breakpoints must be positive")
        if any(right <= left for left, right in zip(self.breakpoints, self.breakpoints[1:])):
            raise ValueError("interpolation breakpoints must increase")

@dataclass(frozen=True)
class EndfTab1:
    header: EndfControlRecord
    interpolation: EndfInterpolation
    x: tuple[float, ...]
    y: tuple[float, ...]

@dataclass(frozen=True)
class EndfTab2:
    header: EndfControlRecord
    interpolation: EndfInterpolation

@dataclass(frozen=True)
class EndfList:
    header: EndfControlRecord
    values: tuple[float, ...]

@dataclass(frozen=True)
class EndfMaterialHeader:
    mat: int
    za: float
    awr: float
    projectile_awr: float
    material_emax_eV: float
    library_release: int
    sublibrary: int
    format_version: int
    comments: tuple[str, ...]

@dataclass(frozen=True)
class EndfMf3Section:
    mat: int
    mt: int
    za: float
    awr: float
    mass_difference_Q_eV: float
    reaction_Q_eV: float
    breakup_flag: int
    interpolation: EndfInterpolation
    incident_energy_eV: tuple[float, ...]
    cross_section_barn: tuple[float, ...]

@dataclass(frozen=True)
class EndfLaw2Knot:
    incident_energy_eV: float
    lang: int
    values: tuple[float, ...]
    item_count: int

@dataclass(frozen=True)
class EndfMf6Product:
    zap: int
    awp: float
    lip: int
    law: int
    yield_interpolation: EndfInterpolation
    yield_incident_energy_eV: tuple[float, ...]
    yield_values: tuple[float, ...]
    angular_interpolation: EndfInterpolation | None
    angular_knots: tuple[EndfLaw2Knot, ...]

@dataclass(frozen=True)
class EndfMf6Section:
    mat: int
    mt: int
    za: float
    awr: float
    reference_frame: int
    products: tuple[EndfMf6Product, ...]

def parse_endf_float(field: str) -> float:
    text = field.strip()
    if not text:
        return 0.0
    normalized = text.replace("D", "E").replace("d", "e")
    if "e" in normalized.lower():
        return float(normalized)
    if re.fullmatch(r"[+-]?\d+(?:\.\d*)?", normalized):
        return float(normalized)
    match = _ENDF_FLOAT_RE.fullmatch(normalized)
    if match is None:
        raise ValueError(f"invalid ENDF floating point field {field!r}")
    return float(f"{match.group(1)}e{match.group(2)}")

def parse_endf_int(field: str) -> int:
    text = field.strip()
    return int(text) if text else 0

def parse_control_record(line: str) -> EndfControlRecord:
    if len(line) < 75:
        raise ValueError("ENDF record must contain at least 75 columns")
    fields = tuple(line[index : index + 11] for index in range(0, 66, 11))

    return EndfControlRecord(
        c1=parse_endf_float(fields[0]),
        c2=parse_endf_float(fields[1]),
        l1=parse_endf_int(fields[2]),
        l2=parse_endf_int(fields[3]),
        n1=parse_endf_int(fields[4]),
        n2=parse_endf_int(fields[5]),
        mat=int(line[66:70]),
        mf=int(line[70:72]),
        mt=int(line[72:75]),
        ns=int(line[75:80]),
    )

def parse_data_fields(line: str) -> tuple[float, ...]:
    if len(line) < 66:
        raise ValueError("ENDF data record must contain 66 data columns")
    return tuple(parse_endf_float(line[index : index + 11]) for index in range(0, 66, 11))

def read_endf_lines(path: str | Path) -> tuple[str, ...]:
    return tuple(Path(path).read_text(encoding="ascii").splitlines())

def section_lines(path: str | Path, *, mat: int, mf: int, mt: int) -> tuple[str, ...]:
    selected: list[str] = []
    for line in read_endf_lines(path):
        if len(line) < 75:
            continue
        try:
            line_mat = int(line[66:70])
            line_mf = int(line[70:72])
            line_mt = int(line[72:75])
        except ValueError:
            continue
        if (line_mat, line_mf, line_mt) == (mat, mf, mt):
            selected.append(line)
    if not selected:
        raise ValueError(f"ENDF section MAT={mat} MF={mf} MT={mt} is absent")
    
    return tuple(selected)

def _take_values(lines: Sequence[str], index: int, count: int) -> tuple[tuple[float, ...], int]:
    values: list[float] = []
    while len(values) < count:
        if index >= len(lines):
            raise ValueError("ENDF section ended before the requested values were read")
        values.extend(parse_data_fields(lines[index]))
        index += 1

    return tuple(values[:count]), index

def parse_tab1(lines: Sequence[str], index: int) -> tuple[EndfTab1, int]:
    header = parse_control_record(lines[index])
    index += 1
    interpolation_values, index = _take_values(lines, index, 2 * header.n1)
    interpolation = EndfInterpolation(breakpoints=tuple(int(round(value)) for value in interpolation_values[0::2]), laws=tuple(int(round(value)) for value in interpolation_values[1::2]))
    xy_values, index = _take_values(lines, index, 2 * header.n2)
    table = EndfTab1(header=header, interpolation=interpolation, x=tuple(xy_values[0::2]), y=tuple(xy_values[1::2]))

    return table, index

def parse_tab2(lines: Sequence[str], index: int) -> tuple[EndfTab2, int]:
    header = parse_control_record(lines[index])
    index += 1
    interpolation_values, index = _take_values(lines, index, 2 * header.n1)
    table = EndfTab2(
        header=header,
        interpolation=EndfInterpolation(breakpoints=tuple(int(round(value)) for value in interpolation_values[0::2]), laws=tuple(int(round(value)) for value in interpolation_values[1::2])),)
    
    return table, index

def parse_list(lines: Sequence[str], index: int) -> tuple[EndfList, int]:
    header = parse_control_record(lines[index])
    index += 1
    values, index = _take_values(lines, index, header.n1)

    return EndfList(header=header, values=values), index

def parse_material_header(path: str | Path, *, mat: int) -> EndfMaterialHeader:
    lines = section_lines(path, mat=mat, mf=1, mt=451)
    if len(lines) < 4:
        raise ValueError("MF=1 MT=451 section is incomplete")
    first = parse_control_record(lines[0])
    third = parse_control_record(lines[2])
    fourth = parse_control_record(lines[3])
    comment_count = fourth.n1
    if len(lines) < 4 + comment_count:
        raise ValueError("MF=1 MT=451 comment block is incomplete")
    comments = tuple(line[:66].rstrip() for line in lines[4 : 4 + comment_count])

    return EndfMaterialHeader(mat=mat, za=first.c1, awr=first.c2, projectile_awr=third.c1, material_emax_eV=third.c2, library_release=third.l1, sublibrary=third.n1, format_version=third.n2, comments=comments)

def parse_mf3_section(path: str | Path, *, mat: int, mt: int) -> EndfMf3Section:
    lines = section_lines(path, mat=mat, mf=3, mt=mt)
    head = parse_control_record(lines[0])
    table, index = parse_tab1(lines, 1)
    if index != len(lines):
        raise ValueError("unexpected trailing records in MF=3 section")
    if table.interpolation.breakpoints[-1] != len(table.x):
        raise ValueError("MF=3 interpolation does not cover all points")
    
    return EndfMf3Section(
        mat=mat,
        mt=mt,
        za=head.c1,
        awr=head.c2,
        mass_difference_Q_eV=table.header.c1,
        reaction_Q_eV=table.header.c2,
        breakup_flag=table.header.l1,
        interpolation=table.interpolation,
        incident_energy_eV=table.x,
        cross_section_barn=table.y,
    )

def parse_mf6_section(path: str | Path, *, mat: int, mt: int) -> EndfMf6Section:
    lines = section_lines(path, mat=mat, mf=6, mt=mt)
    head = parse_control_record(lines[0])
    index = 1
    products: list[EndfMf6Product] = []
    for _ in range(head.n1):
        yield_table, index = parse_tab1(lines, index)
        angular_interpolation: EndfInterpolation | None = None
        angular_knots: tuple[EndfLaw2Knot, ...] = ()
        if yield_table.header.l2 == 2:
            angular_table, index = parse_tab2(lines, index)
            knots: list[EndfLaw2Knot] = []
            for _ in range(angular_table.header.n2):
                angular_list, index = parse_list(lines, index)
                knots.append(EndfLaw2Knot(incident_energy_eV=angular_list.header.c2, lang=angular_list.header.l1, values=angular_list.values, item_count=angular_list.header.n2))
            angular_interpolation = angular_table.interpolation
            angular_knots = tuple(knots)
        elif yield_table.header.l2 != 4:
            raise ValueError(f"unsupported MF=6 LAW={yield_table.header.l2} for MAT={mat} MT={mt}")
        products.append(EndfMf6Product(
                zap=int(round(yield_table.header.c1)),
                awp=yield_table.header.c2,
                lip=yield_table.header.l1,
                law=yield_table.header.l2,
                yield_interpolation=yield_table.interpolation,
                yield_incident_energy_eV=yield_table.x,
                yield_values=yield_table.y,
                angular_interpolation=angular_interpolation,
                angular_knots=angular_knots,
            ))
    if index != len(lines):
        raise ValueError("unexpected trailing records in MF=6 section")
    
    return EndfMf6Section(mat=mat, mt=mt, za=head.c1, awr=head.c2, reference_frame=head.l2, products=tuple(products))