"""Port meshes and local surface trims (NBI tubes + radial blanket inlet/outlet).

NBI ±45° tubes and FFH radial IO pipes are built here (same place
``PortOpenings.arms`` already lived). Host components are trimmed with implicit
boundaries so cut faces show real openings; the port walls supply the bore.
"""
from __future__ import annotations

import math
from typing import Any

import numpy as np
import pyvista as pv

from src.paratan.viewer.components import ComponentMesh, color_for_material, component_label


def openings_for(data: dict):
    """Return a port-opening helper for the device, or ``None`` if nothing to trim."""
    helpers = []
    ports = data.get("ports") or {}
    if ports.get("ports_present"):
        vv = data.get("vacuum_vessel") or {}
        if "central_cell" in vv and "end_plug" in vv:
            helpers.append(TandemPortOpenings(data))
        else:
            helpers.append(PortOpenings(data))
    io = data.get("blanket_inlet_outlet") or {}
    if io.get("present"):
        helpers.append(BlanketInletOutletOpenings(data))
    if not helpers:
        return None
    if len(helpers) == 1:
        return helpers[0]
    return CombinedOpenings(*helpers)


def _port_tube_distance(points, radius, limit, z0=0.0):
    """Signed distance to the +/-45° tube pair centered at ``z0``."""
    x, y, z = points.T
    zr = z - z0
    d = np.minimum(np.hypot(y, (x - zr) / np.sqrt(2)), np.hypot(y, (x + zr) / np.sqrt(2)))
    return np.maximum(d - radius, np.abs(x) - limit)


class PortOpenings:
    def __init__(self, data):
        ports = data["ports"]
        self.radius = float(ports["inner_radius"]) + sum(ports.get("port_layers_thicknesses", []))
        self.limit = float(ports.get("port_x_limits", 250))
        vv = data["vacuum_vessel"]
        mid = float(vv.get("axial_midplane", 0))
        a, b = float(vv["outer_axial_length"]) / 2, float(vv["central_axial_length"]) / 2
        self.z = np.array([mid-b-float(vv["left_bottleneck_length"]), mid-b, mid-a,
                           mid+a, mid+b, mid+b+float(vv["right_bottleneck_length"])])
        rc, rb = float(vv["central_radius"]), float(vv["bottleneck_radius"])
        self.r = np.array([rb, rb, rc, rc, rb, rb])
        self.stepped = vv.get("geometry_style", "conical") != "conical"
        self.mid, self.half, self.rc, self.rb = mid, b, rc, rb
        self.plasma_lipschitz = max(1., np.hypot(1., abs(rc-rb) / max(b-a, 1e-6)))

    def outside_ports(self, points):
        return _port_tube_distance(points, self.radius, self.limit, z0=self.mid)

    def outside_plasma(self, points):
        x, y, z = points.T
        if self.stepped:
            # Union of central and bottle cylinders, also covering the step.
            central = np.maximum(np.hypot(x, y)-self.rc, np.abs(z-self.mid)-self.half)
            bottle = np.maximum.reduce([np.hypot(x, y)-self.rb, self.z[0]-z, z-self.z[-1]])
            return np.minimum(central, bottle)
        return np.maximum.reduce([np.hypot(x,y)-np.interp(z,self.z,self.r), self.z[0]-z, z-self.z[-1]])

    def trim(self, mesh, component):
        if component.group == "ports":
            if component.meta.get("plasma_clipped"):
                return mesh
            return trim_surface(mesh, self.outside_plasma, 1. if self.stepped else self.plasma_lipschitz,
                                classifier=None if self.stepped else self.plasma_bounds)
        if component.group in {"CC", "FW", "VV", "LF", "EP"} and component.name != "vacuum_vessel_plasma":
            return trim_surface(mesh, self.outside_ports, classifier=self.port_bounds)
        return mesh

    def port_bounds(self, triangles):
        x,y,z = triangles.transpose(2,0,1)
        zr = z - self.mid
        bounds = [radial_bounds(np.stack([y,(x-sign*zr)/np.sqrt(2)],axis=2)) for sign in (-1,1)]
        outside = (bounds[0][0] > self.radius) & (bounds[1][0] > self.radius)
        outside |= (x.min(axis=1) > self.limit) | (x.max(axis=1) < -self.limit)
        inside = ((bounds[0][1] < self.radius) | (bounds[1][1] < self.radius)) & (np.abs(x).max(axis=1) < self.limit)
        return outside, inside

    def plasma_bounds(self, triangles):
        low, high = radial_bounds(triangles[:,:,:2])
        zlo,zhi = triangles[:,:,2].min(axis=1),triangles[:,:,2].max(axis=1)
        radii = np.stack([np.interp(zlo,self.z,self.r),np.interp(zhi,self.z,self.r)],axis=1)
        rmin,rmax = radii.min(axis=1),radii.max(axis=1)
        for z,r in zip(self.z,self.r):
            crosses=(zlo<=z)&(zhi>=z)
            rmin=np.where(crosses,np.minimum(rmin,r),rmin)
            rmax=np.where(crosses,np.maximum(rmax,r),rmax)
        return ((low > rmax) | (zlo > self.z[-1]) | (zhi < self.z[0]),
                (high < rmin) & (zlo >= self.z[0]) & (zhi <= self.z[-1]))

    def arms(self, inner, outer, angle, resolution=64):
        """Generate only the two exterior tube arms, with exact boundary rings."""
        tilt=np.radians(angle)
        axis=np.array([np.sin(tilt),0.,np.cos(tilt)])
        radial=np.array([np.cos(tilt),0.,-np.sin(tilt)])
        phi=np.arange(resolution)*2*np.pi/resolution
        circle=np.cos(phi)[:,None]*radial+np.sin(phi)[:,None]*np.array([0.,1.,0.])
        rings=[]
        faces=[]
        for sign in (-1,1):
            for radius in (inner,outer):
                base=radius*circle
                if (self.outside_plasma(base)>0).any():
                    raise ValueError("Port radius/position exceeds the plasma interior; exterior-arm preview is unsupported")
                lo=np.zeros(resolution)
                hi=np.full(resolution,(self.limit+outer*abs(np.cos(tilt)))/abs(np.sin(tilt))+1)
                for _ in range(32):
                    mid=(lo+hi)/2
                    inside=self.outside_plasma(base+sign*mid[:,None]*axis)<0
                    lo=np.where(inside,mid,lo)
                    hi=np.where(inside,hi,mid)
                start=base+sign*((lo+hi)/2)[:,None]*axis
                end_t=(sign*self.limit*np.sign(axis[0])-base[:,0])/axis[0]
                end=base+end_t[:,None]*axis
                rings.extend([start,end])
            # Four rings: inner-start, inner-end, outer-start, outer-end.
            offset=(len(rings)-4)*resolution
            for a,b in ((0,1),(1,3),(3,2),(2,0)):
                for i in range(resolution):
                    j=(i+1)%resolution
                    faces.extend([4,offset+a*resolution+i,offset+a*resolution+j,
                                  offset+b*resolution+j,offset+b*resolution+i])
        mesh=pv.PolyData(np.vstack(rings),np.array(faces)).triangulate().clean()
        return mesh.compute_normals(auto_orient_normals=True,consistent_normals=True,split_vertices=False)

    def trim_outline(self, lines):
        """Retain design outlines, with short gaps where they cross openings."""
        points=[]
        segments=[]
        cells=lines.lines
        cursor=0
        while cursor<len(cells):
            count=int(cells[cursor])
            ids=cells[cursor+1:cursor+1+count]
            for a,b in zip(lines.points[ids[:-1]],lines.points[ids[1:]]):
                n=max(1,int(np.ceil(np.linalg.norm(b-a)/6)))
                row=a+np.linspace(0,1,n+1)[:,None]*(b-a)
                offset=sum(len(p) for p in points)
                points.append(row)
                segments.extend([n+1,*range(offset,offset+n+1)])
            cursor+=count+1
        result=pv.PolyData(np.vstack(points),lines=np.array(segments))
        result.point_data["port_trim"]=self.outside_ports(result.points)
        result=result.clip_scalar(scalars="port_trim",value=0.,invert=False)
        del result.point_data["port_trim"]
        return result


class TandemPortOpenings:
    """End-plug port openings: one +/-45° pair at each end-plug midplane."""

    def __init__(self, data):
        ports = data["ports"]
        self.radius = float(ports["inner_radius"]) + sum(ports.get("port_layers_thicknesses", []))
        self.limit = float(ports.get("port_x_limits", 400))
        vv = data["vacuum_vessel"]
        sep = float(vv["central_cell_end_plug_separation_distance"])
        mid = (
            float(vv["central_cell"]["central_axis_length"]) / 2
            + sep
            + float(vv["end_plug"]["central_axis_length"]) / 2
        )
        self.midplanes = {"left": -mid, "right": mid}
        self.stations = {
            key: PortOpenings(_end_plug_station_yaml(data, key, 0.0))
            for key in ("left", "right")
        }

    def outside_ports(self, points):
        return np.minimum.reduce(
            [
                _port_tube_distance(points, self.radius, self.limit, z0)
                for z0 in self.midplanes.values()
            ]
        )

    def outside_plasma(self, points, side):
        station = self.stations[side]
        shifted = np.asarray(points, dtype=float).copy()
        shifted[:, 2] -= self.midplanes[side]
        return station.outside_plasma(shifted)

    def trim(self, mesh, component):
        if component.group == "ports":
            if component.meta.get("plasma_clipped"):
                return mesh
            side = component.meta.get("plug_side")
            if side not in self.stations:
                return mesh
            station = self.stations[side]
            z0 = self.midplanes[side]
            local = mesh.copy(deep=True)
            local.points[:, 2] -= z0
            trimmed = trim_surface(
                local,
                station.outside_plasma,
                1.0 if station.stepped else station.plasma_lipschitz,
                classifier=None if station.stepped else station.plasma_bounds,
            )
            if trimmed.n_points:
                trimmed.points[:, 2] += z0
            return trimmed
        if component.group in {"FW", "VV", "LF", "EP"} and any(
            token in component.name for token in ("left", "right")
        ):
            return trim_surface(mesh, self.outside_ports, classifier=self.port_bounds)
        return mesh

    def port_bounds(self, triangles):
        outs, inns = [], []
        for z0 in self.midplanes.values():
            x, y, z = triangles.transpose(2, 0, 1)
            zr = z - z0
            bounds = [
                radial_bounds(np.stack([y, (x - sign * zr) / np.sqrt(2)], axis=2))
                for sign in (-1, 1)
            ]
            outside = (bounds[0][0] > self.radius) & (bounds[1][0] > self.radius)
            outside |= (x.min(axis=1) > self.limit) | (x.max(axis=1) < -self.limit)
            inside = (
                (bounds[0][1] < self.radius) | (bounds[1][1] < self.radius)
            ) & (np.abs(x).max(axis=1) < self.limit)
            outs.append(outside)
            inns.append(inside)
        return np.logical_and.reduce(outs), np.logical_or.reduce(inns)

    def arms(self, side, inner, outer, angle, resolution=64):
        mesh = self.stations[side].arms(inner, outer, angle, resolution=resolution)
        mesh = mesh.copy(deep=True)
        mesh.points[:, 2] += self.midplanes[side]
        return mesh

    def trim_outline(self, lines):
        points = []
        segments = []
        cells = lines.lines
        cursor = 0
        while cursor < len(cells):
            count = int(cells[cursor])
            ids = cells[cursor + 1 : cursor + 1 + count]
            for a, b in zip(lines.points[ids[:-1]], lines.points[ids[1:]]):
                n = max(1, int(np.ceil(np.linalg.norm(b - a) / 6)))
                row = a + np.linspace(0, 1, n + 1)[:, None] * (b - a)
                offset = sum(len(p) for p in points)
                points.append(row)
                segments.extend([n + 1, *range(offset, offset + n + 1)])
            cursor += count + 1
        result = pv.PolyData(np.vstack(points), lines=np.array(segments))
        result.point_data["port_trim"] = self.outside_ports(result.points)
        result = result.clip_scalar(scalars="port_trim", value=0.0, invert=False)
        del result.point_data["port_trim"]
        return result


def _end_plug_station_yaml(data: dict, side: str, local_midplane: float) -> dict:
    """Simple-mirror-shaped VV YAML for one end plug, in a local midplane frame."""
    vv = data["vacuum_vessel"]
    ep = vv["end_plug"]
    sep = float(vv["central_cell_end_plug_separation_distance"])
    if side == "left":
        left_len, right_len = float(vv["end_axial_distance"]), 0.9 * sep
    else:
        left_len, right_len = 0.9 * sep, float(vv["end_axial_distance"])
    return {
        "ports": data["ports"],
        "vacuum_vessel": {
            "axial_midplane": local_midplane,
            "outer_axial_length": ep["outer_axial_length"],
            "central_axial_length": ep["central_axis_length"],
            "left_bottleneck_length": left_len,
            "right_bottleneck_length": right_len,
            "central_radius": ep["central_radius"],
            "bottleneck_radius": vv["bottleneck radius"],
            "geometry_style": "conical",
        },
    }


def radial_bounds(triangles):
    """Exact min/max distance to origin for projected triangles (also slivers)."""
    a=triangles
    b=np.roll(a,-1,axis=1)
    edge=b-a
    t=np.clip(-np.sum(a*edge,axis=2)/np.maximum(np.sum(edge*edge,axis=2),1e-20),0,1)
    minimum=np.linalg.norm(a+t[:,:,None]*edge,axis=2).min(axis=1)
    cross=a[:,:,0]*b[:,:,1]-a[:,:,1]*b[:,:,0]
    area=cross.sum(axis=1)
    contains=((cross>=-1e-10).all(axis=1)|(cross<=1e-10).all(axis=1)) & (np.abs(area)>1e-10)
    return np.where(contains,0,minimum),np.linalg.norm(a,axis=2).max(axis=1)


def trim_surface(mesh, distance, lipschitz=1., edge_size=6., classifier=None):
    """Keep distance >= 0, resolving crossings only within a bounded local band.

    A Lipschitz bound catches holes even when all three original vertices are
    outside. Longest-edge bisection avoids uniformly refining whole components.
    """
    triangles = mesh.triangulate()
    indices = triangles.faces.reshape(-1,4)[:,1:]
    pending = triangles.points[indices].astype(float)
    normals = triangles.point_data.get("Normals")
    pending_normals = np.asarray(normals)[indices] if normals is not None else None
    kept, kept_normals = [], []
    changed = False
    for _ in range(24):
        if not len(pending):
            break
        centers = pending.mean(axis=1)
        extent = np.linalg.norm(pending-centers[:,None,:],axis=2).max(axis=1)
        values = distance(centers)
        outside = values > lipschitz*extent + 1e-7
        inside = values < -lipschitz*extent - 1e-7
        if classifier is not None:
            outside, inside = classifier(pending)
        lengths = np.stack([np.linalg.norm(pending[:,1]-pending[:,0],axis=1),
                            np.linalg.norm(pending[:,2]-pending[:,1],axis=1),
                            np.linalg.norm(pending[:,0]-pending[:,2],axis=1)],axis=1)
        small = lengths.max(axis=1) <= edge_size
        accept = outside | (small & ~inside)
        if accept.any():
            kept.append(pending[accept])
            if pending_normals is not None:
                kept_normals.append(pending_normals[accept])
        changed |= bool(inside.any() or (~outside).any())
        split = ~(accept | inside)
        if not split.any():
            pending = pending[:0]
            break
        p = pending[split]
        edge = lengths[split].argmax(axis=1)
        # Rotate each triangle so its longest edge is the first edge.
        order = (edge[:,None]+np.arange(3)) % 3
        p = np.take_along_axis(p,order[:,:,None],axis=1)
        mid = (p[:,0]+p[:,1])/2
        pending = np.concatenate([np.stack([p[:,0],mid,p[:,2]],axis=1),
                                  np.stack([mid,p[:,1],p[:,2]],axis=1)])
        if pending_normals is not None:
            n = np.take_along_axis(pending_normals[split],order[:,:,None],axis=1)
            nm = (n[:,0]+n[:,1])/2
            pending_normals = np.concatenate([np.stack([n[:,0],nm,n[:,2]],axis=1),
                                              np.stack([nm,n[:,1],n[:,2]],axis=1)])
        if len(pending) > 150_000:
            raise ValueError("Port opening exceeds the local mesh budget")
    if len(pending):
        raise ValueError("Port opening did not converge within the local mesh budget")
    if not changed:
        return mesh
    if not kept:
        return pv.PolyData()
    points = np.concatenate(kept).reshape(-1,3)
    faces = np.column_stack([np.full(len(points)//3,3),np.arange(len(points)).reshape(-1,3)]).ravel()
    result = pv.PolyData(points,faces)
    if kept_normals:
        result.GetPointData().SetNormals(pv.convert_array(np.concatenate(kept_normals).reshape(-1,3).astype(np.float32),name="Normals"))
    result.point_data["port_trim"] = distance(points)
    result = result.clip_scalar(scalars="port_trim",value=0.,invert=False)
    del result.point_data["port_trim"]
    return result


# ---------------------------------------------------------------------------
# Radial blanket inlet / outlet (FFH viz case)
# ---------------------------------------------------------------------------


def angular_distance_deg(a: float, b: float) -> float:
    """Smallest absolute angular distance between two angles [deg]."""
    return abs((a - b + 180.0) % 360.0 - 180.0)


def resolve_port_angles(
    n_ports: int,
    angles_deg: list[float] | None,
    nbi_forbidden_angles_deg: list[float] | tuple[float, ...] = (0.0, 180.0),
    nbi_avoid_margin_deg: float = 25.0,
) -> list[float]:
    """Resolve common inlet/outlet azimuths (degrees in [0, 360))."""
    if n_ports < 1:
        raise ValueError("n_ports must be >= 1.")
    spacing = 360.0 / n_ports
    forbidden = [float(a) % 360.0 for a in nbi_forbidden_angles_deg]

    if angles_deg is not None:
        if len(angles_deg) != n_ports:
            raise ValueError(
                f"angles_deg has {len(angles_deg)} values, but n_ports = {n_ports}."
            )
        angles = [float(a) % 360.0 for a in angles_deg]
    else:
        step_deg = 0.5
        n_scan = max(1, int(round(spacing / step_deg)))
        best_offset = 0.0
        best_clearance = -1.0
        for i in range(n_scan):
            offset = i * spacing / n_scan
            candidate = [(offset + j * spacing) % 360.0 for j in range(n_ports)]
            if forbidden:
                clearance = min(
                    angular_distance_deg(a, f) for a in candidate for f in forbidden
                )
            else:
                clearance = 180.0
            if clearance > best_clearance:
                best_clearance = clearance
                best_offset = offset
        angles = [(best_offset + j * spacing) % 360.0 for j in range(n_ports)]

    conflicts = []
    for angle in angles:
        for forbidden_a in forbidden:
            d = angular_distance_deg(angle, forbidden_a)
            if d < nbi_avoid_margin_deg:
                conflicts.append((angle, forbidden_a, d))
    if conflicts:
        lines = [
            "Blanket inlet/outlet angle conflicts with an NBI exclusion direction:"
        ]
        for angle, forbidden_a, d in conflicts:
            lines.append(
                f"  port {angle:.1f} deg vs NBI {forbidden_a:.1f} deg "
                f"(separation {d:.1f} deg < required {nbi_avoid_margin_deg:.1f} deg)"
            )
        raise ValueError("\n".join(lines))

    for i in range(len(angles)):
        for j in range(i + 1, len(angles)):
            if angular_distance_deg(angles[i], angles[j]) < 1.0e-6:
                raise ValueError("Duplicate angles_deg values detected.")
    return angles


def _pipe_frame(angle_deg: float, z_center_cm: float):
    """Orthonormal frame: axis along radial pipe, two perpendiculars for the circle."""
    theta = math.radians(angle_deg)
    axis = np.array([math.cos(theta), math.sin(theta), 0.0])
    radial = np.array([0.0, 0.0, 1.0])
    tangential = np.cross(axis, radial)
    tangential /= np.linalg.norm(tangential)
    return axis, radial, tangential, z_center_cm


def _circle_ring(
    axis: np.ndarray,
    radial: np.ndarray,
    tangential: np.ndarray,
    origin: np.ndarray,
    s: float,
    radius: float,
    resolution: int,
) -> np.ndarray:
    phi = np.arange(resolution) * 2.0 * np.pi / resolution
    center = origin + s * axis
    return (
        center
        + radius * np.cos(phi)[:, None] * radial
        + radius * np.sin(phi)[:, None] * tangential
    )


def _quad_tube(rings: list[np.ndarray], closed_caps: bool = False) -> pv.PolyData:
    """Connect consecutive rings with quads; optional end-disk caps on the first/last ring."""
    n = len(rings[0])
    points = np.vstack(rings)
    faces: list[int] = []
    for k in range(len(rings) - 1):
        a, b = k * n, (k + 1) * n
        for i in range(n):
            j = (i + 1) % n
            faces.extend([4, a + i, a + j, b + j, b + i])
    if closed_caps:
        base = len(points)
        points = np.vstack([points, rings[0].mean(axis=0), rings[-1].mean(axis=0)])
        c0, c1 = base, base + 1
        for i in range(n):
            j = (i + 1) % n
            faces.extend([3, c0, i, j])
            last = (len(rings) - 1) * n
            faces.extend([3, c1, last + j, last + i])
    mesh = pv.PolyData(points, np.asarray(faces, dtype=np.int64))
    return mesh.triangulate().compute_normals(
        auto_orient_normals=True, consistent_normals=True, split_vertices=False
    )


def _radial_solid_pipe(
    angle_deg: float,
    z_center_cm: float,
    r_start_cm: float,
    r_end_cm: float,
    radius_cm: float,
    resolution: int,
) -> pv.PolyData:
    """Filled finite cylinder directed radially outward in the xy-plane."""
    if r_end_cm <= r_start_cm:
        raise ValueError("pipe outer end must be beyond pipe_inner_end_radius_cm")
    if radius_cm <= 0.0:
        raise ValueError("pipe radius must be positive")
    axis, radial, tangential, z0 = _pipe_frame(angle_deg, z_center_cm)
    origin = np.array([0.0, 0.0, z0])
    res = max(16, resolution)
    rings = [
        _circle_ring(axis, radial, tangential, origin, r_start_cm, radius_cm, res),
        _circle_ring(axis, radial, tangential, origin, r_end_cm, radius_cm, res),
    ]
    return _quad_tube(rings, closed_caps=True)


def _radial_annular_pipe(
    angle_deg: float,
    z_center_cm: float,
    r_start_cm: float,
    r_end_cm: float,
    r_inner_cm: float,
    r_outer_cm: float,
    resolution: int,
) -> pv.PolyData:
    """Hollow pipe wall (outer barrel + inner bore + annular end caps)."""
    if r_outer_cm <= r_inner_cm:
        raise ValueError("pipe wall needs r_outer > r_inner")
    axis, radial, tangential, z0 = _pipe_frame(angle_deg, z_center_cm)
    origin = np.array([0.0, 0.0, z0])
    res = max(16, resolution)
    ro0 = _circle_ring(axis, radial, tangential, origin, r_start_cm, r_outer_cm, res)
    ro1 = _circle_ring(axis, radial, tangential, origin, r_end_cm, r_outer_cm, res)
    ri0 = _circle_ring(axis, radial, tangential, origin, r_start_cm, r_inner_cm, res)
    ri1 = _circle_ring(axis, radial, tangential, origin, r_end_cm, r_inner_cm, res)
    n = res
    points = np.vstack([ro0, ro1, ri0, ri1])
    faces: list[int] = []
    for i in range(n):
        j = (i + 1) % n
        faces.extend([4, i, j, n + j, n + i])
    o = 2 * n
    for i in range(n):
        j = (i + 1) % n
        faces.extend([4, o + i, o + n + i, o + n + j, o + j])
    for base_o, base_i in ((0, 2 * n), (n, 3 * n)):
        for i in range(n):
            j = (i + 1) % n
            faces.extend([4, base_o + i, base_o + j, base_i + j, base_i + i])
    mesh = pv.PolyData(points, np.asarray(faces, dtype=np.int64))
    return mesh.triangulate().compute_normals(
        auto_orient_normals=True, consistent_normals=True, split_vertices=False
    )


def build_blanket_inlet_outlet_meshes(
    data: dict[str, Any],
    *,
    n_theta: int = 48,
) -> list[ComponentMesh]:
    """YAML ``blanket_inlet_outlet`` → wall/flow ComponentMesh list (group IO)."""
    cfg = data.get("blanket_inlet_outlet") or {}
    if not cfg.get("present", False):
        return []

    n_ports = int(cfg["n_ports"])
    angles = resolve_port_angles(
        n_ports,
        cfg.get("angles_deg"),
        cfg.get("nbi_forbidden_angles_deg", (0.0, 180.0)),
        float(cfg.get("nbi_avoid_margin_deg", 25.0)),
    )
    z_min = float(cfg["blanket_z_min_cm"])
    z_max = float(cfg["blanket_z_max_cm"])
    r_start = float(cfg["pipe_inner_end_radius_cm"])
    r_end = float(cfg["blanket_outer_radius_cm"]) + float(cfg["pipe_outer_extension_cm"])
    r_flow = float(cfg["pipe_inner_radius_cm"])
    t_wall = float(cfg["pipe_wall_thickness_cm"])
    if r_flow <= 0.0 or t_wall <= 0.0:
        raise ValueError("pipe_inner_radius_cm and pipe_wall_thickness_cm must be > 0")
    r_wall = r_flow + t_wall
    inlet_z = z_min + r_wall
    outlet_z = z_max - r_wall
    wall_mat = str(cfg.get("wall_material", "stainless"))
    flow_mat = str(cfg.get("flow_material", "eutectic_breeding_material"))
    flow_opacity = float(cfg.get("flow_opacity", 0.85))
    resolution = max(24, n_theta // 2)

    components: list[ComponentMesh] = []
    for i, angle in enumerate(angles):
        for kind, zc in (("inlet", inlet_z), ("outlet", outlet_z)):
            idx = i + 1
            wall_mesh = _radial_annular_pipe(
                angle, zc, r_start, r_end, r_flow, r_wall, resolution
            )
            flow_mesh = _radial_solid_pipe(angle, zc, r_start, r_end, r_flow, resolution)
            for suffix, material, mesh, opacity in (
                ("wall", wall_mat, wall_mesh, 1.0),
                ("flow", flow_mat, flow_mesh, flow_opacity),
            ):
                name = f"{kind}_{idx}_{suffix}"
                meta = {
                    "kind": kind,
                    "index": i,
                    "angle_deg": angle,
                    "z0_cm": zc,
                    "r_inner_cm": 0.0 if suffix == "flow" else r_flow,
                    "r_outer_cm": r_flow if suffix == "flow" else r_wall,
                    "pipe_r_start_cm": r_start,
                    "pipe_r_end_cm": r_end,
                }
                components.append(
                    ComponentMesh(
                        name=name,
                        group="IO",
                        material=material,
                        mesh=mesh,
                        meta=meta,
                        color=color_for_material(material),
                        opacity=opacity,
                        profile=None,
                        label=component_label(name, "IO", material, meta),
                    )
                )
    return components


def pipe_specs(data: dict[str, Any]) -> list[dict[str, float]]:
    """Outer-pipe cylinders used to carve openings in blanket/VV meshes."""
    cfg = data.get("blanket_inlet_outlet") or {}
    if not cfg.get("present", False):
        return []
    angles = resolve_port_angles(
        int(cfg["n_ports"]),
        cfg.get("angles_deg"),
        cfg.get("nbi_forbidden_angles_deg", (0.0, 180.0)),
        float(cfg.get("nbi_avoid_margin_deg", 25.0)),
    )
    r_flow = float(cfg["pipe_inner_radius_cm"])
    r_wall = r_flow + float(cfg["pipe_wall_thickness_cm"])
    r_start = float(cfg["pipe_inner_end_radius_cm"])
    r_end = float(cfg["blanket_outer_radius_cm"]) + float(cfg["pipe_outer_extension_cm"])
    inlet_z = float(cfg["blanket_z_min_cm"]) + r_wall
    outlet_z = float(cfg["blanket_z_max_cm"]) - r_wall
    specs = []
    for angle in angles:
        for zc in (inlet_z, outlet_z):
            specs.append(
                {
                    "angle_deg": angle,
                    "z0_cm": zc,
                    "radius_cm": r_wall,
                    "r_start_cm": r_start,
                    "r_end_cm": r_end,
                }
            )
    return specs


class BlanketInletOutletOpenings:
    """Trim host meshes where radial inlet/outlet outer cylinders penetrate."""

    def __init__(self, data: dict[str, Any]):
        self.specs = pipe_specs(data)

    def outside_pipes(self, points: np.ndarray) -> np.ndarray:
        """Signed distance to the union of outer IO pipe cylinders (negative inside)."""
        if not self.specs:
            return np.full(len(points), np.inf)
        x, y, z = points.T
        d_union = np.full(len(points), np.inf)
        for spec in self.specs:
            theta = math.radians(spec["angle_deg"])
            ux, uy = math.cos(theta), math.sin(theta)
            s = x * ux + y * uy
            cx, cy = s * ux, s * uy
            radial = np.hypot(np.hypot(x - cx, y - cy), z - spec["z0_cm"])
            along = np.maximum(spec["r_start_cm"] - s, s - spec["r_end_cm"])
            d = np.maximum(radial - spec["radius_cm"], along)
            d_union = np.minimum(d_union, d)
        return d_union

    def pipe_bounds(self, triangles: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        pts = triangles.reshape(len(triangles), -1, 3)
        centers = pts.mean(axis=1)
        d = self.outside_pipes(centers)
        outside = (
            d > self.specs[0]["radius_cm"] + 20.0
            if self.specs
            else np.ones(len(d), bool)
        )
        inside = d < -1.0
        return outside, inside

    def trim(self, mesh: pv.PolyData, component: ComponentMesh) -> pv.PolyData:
        if not self.specs:
            return mesh
        if component.group == "IO":
            return mesh
        if component.group in {"CC", "FW", "VV", "LF"} and component.name != "vacuum_vessel_plasma":
            return trim_surface(mesh, self.outside_pipes, classifier=self.pipe_bounds)
        return mesh

    def trim_outline(self, lines: pv.PolyData) -> pv.PolyData:
        if not self.specs or lines is None or not lines.n_points:
            return lines
        result = lines.copy(deep=True)
        result.point_data["port_trim"] = self.outside_pipes(result.points)
        result = result.clip_scalar(scalars="port_trim", value=0.0, invert=False)
        if "port_trim" in result.point_data:
            del result.point_data["port_trim"]
        return result


class CombinedOpenings:
    """Apply NBI and blanket-IO trims in sequence."""

    def __init__(self, *helpers: Any):
        self.helpers = [h for h in helpers if h is not None]

    def trim(self, mesh: pv.PolyData, component: ComponentMesh) -> pv.PolyData:
        for helper in self.helpers:
            mesh = helper.trim(mesh, component)
        return mesh

    def trim_outline(self, lines: pv.PolyData) -> pv.PolyData:
        for helper in self.helpers:
            if hasattr(helper, "trim_outline"):
                lines = helper.trim_outline(lines)
        return lines

