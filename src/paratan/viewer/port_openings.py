"""Small, local surface trims for port openings; no global mesh subdivision.

Only triangles near an implicit boundary receive extra vertices. The same trim
is applied after section caps are built, so cut faces also have real openings.
The port's own wall surfaces supply the exposed tunnel walls.
"""
from __future__ import annotations

import numpy as np
import pyvista as pv


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
        x, y, z = points.T
        # Both +/-45 degree tube axes lie in the x-z plane, through the origin.
        d = np.minimum(np.hypot(y, (x-z)/np.sqrt(2)), np.hypot(y, (x+z)/np.sqrt(2)))
        return np.maximum(d-self.radius, np.abs(x)-self.limit)

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
        if component.group in {"CC", "FW", "VV", "LF"} and component.name != "vacuum_vessel_plasma":
            return trim_surface(mesh, self.outside_ports, classifier=self.port_bounds)
        return mesh

    def port_bounds(self, triangles):
        x,y,z = triangles.transpose(2,0,1)
        bounds = [radial_bounds(np.stack([y,(x-sign*z)/np.sqrt(2)],axis=2)) for sign in (-1,1)]
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
