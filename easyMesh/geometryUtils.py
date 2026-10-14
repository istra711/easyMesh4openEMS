from CSXCAD import CSPrimitives, CSProperties
import openEMS
import numpy as np

def process_polygon(automesher, polygon, x, y, z, x_edges, y_edges, diagonal_edges, mesh_data):
    # Check if the input is a list of primitives
    if isinstance(polygon, list):
        # Process each primitive
        for prim in polygon:
            xedges, yedges = [], []
            process_primitive(automesher, prim, x, y, xedges, yedges, diagonal_edges)
            xedges.sort(key=lambda edge: edge[0])
            yedges.sort(key=lambda edge: edge[0])
            # Ensure no duplicate edges are added
            i = 0
            while i < len(xedges) - 1:
                if xedges[i][0] == xedges[i+1][0]:
                    del xedges[i]
                else:
                    i += 1
            i = 0
            while i < len(yedges) - 1:
                if yedges[i][0] == yedges[i+1][0]:
                    del yedges[i]
                else:
                    i += 1
            x_edges.extend(xedges)
            y_edges.extend(yedges)
        # Collect z-coordinates from the polygon
        z.extend(collect_z_coordinates(polygon))

        # Process z-coordinates and add them to the mesh_data list
        process_z_coordinates(automesher, z, mesh_data)

    else:
        # If the polygon is a single primitive, get its x and y coordinates
        # prim_x_coords, prim_y_coords = polygon.GetCoords()[0], polygon.GetCoords()[1]
        # x.extend(prim_x_coords)
        # y.extend(prim_y_coords)
        # Process the single polygon to extract edges and coordinates
        process_single_polygon(automesher, polygon, x, y, x_edges, y_edges, diagonal_edges)

        # Collect z-coordinates from the single polygon
        z.extend(collect_z_coordinates([polygon]))

        # Process z-coordinates and add them to the mesh_data
        process_z_coordinates(automesher, z, mesh_data)
    # Remove lines from mesh_data[0] or mesh_data[1] if they match x_edges[0] or y_edges[0]
    for x_edge in x_edges:
        if x_edge[0] in mesh_data[0]:
            mesh_data[0].remove(x_edge[0])
    for y_edge in y_edges:
        if y_edge[0] in mesh_data[1]:
            mesh_data[1].remove(y_edge[0])
    z = [(i[0], None, None, i[1], False) for i in z]


def get_unique_edges(edges):
        unique_edges = [(edge[0], edge[3]) for edge in edges]
        unique_edges = list(set(unique_edges))
        unique_edges = list({edge[0]: edge for edge in unique_edges}.values())
        unique_edges.sort(key=lambda x: x[0])
        return unique_edges 

def tranfer_box_to_polygon(box):
    start = np.fmin(box.GetStart(), box.GetStop())
    stop = np.fmax(box.GetStart(), box.GetStop())
    x_coords = [start[0], stop[0], stop[0], start[0], start[0]]
    y_coords = [start[1], start[1], stop[1], stop[1], start[1]]
    z_coords = [float(start[2]), float(stop[2])]
    return x_coords, y_coords, z_coords

def transfer_port_to_polygon(start, stop):
    port_coords_x = [start[0], stop[0], stop[0], start[0], start[0]]
    port_coords_y = [start[1], start[1], stop[1], stop[1], start[1]]
    port_coords_z = [start[2], stop[2]]
    return port_coords_x, port_coords_y, port_coords_z


def _is_polyhedron_type(prim) -> bool:
    """True for CSPrimPolyhedron / CSPrimPolyhedronReader (no GetCoords API)."""
    if not hasattr(prim, 'GetType'):
        return False
    try:
        tn = (prim.GetTypeName() or '').lower() if hasattr(prim, 'GetTypeName') else ''
        if 'polyhedron' in tn:
            return True
    except Exception:
        pass
    try:
        t = prim.GetType()
    except Exception:
        return False
    for n in ('POLYHEDRON', 'POLYHEDRONREADER'):
        v = getattr(CSPrimitives.PrimitiveType, n, None)
        if v is None:
            continue
        if t == v:
            return True
        try:
            if int(t) == int(v):
                return True
        except Exception:
            pass
    return False


def _prim_aabb_xy_z(prim):
    """Return (x_coords, y_coords, z_min, z_max) from a polyhedron AABB."""
    if not hasattr(prim, 'GetBoundBox'):
        return None
    bb = prim.GetBoundBox()
    if bb is None:
        return None
    try:
        lo, hi = bb[0], bb[1]
        xmin, ymin, zmin = float(lo[0]), float(lo[1]), float(lo[2])
        xmax, ymax, zmax = float(hi[0]), float(hi[1]), float(hi[2])
    except Exception:
        return None
    if not all(np.isfinite([xmin, ymin, zmin, xmax, ymax, zmax])):
        return None
    x_coords = [xmin, xmax, xmax, xmin, xmin]
    y_coords = [ymin, ymin, ymax, ymax, ymin]
    return x_coords, y_coords, zmin, zmax


def _handle_polyhedron(automesher, prim, x, y, x_edges, y_edges, diagonal_edges):
    data = _prim_aabb_xy_z(prim)
    if data is None:
        return
    x_coords, y_coords, _zmin, _zmax = data
    x.extend(x_coords)
    y.extend(y_coords)
    collect_edges(x_coords, y_coords, prim, x_edges, y_edges, diagonal_edges)


def process_primitive(automesher, prim, x, y, x_edges, y_edges, diagonal_edges):
    if not hasattr(prim, 'GetType'):
        port_coords_x, port_coords_y, port_coords_z = transfer_port_to_polygon(prim.start, prim.stop)
        x.extend(port_coords_x)
        y.extend(port_coords_y)
        collect_edges(port_coords_x, port_coords_y, prim, x_edges, y_edges, diagonal_edges)
        return
    ptype = prim.GetType()
    if ptype == CSPrimitives.PrimitiveType.BOX:
        box_coords_x, box_coords_y, box_coords_z = tranfer_box_to_polygon(prim)
        x.extend(box_coords_x)
        y.extend(box_coords_y)
        collect_edges(box_coords_x, box_coords_y, prim, x_edges, y_edges, diagonal_edges)
        return
    if _is_polyhedron_type(prim):
        _handle_polyhedron(automesher, prim, x, y, x_edges, y_edges, diagonal_edges)
        return
    if ptype == CSPrimitives.PrimitiveType.CYLINDER:
        start = prim.GetStart()
        stop = prim.GetStop()
        radius = prim.GetRadius()
        # Contribute BOTH end discs: a cylinder along x/y would otherwise only
        # register the start end and the mesh (and with it the FDTD domain)
        # stops inside the wire. For z-cylinders start/stop share x/y, so the
        # extra entries are duplicates and get de-duplicated downstream.
        for end in (start, stop):
            x_vals = [end[0] - radius, end[0] + radius]
            y_vals = [end[1] - radius, end[1] + radius]
            x.extend(x_vals)
            y.extend(y_vals)
            x_edges.append([x_vals[0], y_vals[0], y_vals[1], prim, False])
            x_edges.append([x_vals[1], y_vals[0], y_vals[1], prim, False])
            y_edges.append([y_vals[0], x_vals[0], x_vals[1], prim, False])
            y_edges.append([y_vals[1], x_vals[0], x_vals[1], prim, False])
        automesher.found_circles.append([[start[0]-radius,start[0]+radius],[start[1]-radius, start[1]+radius]])
        return
    if ptype == CSPrimitives.PrimitiveType.CYLINDRICALSHELL:
        start = prim.GetStart()
        stop = prim.GetStop()
        radius = prim.GetRadius()
        thickness = prim.GetShellWidth()
        outer = radius + thickness / 2
        for end in (start, stop):
            x.extend([end[0] - outer, end[0] + outer])
            y.extend([end[1] - outer, end[1] + outer])
            for rr in (radius, outer):
                x_vals = [end[0] - rr, end[0] + rr]
                y_vals = [end[1] - rr, end[1] + rr]
                x_edges.append([x_vals[0], y_vals[0], y_vals[1], prim, False])
                x_edges.append([x_vals[1], y_vals[0], y_vals[1], prim, False])
                y_edges.append([y_vals[0], x_vals[0], x_vals[1], prim, False])
                y_edges.append([y_vals[1], x_vals[0], x_vals[1], prim, False])
        automesher.found_circles.append([[start[0]-radius,start[0]+radius],[start[1]-radius, start[1]+radius]])
        automesher.found_circles.append([[start[0]-radius-thickness/2,start[0]+radius+thickness/2],[start[1]-radius-thickness/2, start[1]+radius+thickness/2]])
        return
    if hasattr(prim, 'GetCoords'):
        xx, yy = prim.GetCoords()[0], prim.GetCoords()[1]
        x.extend(xx)
        y.extend(yy)
        if xx[-1] != xx[0] or yy[-1] != yy[0]:
            xx = np.append(xx, xx[0])
            yy = np.append(yy, yy[0])
        collect_edges(xx, yy, prim, x_edges, y_edges, diagonal_edges)
    else:
        _handle_polyhedron(automesher, prim, x, y, x_edges, y_edges, diagonal_edges)


def collect_edges(x_coords, y_coords, prim, x_edges, y_edges, diagonal_edges):
    for i in range(len(x_coords) - 1):
        if x_coords[i] != x_coords[i + 1] and y_coords[i] != y_coords[i + 1]:
            diagonal_edges.append([x_coords[i], x_coords[i + 1], y_coords[i], y_coords[i + 1], prim])
        if x_coords[i] == x_coords[i + 1]:
            x_edges.append([x_coords[i], y_coords[i], y_coords[i + 1], prim, False])
        if y_coords[i] == y_coords[i + 1]:
            y_edges.append([y_coords[i], x_coords[i], x_coords[i + 1], prim, False])

def collect_z_coordinates(polygon):
    z = []
    for prim in polygon:
        if not hasattr(prim, 'GetType'):
            continue
        if _is_polyhedron_type(prim):
            data = _prim_aabb_xy_z(prim)
            if data is not None:
                z.append((data[2], prim))
                z.append((data[3], prim))
            continue
        if prim.GetType() == CSPrimitives.PrimitiveType.BOX:
            box_coords_z = tranfer_box_to_polygon(prim)[2]
            z.append((box_coords_z[0], prim))
            z.append((box_coords_z[1], prim))
        elif prim.GetType() == CSPrimitives.PrimitiveType.CYLINDER:
            z.append((prim.GetStart()[2], prim))
            z.append((prim.GetStop()[2], prim))
        elif prim.GetType() == CSPrimitives.PrimitiveType.CYLINDRICALSHELL:
            z.append((prim.GetStart()[2], prim))
            z.append((prim.GetStop()[2], prim))
        else:
            try:
                z.append((prim.GetElevation(), prim))
            except Exception:  # noqa: BLE001
                data = _prim_aabb_xy_z(prim)
                if data is not None:
                    z.append((data[2], prim))
                    z.append((data[3], prim))
            if prim.GetType() == CSPrimitives.PrimitiveType.LINPOLY:
                try:
                    z.append((prim.GetElevation() + prim.GetLength(), prim))
                except Exception:  # noqa: BLE001
                    pass
    z = list(set(z))
    z.sort(key=lambda x: x[0])
    return z

def process_z_coordinates(automesher, z, mesh_data):
    for z_val, prim in z:
        dirs = automesher.primitives_mesh_setup.get(prim, {}).get('dirs') or \
        automesher.properties_mesh_setup.get(prim.GetProperty(), {}).get('dirs') or \
        automesher.global_mesh_setup.get('dirs')
        if dirs is not None and 'z' in dirs:
            mesh_data[2].append(z_val)            

def process_single_polygon(automesher, prim, x, y, x_edges, y_edges, diagonal_edges):
    process_primitive(automesher, prim, x, y, x_edges, y_edges, diagonal_edges)


def detect_all_circles_in_polygon(automesher, polygon, min_points=20, tolerance=0.01):
    
    if automesher.global_mesh_setup.get('use_circle_detection', False) is False:
        return []
    if isinstance(polygon, list):
        coords = [prim.GetCoords() for prim in polygon if hasattr(prim, 'GetCoords')]
        if not coords:
            return []
        x_coords = np.concatenate([coord[0] for coord in coords])
        y_coords = np.concatenate([coord[1] for coord in coords])
    else:                
        coords = polygon.GetCoords() if hasattr(polygon, 'GetCoords' ) else None
        if coords is not None:
            x_coords = np.array(coords[0])
            y_coords = np.array(coords[1])
        else:
            return []
    N = len(x_coords)

    found_segments = []
    used_indices = set()

    for start in range(N - min_points + 1):
        for length in range(min_points, N - start + 1):
            indices = tuple(range(start, start + length))

            # Falls diese Punkte schon Teil eines gefundenen Kreises sind, überspringen
            if any(i in used_indices for i in indices):
                continue

            sub_x = x_coords[list(indices)]
            sub_y = y_coords[list(indices)]

            # Mittelpunkt-Näherung
            x0, y0 = np.mean(sub_x), np.mean(sub_y)
            radii = np.sqrt((sub_x - x0)**2 + (sub_y - y0)**2)
            mean_radius = np.mean(radii)
            deviation = np.abs(radii - mean_radius) / mean_radius

            if np.all(deviation < tolerance):
                found_segments.append((sub_x, sub_y))
                used_indices.update(indices)
                break  # Nicht überlappend: nächster Startpunkt
    if found_segments:
        automesher.found_circles.extend(found_segments)

def fit_circle_least_squares(x, y):
    """
    Least-squares Kreis-Fit nach algebraischer Methode.

        np.c_[np.array([1,2,3]), np.array([4,5,6])]:
            array([[1, 4],
                    [2, 5],
                    [3, 6]])
    """
    A = np.c_[2*x, 2*y, np.ones(len(x))]
    b = x**2 + y**2
    sol, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
    xc, yc = sol[0], sol[1]
    r = np.sqrt(sol[2] + xc**2 + yc**2)
    return xc, yc, r

def detect_all_arcs_in_polygon(automesher, polygon, min_points=10, tolerance=0.03, max_radius_dev=0.03):

    if automesher.global_mesh_setup.get('use_arc_detection', False) is False:
        return []
    if isinstance(polygon, list):
        coords = [prim.GetCoords() for prim in polygon if hasattr(prim, 'GetCoords')]
        if not coords:
            return []
        x_coords = np.concatenate([coord[0] for coord in coords])
        y_coords = np.concatenate([coord[1] for coord in coords])
    else:
        coords = polygon.GetCoords() if hasattr(polygon, 'GetCoords') else None
        if coords is None:
            return []
        x_coords, y_coords = np.array(coords[0]), np.array(coords[1])
    N = len(x_coords)

    found_arcs = []
    used_indices = set()

    for start in range(N - min_points + 1):
        for length in range(min_points, N - start + 1):
            indices = tuple(range(start, start + length))

            if any(i in used_indices for i in indices):
                continue

            sub_x = x_coords[list(indices)]
            sub_y = y_coords[list(indices)]

            try:
                xc, yc, r = fit_circle_least_squares(sub_x, sub_y)
            except Exception:
                continue

            distances = np.sqrt((sub_x - xc)**2 + (sub_y - yc)**2)
            mean_radius = np.mean(distances)
            deviation = np.abs(distances - mean_radius) / mean_radius

            if np.max(deviation) < max_radius_dev:
                found_arcs.append((sub_x, sub_y))
                used_indices.update(indices)
                break  # Nicht überlappen, nächster Start

    return found_arcs


def is_circle_in_polygon(polygon, min_points=8, tolerance=0.01, return_segment=False):
    if not hasattr(polygon, 'GetCoords'):
        return (False, None, None) if return_segment else False
    coords = polygon.GetCoords()
    x_coords, y_coords = np.array(coords[0]), np.array(coords[1])
    N = len(x_coords)
    for start in range(N - min_points + 1):
        for length in range(min_points, N - start + 1):
            sub_x = x_coords[start:start+length]
            sub_y = y_coords[start:start+length]
            if len(sub_x) < 3:
                continue
            x0, y0 = np.mean(sub_x), np.mean(sub_y)
            radii = np.sqrt((sub_x - x0) ** 2 + (sub_y - y0) ** 2)
            mean_radius = np.mean(radii)
            deviation = np.abs(radii - mean_radius) / mean_radius
            if np.all(deviation < tolerance):
                if return_segment:
                    return True, sub_x, sub_y
                return True
    return (False, None, None) if return_segment else False


def calc_min_distance( x):
    min_distance = float('inf')
    for i in range(len(x)):
        for j in range(i + 1, len(x)):
            distance = abs(x[i] - x[j])
            if distance > 0 and distance < min_distance:
                min_distance = distance
    return min_distance

def point_in_polygon(polygon, point):
    """
    Raycasting Algorithm to find out whether a point is in a given polygon.
    Performs the even-odd-rule Algorithm to find out whether a point is in a given polygon.
    This runs in O(n) where n is the number of edges of the polygon.
    *
    :param polygon: an array representation of the polygon where polygon[i][0] is the x Value of the i-th point and polygon[i][1] is the y Value.
    :param point:   an array representation of the point where point[0] is its x Value and point[1] is its y Value
    :return: whether the point is in the polygon (not on the edge, just turn < into <= and > into >= for that)
    """

    # A point is in a polygon if a line from the point to infinity crosses the polygon an odd number of times
    odd = False
    # For each edge (In this case for each point of the polygon and the previous one)
    i = 0
    j = len(polygon[0]) - 1
    while i < len(polygon[0]) - 1:
        i = i + 1
        # If a line from the point into infinity crosses this edge
        # One point needs to be above, one below our y coordinate
        # ...and the edge doesn't cross our Y corrdinate before our x coordinate (but between our x coordinate and infinity)

        if (((polygon[1][i] > point[1]) != (polygon[1][j] > point[1])) and (point[0] < ((polygon[0][j] - polygon[0][i]) * (point[1] - polygon[1][i]) / (polygon[1][j] - polygon[1][i])) +polygon[0][i])):                # Invert odd
            odd = not odd
        j = i
    # If the number of crossings was odd, the point is in the polygon
    return odd

def yline_in_polygon(polygon, x_point, y_start, y_end):

    loop = np.linspace(y_start, y_end, 10)
    loop = loop[1:-1] 

    for y_val in loop:
        point = [x_point, y_val]
        if not point_in_polygon(polygon, point):
            return False

    return True

def xline_in_polygon(polygon, x_start, x_end, y_point):

    loop = np.linspace(x_start, x_end, 10)
    loop = loop[1:-1] 

    for x_val in loop:
        point = [x_val, y_point]
        if not point_in_polygon(polygon, point):
            return False

    return True

def metal_edge(automesher, edges, x_coords, y_coords, mesh_data, direction):

    min_distance_x = calc_min_distance(x_coords)
    min_distance_y = calc_min_distance(y_coords)
    edges_to_add = []
    edges_to_remove = []
    mer = np.array([-1.0, 2.0]) / 3 * automesher.max_res 
    if automesher.global_mesh_setup.get('smooth_metal_edge', False) == 'one_third_two_thirds':
        # Pre-pass: densify at ALL cylinder end discs regardless of edge pairing.
        # Thin wires (e.g. Yagi elements, r=1 mm) have interleaved edges from
        # multiple directors; consecutive edge pairs have tiny gaps that fail
        # the [max_res, mesh_res] check. Handle cylinders separately.
        _seen_cyl = set()
        for edge in edges:
            prim = edge[3]
            _pid = id(prim)
            if _pid in _seen_cyl:
                continue
            if not (hasattr(prim, 'GetType') and prim.GetType() == CSPrimitives.PrimitiveType.CYLINDER):
                continue
            if not (hasattr(prim, 'GetProperty') and isinstance(prim.GetProperty(), CSProperties.CSPropMetal)):
                continue
            _seen_cyl.add(_pid)
            _radius = prim.GetRadius()
            _start = prim.GetStart()
            _stop = prim.GetStop()
            _target = automesher.cylinder_end_lines[0 if direction == 'x' else 1]
            for _end_pt in (_start, _stop):
                _ex = _end_pt[0] if direction == 'x' else _end_pt[1]
                _e1 = _ex - _radius
                _e2 = _ex + _radius
                _mid = (_e1 + _e2) / 2.0
                _third = (_e2 - _e1) / 3.0
                for _v in (_e1, _e1 + _third, _mid, _e2 - _third, _e2):
                    if _v not in _target:
                        _target.append(_v)
        for i in range(len(edges)-1):  # Iterate over a copy of the list
            edge = edges[i]
            next_edge = edges[i + 1]
            # print(f"Processing edge {i}: {edge} with next edge {next_edge}")
            if hasattr(edge[3], 'GetProperty') and isinstance(edge[3].GetProperty(), CSProperties.CSPropMetal) or \
                (hasattr(edge[3],'priority') and isinstance(edge[3],openEMS.ports.MSLPort)):
                if (hasattr(edge[3],'priority') and not isinstance(edge[3],openEMS.ports.MSLPort)):
                    continue
                if edge[0] == next_edge[0] :
                    continue
                if edge[0] == min(edges, key=lambda e: e[0])[0]:
                    continue
                edges_with_same_x = [e for e in edges if e[0] == edge[0] and (hasattr(e[3], 'priority') and isinstance(e[3], openEMS.ports.MSLPort)) and not e[3] == edge[3]]
                next_edges_with_same_x = [e for e in edges if e[0] == next_edge[0] and (hasattr(e[3], 'priority') and isinstance(e[3], openEMS.ports.MSLPort)) and not e[3] == next_edge[3]]
                x, y, x_edges, y_edges, diagonal_edges = [], [], [], [], []
                process_primitive(automesher, edge[3], x, y, x_edges, y_edges, diagonal_edges)
                coords = [x, y]
                _gap = abs(np.diff([edge[0], next_edge[0]]))
                _min_gap = min(automesher.max_res, 1.5)
                if _gap <= automesher.mesh_res and _gap >= _min_gap:
                        y1, y2 = edge[1], edge[2]
                        y1_next, y2_next = next_edge[1], next_edge[2]
                        if (y1 <= y1_next <= y2 or y1 >= y1_next >= y2 or
                            y1 <= y2_next <= y2 or y1 >= y2_next >= y2 or
                            y1_next <= y1 <= y2_next or y1_next >= y1 >= y2_next or
                            y1_next <= y2 <= y2_next or y1_next >= y2 >= y2_next):
                            new_max_res = np.max(np.diff(np.linspace(edge[0], next_edge[0], automesher.num_lines)))
                            mer = np.array([-new_max_res/2, new_max_res])

                            # print(f"New max resolution: {new_max_res}")
                            # xlines = SmoothMeshLines([edge[0], next_edge[0]], new_max_res)
                            # print(f"Adding xlines: {xlines}")
                            # mesh_data.extend(xlines)
                            if direction == 'x':
                                condition1 = yline_in_polygon(coords, edge[0]+min_distance_x/2, edge[1], edge[2]) and not yline_in_polygon(coords, edge[0]-min_distance_x/2, edge[1], edge[2])
                                condition2 = yline_in_polygon(coords, edge[0]-min_distance_x/2, edge[1], edge[2]) and yline_in_polygon(coords, edge[0]+min_distance_x/2, edge[1], edge[2])
                                condition3 = yline_in_polygon(coords, next_edge[0]+min_distance_x/2, next_edge[1], next_edge[2]) and not yline_in_polygon(coords, next_edge[0]-min_distance_x/2, next_edge[1], next_edge[2])
                                condition4 = yline_in_polygon(coords, next_edge[0]-min_distance_x/2, next_edge[1], next_edge[2]) and yline_in_polygon(coords, next_edge[0]+min_distance_x/2, next_edge[1], next_edge[2])
                            if direction == 'y':
                                condition1 = xline_in_polygon(coords, edge[1], edge[2], edge[0]+min_distance_y/2) and not xline_in_polygon(coords, edge[1], edge[2], edge[0]-min_distance_y/2)
                                condition2 = xline_in_polygon(coords, edge[1], edge[2], edge[0]-min_distance_y/2) and xline_in_polygon(coords, edge[1], edge[2], edge[0]+min_distance_y/2)
                                condition3 = xline_in_polygon(coords, next_edge[1], next_edge[2], next_edge[0]+min_distance_y/2) and not xline_in_polygon(coords, next_edge[1], next_edge[2], next_edge[0]-min_distance_y/2)
                                condition4 = xline_in_polygon(coords, next_edge[1], next_edge[2], next_edge[0]-min_distance_y/2) and xline_in_polygon(coords, next_edge[1], next_edge[2], next_edge[0]+min_distance_y/2)
                            if condition1:
                                # mesh_data_in_range = [mesh_data for mesh_data in mesh_data if edge[0]-mer[1] <= mesh_data <= edge[0]-mer[0]]
                                # if not mesh_data_in_range:
                                #     edges_to_add.append([edge[0]-mer[1], edge[1], edge[2], edge[3], True])
                                #     edges_to_add.append([edge[0]-mer[0], edge[1], edge[2], edge[3], True])
                                #     edges_to_remove.append(edge)
                                # else:
                                    # mesh_data = [h for h in mesh_data if h not in mesh_data_in_range]
                                    if edges_with_same_x:
                                        continue
                                    mesh_data.append(edge[0]-mer[1])
                                    mesh_data.append(edge[0]-mer[0])
                                    edges_to_add.append([edge[0]-mer[1], edge[1], edge[2], edge[3], True])
                                    edges_to_add.append([edge[0]-mer[0], edge[1], edge[2], edge[3], True])
                                    edges_to_add.append([edge[0], edge[1], edge[2], edge[3], True])
                                    edges_to_remove.append(edge)
                                    edges[i][4] = True
                            elif condition2:
                                continue
                            else:
                                # mesh_data_in_range = [mesh_data for mesh_data in mesh_data if edge[0]+mer[0] <= mesh_data <= edge[0]+mer[1]]
                                # if not mesh_data_in_range:
                                #     edges_to_add.append([edge[0]+mer[0], edge[1], edge[2], edge[3], True])
                                #     edges_to_add.append([edge[0]+mer[1], edge[1], edge[2], edge[3], True])
                                #     edges_to_remove.append(edge)
                                # else:
                                    # mesh_data = [h for h in mesh_data if h not in mesh_data_in_range]
                                    if edges_with_same_x:
                                        continue
                                    mesh_data.append(edge[0]+mer[0])
                                    mesh_data.append(edge[0]+mer[1])
                                    edges_to_add.append([edge[0]+mer[0], edge[1], edge[2], edge[3], True])
                                    edges_to_add.append([edge[0]+mer[1], edge[1], edge[2], edge[3], True])
                                    edges_to_add.append([edge[0], edge[1], edge[2], edge[3], True])
                                    edges_to_remove.append(edge)
                                    edges[i][4] = True
                            if condition3:
                                # mesh_data_in_range = [mesh_data for mesh_data in mesh_data if next_edge[0]-mer[1] <= mesh_data <= next_edge[0]-mer[0]]
                                # if not mesh_data_in_range:
                                #     edges_to_add.append([next_edge[0]-mer[1], next_edge[1], next_edge[2], next_edge[3], True])
                                #     edges_to_add.append([next_edge[0]-mer[0], next_edge[1], next_edge[2], next_edge[3], True])
                                #     edges_to_remove.append(next_edge)
                                # else:
                                    # mesh_data = [h for h in mesh_data if h not in mesh_data_in_range]
                                    if next_edges_with_same_x:
                                        continue
                                    mesh_data.append(next_edge[0]-mer[1])
                                    mesh_data.append(next_edge[0]-mer[0])
                                    edges_to_add.append([next_edge[0]-mer[1], next_edge[1], next_edge[2], next_edge[3], True])
                                    edges_to_add.append([next_edge[0]-mer[0], next_edge[1], next_edge[2], next_edge[3], True])
                                    edges_to_add.append([next_edge[0], next_edge[1], next_edge[2], next_edge[3], True])
                                    edges_to_remove.append(next_edge)
                                    edges[i + 1][4] = True
                            elif condition4:
                                continue
                            else:
                                # mesh_data_in_range = [mesh_data for mesh_data in mesh_data if next_edge[0]+mer[0] <= mesh_data <= next_edge[0]+mer[1]]
                                # if not mesh_data_in_range:
                                #     edges_to_add.append([next_edge[0]+mer[0], next_edge[1], next_edge[2], next_edge[3], True])
                                #     edges_to_add.append([next_edge[0]+mer[1], next_edge[1], next_edge[2], next_edge[3], True])
                                #     edges_to_remove.append(next_edge)
                                # else:
                                    # mesh_data = [h for h in mesh_data if h not in mesh_data_in_range]
                                    if next_edges_with_same_x:
                                        continue
                                    mesh_data.append(next_edge[0]+mer[0])
                                    mesh_data.append(next_edge[0]+mer[1])
                                    edges_to_add.append([next_edge[0]+mer[0], next_edge[1], next_edge[2], next_edge[3], True])
                                    edges_to_add.append([next_edge[0]+mer[1], next_edge[1], next_edge[2], next_edge[3], True])
                                    edges_to_add.append([next_edge[0], next_edge[1], next_edge[2], next_edge[3], True])
                                    edges_to_remove.append(next_edge)
                                    edges[i + 1][4] = True

                if direction == 'x':
                    condition1 = yline_in_polygon(coords, edge[0]+min_distance_x/2, edge[1], edge[2]) and not yline_in_polygon(coords, edge[0]-min_distance_x/2, edge[1], edge[2])
                    condition2 = yline_in_polygon(coords, edge[0]-min_distance_x/2, edge[1], edge[2]) and yline_in_polygon(coords, edge[0]+min_distance_x/2, edge[1], edge[2])
                if direction == 'y':
                    condition1 = xline_in_polygon(coords, edge[1], edge[2], edge[0]+min_distance_y/2) and not xline_in_polygon(coords, edge[1], edge[2], edge[0]-min_distance_y/2)
                    condition2 = xline_in_polygon(coords, edge[1], edge[2], edge[0]-min_distance_y/2) and xline_in_polygon(coords, edge[1], edge[2], edge[0]+min_distance_y/2)
                if condition1:
                    # mesh_data_in_range = [mesh_data for mesh_data in mesh_data if edge[0]-mer[1] <= mesh_data <= edge[0]-mer[0]]
                    # if not mesh_data_in_range:
                    #     edges_to_add.append([edge[0]-mer[1], edge[1], edge[2], edge[3], True])
                    #     edges_to_add.append([edge[0]-mer[0], edge[1], edge[2], edge[3], True])
                    #     edges_to_remove.append(edge)
                    # else:
                        # mesh_data = [h for h in mesh_data if h not in mesh_data_in_range]
                        mesh_data.append(edge[0]-mer[1])
                        mesh_data.append(edge[0]-mer[0])
                        edges_to_add.append([edge[0]-mer[1], edge[1], edge[2], edge[3], True])
                        edges_to_add.append([edge[0]-mer[0], edge[1], edge[2], edge[3], True])
                        edges_to_add.append([edge[0], edge[1], edge[2], edge[3], True])
                        edges_to_remove.append(edge)
                        edges[i][4] = True
                elif condition2:
                    continue
                else:
                    # mesh_data_in_range = [mesh_data for mesh_data in mesh_data if edge[0]+mer[0] <= mesh_data <= edge[0]+mer[1]]
                    # if not mesh_data_in_range:
                    #     edges_to_add.append([edge[0]+mer[0], edge[1], edge[2], edge[3], True])
                    #     edges_to_add.append([edge[0]+mer[1], edge[1], edge[2], edge[3], True])
                    #     edges_to_remove.append(edge)
                    # else:
                        # mesh_data = [h for h in mesh_data if h not in mesh_data_in_range]
                        mesh_data.append(edge[0]+mer[0])
                        mesh_data.append(edge[0]+mer[1])
                        edges_to_add.append([edge[0]+mer[0], edge[1], edge[2], edge[3], True])
                        edges_to_add.append([edge[0]+mer[1], edge[1], edge[2], edge[3], True])
                        edges_to_add.append([edge[0], edge[1], edge[2], edge[3], True])
                        edges_to_remove.append(edge)
                        edges[i][4] = True
        # edges.extend(edges_to_add)
        # for edge in edges_to_remove:
        #     if edge in edges:  
        #         edges.remove(edge)
        # if direction == 'x':
        #     min_distance_x = calc_min_distance(x)
        # if direction == 'y':
        #     min_distance_y = calc_min_distance(y)
        # for i in range(len(edges) - 1):
        #     if direction == 'x':
        #         condition1 = yline_in_polygon(coords, edges[i][0]+min_distance_x/2, edges[i][1], edges[i][2]) and not yline_in_polygon(coords, edges[i][0]-min_distance_x/2, edges[i][1], edges[i][2])
        #         condition2 = yline_in_polygon(coords, edges[i][0]-min_distance_x/2, edges[i][1], edges[i][2]) and yline_in_polygon(coords, edges[i][0]+min_distance_x/2, edges[i][1], edges[i][2])
        #     if direction == 'y':
        #         condition1 = xline_in_polygon(coords, edges[i][1], edges[i][2], edges[i][0]+min_distance_y/2) and not xline_in_polygon(coords, edges[i][1], edges[i][2], edges[i][0]-min_distance_y/2)
        #         condition2 = xline_in_polygon(coords, edges[i][1], edges[i][2], edges[i][0]-min_distance_y/2) and xline_in_polygon(coords, edges[i][1], edges[i][2], edges[i][0]+min_distance_y/2)
        #         if i > 0 and abs(edges[i][0] - edges[i + 1][0]) > mesh_res and abs(edges[i][0] - edges[i - 1][0]) > mesh_res:
        #             if condition1:
        #                 mesh_data_in_range =  [mesh_data for mesh_data in mesh_data if edges[i][0]-mer[1] <= mesh_data <= edges[i][0]-mer[0]]
        #                 if not mesh_data_in_range:
        #                     mesh_data.append(edges[i][0]-mer[1])
        #                     mesh_data.append(edges[i][0]-mer[0])
        #                 else:
        #                     mesh_data = [h for h in mesh_data if h not in mesh_data_in_range]
        #                     mesh_data.append(edges[i][0]-mer[1])
        #                     mesh_data.append(edges[i][0]-mer[0])
        #             elif condition2:
        #                 continue
        #             else:
        #                 mesh_data_in_range =  [mesh_data for mesh_data in mesh_data if edges[i][0]+mer[0] <= mesh_data <= edges[i][0]+mer[1]]
        #                 if not mesh_data_in_range:
        #                     mesh_data.append(edges[i][0]+mer[0])
        #                     mesh_data.append(edges[i][0]+mer[1])
        #                 else:
        #                     mesh_data = [h for h in mesh_data if h not in mesh_data_in_range]
        #                     mesh_data.append(edges[i][0]+mer[0])
        #                     mesh_data.append(edges[i][0]+mer[1])
    elif automesher.global_mesh_setup.get('smooth_metal_edge', False) == 'extra_lines':
        for edge in edges:  # Iterate over a copy of the list
            if hasattr(edge[3], 'GetProperty') and isinstance(edge[3].GetProperty(), CSProperties.CSPropMetal) or \
                (hasattr(edge[3],'priority') and isinstance(edge[3],openEMS.ports.MSLPort)):
                if edge[0] == min(edges, key=lambda e: e[0])[0]:
                    continue  
                if edge[0] == max(edges, key=lambda e: e[0])[0]:
                    continue 
                mesh_data.append(edge[0] - automesher.min_cellsize)
                mesh_data.append(edge[0])
                mesh_data.append(edge[0] + automesher.min_cellsize)
