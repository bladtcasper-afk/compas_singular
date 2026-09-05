"""**A picture of the current stage, with no plotting library.**

Neither interpreter this pipeline runs in has matplotlib or Pillow -- not the
``singular312`` conda env and not Rhino 8's CPython, both checked. Rather than
add a dependency that would have to be installed in Rhino's site-env before a
command could draw anything, this module writes the PNG itself: ``zlib`` from
the standard library, a scanline buffer, and Bresenham. It therefore runs
anywhere the rest of the library runs, which is the whole point.

:func:`render_png` is for a model that reads images; :func:`render_svg` is the
same scene as text, which is smaller, sharper and diffable, and is what a human
or an HTML page should get. **A vision API takes PNG or JPEG and not SVG**, so
the two are not interchangeable and both exist.

**What is drawn depends on the stage**, because the useful picture is different
at each: the field stage shows the domain and its separatrices, the coarse stage
shows patches large enough to label individually, and the dense stage shows a
mesh with far too many faces to label and is drawn for its singularities and its
flow instead.

The glyphs are a 3x5 bitmap font covering the digits and ``v``/``e``/``f``,
which is exactly the alphabet :class:`AddressBook` labels use. It is scaled up
rather than smoothed; at scale 2 or 3 it is legible and it costs sixteen
integers.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import struct
import zlib


__all__ = ['Canvas', 'View', 'render_png', 'render_svg', 'scene_for']


WHITE = (255, 255, 255)
BLACK = (0, 0, 0)
WALL = (20, 20, 20)
PATCH_EDGE = (70, 110, 170)
DENSE_EDGE = (150, 165, 185)
SEPARATRIX = (240, 150, 40)
POSITIVE = (200, 40, 40)      # a quad mesh wants valence 3 here
NEGATIVE = (40, 80, 200)      # ... and valence 5 here
POLE = (200, 0, 120)
GUIDE = (150, 60, 190)
LABEL = (90, 90, 90)


#: 3x5, one int per row, bit 2 leftmost. The alphabet AddressBook labels use.
FONT = {
    '0': (0b111, 0b101, 0b101, 0b101, 0b111),
    '1': (0b010, 0b110, 0b010, 0b010, 0b111),
    '2': (0b111, 0b001, 0b111, 0b100, 0b111),
    '3': (0b111, 0b001, 0b111, 0b001, 0b111),
    '4': (0b101, 0b101, 0b111, 0b001, 0b001),
    '5': (0b111, 0b100, 0b111, 0b001, 0b111),
    '6': (0b111, 0b100, 0b111, 0b101, 0b111),
    '7': (0b111, 0b001, 0b001, 0b001, 0b001),
    '8': (0b111, 0b101, 0b111, 0b101, 0b111),
    '9': (0b111, 0b101, 0b111, 0b001, 0b111),
    'v': (0b000, 0b101, 0b101, 0b101, 0b010),
    'e': (0b010, 0b101, 0b111, 0b100, 0b011),
    'f': (0b011, 0b100, 0b110, 0b100, 0b100),
    's': (0b011, 0b100, 0b010, 0b001, 0b110),
    # the rest of the alphabet the stage names need: field, coarse, dense
    'a': (0b000, 0b011, 0b101, 0b101, 0b011),
    'c': (0b000, 0b011, 0b100, 0b100, 0b011),
    'd': (0b001, 0b011, 0b101, 0b101, 0b011),
    'i': (0b010, 0b000, 0b010, 0b010, 0b010),
    'l': (0b110, 0b010, 0b010, 0b010, 0b111),
    'n': (0b000, 0b110, 0b101, 0b101, 0b101),
    'o': (0b000, 0b010, 0b101, 0b101, 0b010),
    'r': (0b000, 0b011, 0b100, 0b100, 0b100),
    '-': (0b000, 0b000, 0b111, 0b000, 0b000),
    '.': (0b000, 0b000, 0b000, 0b000, 0b010),
    ' ': (0b000, 0b000, 0b000, 0b000, 0b000),
}


class View(object):
    """World XY to pixels, aspect preserved, Y flipped, with a margin."""

    def __init__(self, bbox, width, height, margin=24):
        (x0, y0), (x1, y1) = bbox
        dx = max(x1 - x0, 1e-9)
        dy = max(y1 - y0, 1e-9)
        self.scale = min((width - 2.0 * margin) / dx,
                         (height - 2.0 * margin) / dy)
        # centre the drawing in whatever room the aspect ratio left over
        self.ox = (width - self.scale * dx) / 2.0 - self.scale * x0
        self.oy = (height - self.scale * dy) / 2.0 - self.scale * y0
        self.height = height

    def __call__(self, xyz):
        px = self.ox + self.scale * xyz[0]
        py = self.oy + self.scale * xyz[1]
        return int(round(px)), int(round(self.height - py))


class Canvas(object):
    """An RGB pixel buffer that can write itself as a PNG."""

    def __init__(self, width=900, height=900, background=WHITE):
        self.width = int(width)
        self.height = int(height)
        r, g, b = background
        self.buffer = bytearray(bytes((r, g, b)) * (self.width * self.height))

    # ------------------------------------------------------------------

    def pixel(self, x, y, color):
        if 0 <= x < self.width and 0 <= y < self.height:
            i = (y * self.width + x) * 3
            self.buffer[i] = color[0]
            self.buffer[i + 1] = color[1]
            self.buffer[i + 2] = color[2]

    def disc(self, centre, radius, color):
        cx, cy = centre
        r = int(radius)
        rr = r * r
        for dy in range(-r, r + 1):
            for dx in range(-r, r + 1):
                if dx * dx + dy * dy <= rr:
                    self.pixel(cx + dx, cy + dy, color)

    def line(self, a, b, color, width=1):
        """Bresenham, thickened by stamping a disc when ``width`` > 1."""
        x0, y0 = a
        x1, y1 = b
        dx = abs(x1 - x0)
        dy = -abs(y1 - y0)
        sx = 1 if x0 < x1 else -1
        sy = 1 if y0 < y1 else -1
        err = dx + dy
        radius = (width - 1) // 2
        while True:
            if radius > 0:
                self.disc((x0, y0), radius, color)
            else:
                self.pixel(x0, y0, color)
            if x0 == x1 and y0 == y1:
                break
            e2 = 2 * err
            if e2 >= dy:
                err += dy
                x0 += sx
            if e2 <= dx:
                err += dx
                y0 += sy

    def polyline(self, points, color, width=1, close=False):
        pts = list(points)
        if close and len(pts) > 2:
            pts = pts + pts[:1]
        for a, b in zip(pts, pts[1:]):
            self.line(a, b, color, width)

    def text(self, origin, string, color=LABEL, scale=2):
        """Draw ``string`` with its top-left at ``origin``. Unknown chars skipped."""
        x, y = origin
        for char in str(string):
            glyph = FONT.get(char)
            if glyph is None:
                x += 4 * scale
                continue
            for row, bits in enumerate(glyph):
                for col in range(3):
                    if bits & (1 << (2 - col)):
                        for sy in range(scale):
                            for sx in range(scale):
                                self.pixel(x + col * scale + sx,
                                           y + row * scale + sy, color)
            x += 4 * scale

    # ------------------------------------------------------------------

    def png(self):
        """The buffer as PNG bytes. Filter 0 on every scanline, then zlib."""
        raw = bytearray()
        stride = self.width * 3
        for y in range(self.height):
            raw.append(0)                       # filter type: none
            raw.extend(self.buffer[y * stride:(y + 1) * stride])

        def chunk(tag, data):
            out = struct.pack('>I', len(data)) + tag + data
            return out + struct.pack('>I', zlib.crc32(tag + data) & 0xffffffff)

        header = struct.pack('>IIBBBBB', self.width, self.height, 8, 2, 0, 0, 0)
        return (b'\x89PNG\r\n\x1a\n'
                + chunk(b'IHDR', header)
                + chunk(b'IDAT', zlib.compress(bytes(raw), 6))
                + chunk(b'IEND', b''))


# ======================================================================
# what to draw
# ======================================================================

def _bbox(groups):
    xs, ys = [], []
    for points in groups:
        for p in points:
            xs.append(p[0])
            ys.append(p[1])
    if not xs:
        return ((0.0, 0.0), (1.0, 1.0))
    return ((min(xs), min(ys)), (max(xs), max(ys)))


def _mesh_edges(mesh):
    return [[mesh.vertex_coordinates(u), mesh.vertex_coordinates(v)]
            for u, v in mesh.edges()]


def _singular_vertices(mesh):
    """``(point, index)`` for every irregular INTERIOR vertex.

    ``index`` is ``4 - valence``: positive where a quad mesh wants a valence-3
    vertex, negative where it wants valence-5. Same sign convention, and the
    same two colours, as ``framefield.viz``.
    """
    out = []
    for vkey in mesh.vertices():
        if mesh.is_vertex_on_boundary(vkey):
            continue
        valence = len(mesh.vertex_neighbors(vkey))
        if valence != 4:
            out.append((mesh.vertex_coordinates(vkey), 4 - valence))
    return out


def scene_for(session, labels=None):
    """The drawable content of the session's current stage.

    Returns
    -------
    dict
        ``walls``, ``guides``, ``separatrices`` -- lists of point lists;
        ``edges`` -- the mesh, as segments; ``singularities`` --
        ``(point, index)``; ``poles`` -- points; ``labels`` --
        ``(point, text)``; ``stage`` and ``title``.
    """
    d = session.decomposition
    stage = session.stage
    mesh = session.mesh

    scene = {
        'stage': stage,
        'walls': [list(loop) for loop in (d._loops() or [])],
        'guides': [list(g) for g in (session.guides or [])],
        'separatrices': [],
        'edges': [],
        'singularities': [],
        'poles': [],
        'labels': [],
    }

    if stage == 'field':
        scene['separatrices'] = [list(p) for p in (d.polylines or [])]
        scene['title'] = 'field'
        return scene

    scene['edges'] = _mesh_edges(mesh)
    scene['singularities'] = _singular_vertices(mesh)

    face_pole = mesh.attributes.get('face_pole') or {}
    for fkey, vkey in face_pole.items():
        if fkey in list(mesh.faces()):
            try:
                scene['poles'].append(mesh.vertex_coordinates(vkey))
            except KeyError:
                pass

    if labels:
        book = session.book()
        if book is not None:
            for label in labels:
                address = book.address(label)
                if address is not None:
                    scene['labels'].append((list(address), label))

    scene['title'] = stage
    return scene


def render_png(session, path=None, width=900, height=900, labels=None,
               scale=2):
    """**A picture of the current stage as PNG bytes**, written to ``path`` if given.

    Parameters
    ----------
    session : MeshEditSession
    path : str, optional
        Where to write. The bytes are returned either way.
    labels : list[str], optional
        Labels to draw, from :class:`AddressBook`. Pass a handful; drawing a
        thousand makes an unreadable picture and a large file.
    """
    scene = scene_for(session, labels=labels)
    groups = ([[p for seg in scene['edges'] for p in seg]]
              + scene['walls'] + scene['guides'] + scene['separatrices'])
    canvas = Canvas(width, height)
    view = View(_bbox(groups), width, height)

    for loop in scene['walls']:
        canvas.polyline([view(p) for p in loop], WALL, width=3, close=True)
    for curve in scene['separatrices']:
        canvas.polyline([view(p) for p in curve], SEPARATRIX, width=2)
    for guide in scene['guides']:
        canvas.polyline([view(p) for p in guide], GUIDE, width=3)

    edge_color = PATCH_EDGE if scene['stage'] == 'coarse' else DENSE_EDGE
    edge_width = 2 if scene['stage'] == 'coarse' else 1
    for a, b in scene['edges']:
        canvas.line(view(a), view(b), edge_color, edge_width)

    for point, index in scene['singularities']:
        canvas.disc(view(point), 5, POSITIVE if index > 0 else NEGATIVE)
    for point in scene['poles']:
        canvas.disc(view(point), 5, POLE)
    for point, text in scene['labels']:
        px, py = view(point)
        canvas.text((px + 4, py - 3 * scale), text, LABEL, scale)

    canvas.text((10, 10), scene['title'], BLACK, 2)

    data = canvas.png()
    if path:
        with open(path, 'wb') as handle:
            handle.write(data)
    return data


def render_svg(session, path=None, width=900, height=900, labels=None):
    """The same scene as SVG text. Smaller and sharper -- but **not** for a
    vision API, which takes PNG or JPEG only."""
    scene = scene_for(session, labels=labels)
    groups = ([[p for seg in scene['edges'] for p in seg]]
              + scene['walls'] + scene['guides'] + scene['separatrices'])
    view = View(_bbox(groups), width, height)

    def rgb(color):
        return 'rgb({},{},{})'.format(*color)

    parts = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {} {}" '
             'width="{}" height="{}">'.format(width, height, width, height),
             '<rect width="100%" height="100%" fill="white"/>']

    def path_of(points, close=False):
        pts = [view(p) for p in points]
        d = 'M {} {}'.format(*pts[0]) + ''.join(
            ' L {} {}'.format(x, y) for x, y in pts[1:])
        return d + (' Z' if close else '')

    edge_color = PATCH_EDGE if scene['stage'] == 'coarse' else DENSE_EDGE
    edge_width = 1.6 if scene['stage'] == 'coarse' else 0.7
    for a, b in scene['edges']:
        (x0, y0), (x1, y1) = view(a), view(b)
        parts.append('<line x1="{}" y1="{}" x2="{}" y2="{}" stroke="{}" '
                     'stroke-width="{}"/>'.format(x0, y0, x1, y1,
                                                  rgb(edge_color), edge_width))
    for curve in scene['separatrices']:
        if len(curve) > 1:
            parts.append('<path d="{}" fill="none" stroke="{}" stroke-width="1.6"/>'
                         .format(path_of(curve), rgb(SEPARATRIX)))
    for loop in scene['walls']:
        if len(loop) > 1:
            parts.append('<path d="{}" fill="none" stroke="{}" stroke-width="2.4"/>'
                         .format(path_of(loop, close=True), rgb(WALL)))
    for guide in scene['guides']:
        if len(guide) > 1:
            parts.append('<path d="{}" fill="none" stroke="{}" stroke-width="2.4"/>'
                         .format(path_of(guide), rgb(GUIDE)))
    for point, index in scene['singularities']:
        x, y = view(point)
        parts.append('<circle cx="{}" cy="{}" r="4.5" fill="{}"/>'.format(
            x, y, rgb(POSITIVE if index > 0 else NEGATIVE)))
    for point in scene['poles']:
        x, y = view(point)
        parts.append('<circle cx="{}" cy="{}" r="4.5" fill="{}"/>'.format(
            x, y, rgb(POLE)))
    for point, text in scene['labels']:
        x, y = view(point)
        parts.append('<text x="{}" y="{}" font-size="11" fill="{}">{}</text>'
                     .format(x + 5, y - 5, rgb(LABEL), text))

    parts.append('</svg>')
    out = '\n'.join(parts)
    if path:
        with open(path, 'w') as handle:
            handle.write(out)
    return out
