"""**Naming things that renumber. The layer's foundation.**

A person picks a corner with a mouse. Anything else has to NAME one, and the
three key spaces this pipeline runs on are all unstable:

* **vertex keys** are renumbered by ``edit_coarse`` -- it welds and repairs, and
  2 of 8 keys were measured to survive on a disc;
* **strip keys** and **polyedge keys** come from popping ``self.edges()``, so
  they follow vertex insertion order, and every topological edit renumbers them.

So nothing here names anything by key. A thing is named by WHERE IT IS, rounded
to :data:`PRECISION` -- the resolution ``geometric_key`` and ``from_polylines``'
endpoint matching already use, so an address made here still matches after a
round trip through either.

**Labels are minted per snapshot and are not promised to survive an edit.**
:class:`AddressBook` hands out short readable names -- ``v3``, ``e17``, ``f2`` --
because a caller reading a state digest cannot work with a tuple of floats. They
are assigned in sorted order of the rounded coordinate, so **the same geometry
always produces the same labels**: a session that changes nothing sees stable
names, and two runs of the same script address the same corners. What is NOT
promised is that ``v3`` means the same corner after a cut, a strip deletion or a
commit -- it may not exist at all. That honesty is the point:
:meth:`AddressBook.carry` re-resolves the addresses onto a rebuilt mesh and
reports what it lost, rather than silently pointing an old name at a new corner.

**A strip is addressed by its own smallest edge.** A strip has many edges and no
position of its own, so its address is the minimum edge address over the strip
-- deterministic, and independent of which edge a caller happened to pick.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function


__all__ = [
    'PRECISION', 'point_address', 'vertex_address', 'edge_address',
    'face_address', 'strip_address', 'AddressBook',
]


#: Decimals an address is rounded at. The resolution ``geometric_key`` and
#: ``CoarseEditor.PRECISION`` already use.
PRECISION = 3


def point_address(xyz, precision=PRECISION):
    """A point as a hashable, comparable address.

    A tuple of rounded floats rather than a string: it sorts numerically, which
    is what makes :class:`AddressBook`'s labels deterministic, and it round-trips
    back to coordinates without parsing.
    """
    point = list(xyz)
    while len(point) < 3:
        point.append(0.0)
    # ``round(x, 3) + 0.0`` normalises -0.0 to 0.0, so a point on an axis gets
    # the same address whichever side it was computed from.
    return tuple(round(float(v), precision) + 0.0 for v in point[:3])


def vertex_address(mesh, vkey, precision=PRECISION):
    """Where this vertex is."""
    return point_address(mesh.vertex_coordinates(vkey), precision)


def edge_address(mesh, edge, precision=PRECISION):
    """The rounded MIDPOINT of an edge.

    The midpoint rather than the pair of endpoints: it is one point, so edges
    sort the same way vertices and faces do, and ``(u, v)`` and ``(v, u)`` give
    the same answer without having to canonicalise the order.
    """
    u, v = edge
    a = mesh.vertex_coordinates(u)
    b = mesh.vertex_coordinates(v)
    return point_address([(a[i] + b[i]) / 2.0 for i in range(3)], precision)


def face_address(mesh, fkey, precision=PRECISION):
    """The rounded centroid of a face."""
    points = [mesh.vertex_coordinates(v) for v in mesh.face_vertices(fkey)]
    if not points:
        return point_address([0.0, 0.0, 0.0], precision)
    n = float(len(points))
    return point_address([sum(p[i] for p in points) / n for i in range(3)],
                         precision)


def strip_address(mesh, skey, precision=PRECISION):
    """The smallest edge address over the strip. ``None`` if it has no edges.

    Deterministic and independent of which edge a caller picked, so a density
    set through one edge is found again through any other edge of the same
    strip.
    """
    addresses = []
    for u, v in mesh.strip_edges(skey):
        if u == v:
            # A collapsed edge -- a pole's. It has no midpoint worth having and
            # would tie with the vertex itself.
            continue
        addresses.append(edge_address(mesh, (u, v), precision))
    return min(addresses) if addresses else None


class AddressBook(object):
    """Short labels for the vertices, edges and faces of one mesh.

    Built for a mesh as it stands. Labels are assigned in sorted address order,
    so identical geometry always yields identical labels.

    Parameters
    ----------
    mesh : Mesh
    precision : int, optional
    prefixes : tuple, optional
        The three label prefixes, ``('v', 'e', 'f')`` by default.

    Attributes
    ----------
    vertices, edges, faces : dict
        ``label -> address``.
    """

    def __init__(self, mesh, precision=PRECISION, prefixes=('v', 'e', 'f')):
        self.precision = precision
        self.prefixes = tuple(prefixes)
        self.vertices = {}
        self.edges = {}
        self.faces = {}
        self._vertex_keys = {}
        self._edge_keys = {}
        self._face_keys = {}
        self._build(mesh)

    # ------------------------------------------------------------------
    # building
    # ------------------------------------------------------------------

    def _build(self, mesh):
        vp, ep, fp = self.prefixes

        items = sorted((vertex_address(mesh, v, self.precision), v)
                       for v in mesh.vertices())
        for i, (address, vkey) in enumerate(items):
            label = '{}{}'.format(vp, i)
            self.vertices[label] = address
            self._vertex_keys[label] = vkey

        items = sorted((edge_address(mesh, e, self.precision), e)
                       for e in mesh.edges())
        for i, (address, edge) in enumerate(items):
            label = '{}{}'.format(ep, i)
            self.edges[label] = address
            self._edge_keys[label] = tuple(edge)

        items = sorted((face_address(mesh, f, self.precision), f)
                       for f in mesh.faces())
        for i, (address, fkey) in enumerate(items):
            label = '{}{}'.format(fp, i)
            self.faces[label] = address
            self._face_keys[label] = fkey

    # ------------------------------------------------------------------
    # label -> live key, on the mesh this book was built from
    # ------------------------------------------------------------------

    def vertex(self, label):
        """The vertex key for a label, or ``None``."""
        return self._vertex_keys.get(label)

    def edge(self, label):
        """The ``(u, v)`` for a label, or ``None``."""
        return self._edge_keys.get(label)

    def face(self, label):
        """The face key for a label, or ``None``."""
        return self._face_keys.get(label)

    def address(self, label):
        """The geometric address behind a label. This is the part that travels."""
        for table in (self.vertices, self.edges, self.faces):
            if label in table:
                return table[label]
        return None

    def label_of_vertex(self, vkey):
        for label, key in self._vertex_keys.items():
            if key == vkey:
                return label
        return None

    def label_of_edge(self, edge):
        u, v = edge
        for label, key in self._edge_keys.items():
            if key == (u, v) or key == (v, u):
                return label
        return None

    def label_of_face(self, fkey):
        for label, key in self._face_keys.items():
            if key == fkey:
                return label
        return None

    # ------------------------------------------------------------------
    # carrying labels across a rebuild
    # ------------------------------------------------------------------

    def carry(self, mesh, tol=None):
        """Re-resolve these labels onto ``mesh``. ``(new book, report)``.

        The new book is built for ``mesh`` in the normal way -- its labels are
        its own. The report says what happened to THIS book's labels, which is
        what a caller holding a plan across a commit needs:

        ``kept``
            ``{old label: new label}`` for every address still present.
        ``lost``
            Old labels whose address is gone.
        ``moved``
            Old labels whose address survived but under a different label,
            i.e. the geometry is the same and the numbering changed.

        Matching is on the rounded address, with an optional ``tol`` for a
        near-match when a repair nudged a corner. Without ``tol`` it is exact,
        which is what a pure renumbering needs.
        """
        book = AddressBook(mesh, precision=self.precision, prefixes=self.prefixes)
        report = {'kept': {}, 'lost': [], 'moved': []}

        for kind in ('vertices', 'edges', 'faces'):
            old = getattr(self, kind)
            new = getattr(book, kind)
            reverse = {}
            for label, address in new.items():
                reverse.setdefault(address, label)

            for label, address in old.items():
                match = reverse.get(address)
                if match is None and tol:
                    match = _nearest(address, reverse, tol)
                if match is None:
                    report['lost'].append(label)
                    continue
                report['kept'][label] = match
                if match != label:
                    report['moved'].append(label)

        report['lost'].sort()
        report['moved'].sort()
        return book, report

    # ------------------------------------------------------------------

    def counts(self):
        return {'vertices': len(self.vertices), 'edges': len(self.edges),
                'faces': len(self.faces)}

    def __repr__(self):
        c = self.counts()
        return '<AddressBook {} vertices, {} edges, {} faces>'.format(
            c['vertices'], c['edges'], c['faces'])


def _nearest(address, reverse, tol):
    """The label whose address is within ``tol`` of this one, or ``None``."""
    best_label, best_d = None, tol
    for other, label in reverse.items():
        d = max(abs(address[i] - other[i]) for i in range(3))
        if d <= best_d:
            best_label, best_d = label, d
    return best_label
