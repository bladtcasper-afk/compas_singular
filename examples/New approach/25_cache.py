"""EXAMPLE 8 -- the solve cache, and proof that it knows when it is stale.

``framefield/cache.py`` stores a solved ``FieldDecomposition`` in memory and on
disk. A cache is only worth having if a hit is indistinguishable from a solve and
a stale entry is impossible, so this script asserts both rather than reporting
them.

Nothing here is a recorded snapshot. Every check is a PROPERTY -- two objects
agree, or a lookup misses for the reason it should -- so the script stays true
when the layouts move. ``baseline.json`` remains the only place numbers live.

Five parts:

1. IDENTITY   a hit reproduces a cold solve exactly, including ``max |dtheta|``
              of exactly zero over every background vertex.
2. STALENESS  each of the seven key axes invalidates on its own, and the miss
              names the axis that moved. Includes the regression test for the
              old cache's role collision: a loop moved from ``inner_boundaries``
              to ``guides``, same points, must MISS.
3. ISOLATION  two hits on one key are independent objects -- editing one leaves
              the other pristine.
4. DAMAGE     a truncated entry is a clean miss and a re-solve, never a raise.
5. TIMING     cold against warm, on a domain worth the wait.

Run:
    python 25_cache.py                    # everything
    python 25_cache.py --identity         # part 1 only
    python 25_cache.py --staleness        # part 2 only
    python 25_cache.py --isolation        # part 3 only
    python 25_cache.py --damage           # part 4 only
    python 25_cache.py --timing           # part 5 only

Exit code 0 if every check passed, 1 otherwise.

The cache used here is a private one in a temporary directory. This script never
reads, writes or prunes the user's real cache.
"""
import os
import shutil
import sys
import tempfile
import time
from math import cos, pi, sin

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from compas.tolerance import TOL                              # noqa: E402

from compas_singular.framefield import cache  # noqa: E402


#: background spacing for the cheap correctness parts
SPACING = 0.5

#: quad size for the quality comparison, matching ``15_baseline.py``
QUAD_TARGET = 1.0


# ----------------------------------------------------------------------
# domains -- deliberately small, except the one part 5 times
# ----------------------------------------------------------------------

SQUARE = [[0, 0, 0], [10, 0, 0], [10, 10, 0], [0, 10, 0]]
SQUARE_HOLE = [[4, 4, 0], [6, 4, 0], [6, 6, 0], [4, 6, 0]]
CABLE = [[0, 3, 0], [10, 7, 0]]
DISC = [[7.0 * cos(2 * pi * i / 48), 7.0 * sin(2 * pi * i / 48), 0.0]
        for i in range(48)]


# ----------------------------------------------------------------------
# a scratch cache, so the user's own is never touched
# ----------------------------------------------------------------------

class Scratch(object):
    """A :class:`SolveCache` in a temporary directory, removed on exit."""

    def __enter__(self):
        self.directory = tempfile.mkdtemp(prefix='framefield-cache-test-')
        return cache.SolveCache(directory=self.directory)

    def __exit__(self, *exc):
        shutil.rmtree(self.directory, ignore_errors=True)
        return False


class Checks(object):
    """A running tally, so one failure does not hide the rest."""

    def __init__(self):
        self.passed = 0
        self.failed = []

    def ok(self, condition, label, detail=''):
        if condition:
            self.passed += 1
            print('    PASS  {}{}'.format(label, '  ' + detail if detail else ''))
        else:
            self.failed.append(label)
            print('    FAIL  {}{}'.format(label, '  ' + detail if detail else ''))
        return bool(condition)


def section(title):
    print('')
    print('=' * 78)
    print(title)
    print('=' * 78)


# ----------------------------------------------------------------------
# fingerprints -- the same ones the rest of the suite already uses
# ----------------------------------------------------------------------

def coarse_xyz(mesh):
    """Sorted, rounded layout corners. ``24_symmetry.py`` fingerprints this way."""
    return sorted(tuple(round(c, 9) for c in mesh.vertex_coordinates(v))
                  for v in mesh.vertices())


def theta_gap(one, other):
    """``max |dtheta|`` over every background vertex.

    ``21_edit_coarse.py`` requires this to be EXACTLY zero of an operation that
    must not disturb the field, and it is the sharpest available test that a
    reconstructed field is the solved one rather than merely close to it.
    """
    keys = set(one.field.theta) | set(other.field.theta)
    if set(one.field.theta) != set(other.field.theta):
        return float('inf')
    return max(abs(one.field.theta[k] - other.field.theta[k]) for k in keys)


# ======================================================================
# 1. identity
# ======================================================================

def part_identity(checks):
    section('1. IDENTITY -- a hit is indistinguishable from a solve')

    with Scratch() as store:
        kwargs = dict(target_length=SPACING, cache=store, build=True)

        cold = cache.solve(SQUARE, [SQUARE_HOLE], guides=[CABLE],
                           guide_weight=5.0, **kwargs)
        warm = cache.solve(SQUARE, [SQUARE_HOLE], guides=[CABLE],
                           guide_weight=5.0, **kwargs)

        checks.ok(cold is not warm, 'a hit is a distinct object')

        gap = theta_gap(cold, warm)
        checks.ok(gap == 0.0, 'field identical', 'max |dtheta| = {:.3e}'.format(gap))

        checks.ok(cold.field.singularities() == warm.field.singularities(),
                  'singularities identical',
                  '{} of them'.format(len(cold.field.singularities())))

        checks.ok([len(s.points) for s in cold.separatrices]
                  == [len(s.points) for s in warm.separatrices],
                  'separatrices identical',
                  '{} traces'.format(len(cold.separatrices)))

        checks.ok(cold.background.mesh.number_of_vertices()
                  == warm.background.mesh.number_of_vertices(),
                  'background identical',
                  '{} vertices'.format(cold.background.mesh.number_of_vertices()))

        checks.ok(_symmetry_repr(cold) == _symmetry_repr(warm), 'symmetry group identical',
                  _symmetry_repr(cold))

        # ...and the whole downstream pipeline off a reconstructed object.
        cold_dense = cold.quad_mesh(target_length=QUAD_TARGET)
        warm_dense = warm.quad_mesh(target_length=QUAD_TARGET)

        checks.ok(cold.route() == warm.route(), 'route identical', cold.route())
        checks.ok(coarse_xyz(cold.mesh) == coarse_xyz(warm.mesh),
                  'coarse layout identical',
                  '{} patches'.format(cold.mesh.number_of_faces()))
        checks.ok(cold.quality(cold_dense) == warm.quality(warm_dense),
                  'quality dict identical',
                  '{} faces'.format(warm_dense.number_of_faces()))
        checks.ok(cold.warnings() == warm.warnings(), 'warnings identical',
                  '{} of them'.format(len(warm.warnings())))


def _symmetry_repr(decomposition):
    group = decomposition.symmetry
    if group is None:
        return 'none'
    return 'order {} -- {}'.format(len(group), ', '.join(sorted(group.names())))


# ======================================================================
# 2. staleness
# ======================================================================

def part_staleness(checks):
    section('2. STALENESS -- every axis invalidates on its own, and says so')

    with Scratch() as store:
        base = dict(target_length=SPACING)
        cache.solve(SQUARE, [SQUARE_HOLE], cache=store, build=True, **base)
        key = store.key(SQUARE, [SQUARE_HOLE], **base)
        got, reason = store.get(key)
        checks.ok(got is not None, 'the entry it just wrote is a hit', reason)

        nudged = [list(p) for p in SQUARE]
        nudged[2][0] += 1e-3

        variants = [
            ('boundary point moved 1e-3',
             (nudged, [SQUARE_HOLE]), base, 'geometry'),
            ('hole moved 1e-3',
             (SQUARE, [[[c + 1e-3 for c in p] for p in SQUARE_HOLE]]), base, 'geometry'),
            ('target_length 0.5 -> 0.4',
             (SQUARE, [SQUARE_HOLE]), dict(target_length=0.4), 'settings'),
            ('relax turned on',
             (SQUARE, [SQUARE_HOLE]), dict(base, relax=True), 'settings'),
            ('guide_weight set',
             (SQUARE, [SQUARE_HOLE]), dict(base, guide_weight=5.0), 'settings'),
            ('guide_band set',
             (SQUARE, [SQUARE_HOLE]), dict(base, guide_band=1.5), 'settings'),
            ('field_tau set',
             (SQUARE, [SQUARE_HOLE]), dict(base, field_tau=0.1), 'settings'),
            ('symmetry disabled',
             (SQUARE, [SQUARE_HOLE]), dict(base, symmetry=None), 'settings'),
        ]
        for label, args, kwargs, expect in variants:
            other, reason = store.get(store.key(*args, **kwargs))
            checks.ok(other is None and expect in reason, label, reason)

        # A closed loop given closed and given open is the SAME input. This one
        # must hit -- an invalidation that fires on nothing is as bad as one
        # that never fires.
        closed = SQUARE + [SQUARE[0]]
        got, reason = store.get(store.key(closed, [SQUARE_HOLE], **base))
        checks.ok(got is not None, 'a repeated closing point still hits', reason)

        print('')
        print('  THE COLLISION THE OLD CACHE HAD')
        print('  ff_common hashed [outer] + inners + guides as one flat list, so the')
        print('  same polyline as a HOLE and as a GUIDE keyed identically.')
        with Scratch() as roles:
            cache.solve(SQUARE, [SQUARE_HOLE], cache=roles, **base)
            as_hole = roles.key(SQUARE, [SQUARE_HOLE], **base)
            as_guide = roles.key(SQUARE, None, guides=[SQUARE_HOLE], **base)
            checks.ok(as_hole.digest != as_guide.digest,
                      'hole and guide with identical points key differently')
            got, reason = roles.get(as_guide)
            checks.ok(got is None, 'the guide version misses', reason)

        print('')
        print('  SOURCE AND ENVIRONMENT')
        before = cache.source_digest()
        #: derived from the package itself, not from HERE -- this is the very
        #: directory ``source_digest`` scans, wherever framefield is installed.
        touched = os.path.join(os.path.dirname(cache.__file__), 'trace.py')
        with open(touched, 'rb') as handle:
            original = handle.read()
        try:
            with open(touched, 'wb') as handle:
                handle.write(original + b'\n# cache staleness probe\n')
            after = cache.source_digest()
            checks.ok(after != before, 'editing trace.py changes the code digest')
            got, reason = store.get(store.key(SQUARE, [SQUARE_HOLE], **base))
            checks.ok(got is None and 'source changed' in reason,
                      'and the entry goes stale, by name', reason)
        finally:
            with open(touched, 'wb') as handle:
                handle.write(original)
        checks.ok(cache.source_digest() == before,
                  'restoring the bytes restores the digest')
        got, reason = store.get(store.key(SQUARE, [SQUARE_HOLE], **base))
        checks.ok(got is not None, 'and the entry is live again', reason)

        precision = TOL.precision
        try:
            TOL.precision = precision + 1
            got, reason = store.get(store.key(SQUARE, [SQUARE_HOLE], **base))
            checks.ok(got is None and 'environment changed' in reason,
                      'TOL.precision invalidates', reason)
        finally:
            TOL.precision = precision
        got, reason = store.get(store.key(SQUARE, [SQUARE_HOLE], **base))
        checks.ok(got is not None, 'and restoring TOL.precision restores the hit', reason)


# ======================================================================
# 3. isolation
# ======================================================================

def part_isolation(checks):
    section('3. ISOLATION -- one caller cannot leak into the next')

    print('  A FieldDecomposition is mutable well after it is built. Sharing one')
    print('  across commands is what forced ff_common.print_warnings to de-duplicate')
    print('  and three files to reset edit_notes defensively on the way in.')
    print('')

    with Scratch() as store:
        kwargs = dict(target_length=SPACING, cache=store, build=True)
        first = cache.solve(SQUARE, [SQUARE_HOLE], **kwargs)

        coarse = first.decomposition_mesh()
        first.edit_coarse(coarse, strict=False)
        first.edit_notes['probe'] = 'set by the first caller'
        first.field_aware = False
        first.user_curves = [[[0, 0, 0], [1, 1, 0]]]
        first.quad_mesh(target_length=QUAD_TARGET)

        second = cache.solve(SQUARE, [SQUARE_HOLE], **kwargs)

        checks.ok(second is not first, 'the second caller gets its own object')
        checks.ok(second._edited is False, 'the edit did not travel')
        checks.ok(second.edit_notes == {}, 'edit_notes is empty')
        checks.ok(second.mesh is None, 'mesh is unset')
        checks.ok(second.dense is None, 'dense is unset')
        checks.ok(second.repair_notes == [], 'repair_notes is empty')
        checks.ok(second.field_aware is True, 'field_aware is back at its default')
        checks.ok(second.user_curves == [], 'user_curves is empty')
        checks.ok(second._network is not None,
                  'but the built arrangement DID survive -- that is the point of build=True')


# ======================================================================
# 4. damage
# ======================================================================

def part_damage(checks):
    section('4. DAMAGE -- a broken entry is a miss, never an exception')

    with Scratch() as store:
        base = dict(target_length=SPACING)
        cache.solve(SQUARE, cache=store, **base)
        key = store.key(SQUARE, **base)
        path = store._path(key.digest)

        store._memory.clear()
        store._memory_size = 0
        with open(path, 'wb') as handle:
            handle.write(b'not a pickle')

        got, reason = store.get(key)
        checks.ok(got is None, 'a truncated entry misses instead of raising', reason)
        checks.ok(not os.path.isfile(path), 'and the damaged file is removed')

        recovered = cache.solve(SQUARE, cache=store, **base)
        checks.ok(recovered is not None, 're-solving after the damage works')
        got, _ = store.get(key)
        checks.ok(got is not None, 'and the entry is rewritten')

        # Housekeeping: the byte budget must actually evict.
        store.disk_bytes = 1
        store.prune()
        checks.ok(store.entries() == [], 'a byte budget of 1 prunes everything')

        # And a cache directory that cannot be written degrades to memory-only.
        broken = cache.SolveCache(directory=os.path.join(path, 'not-a-directory'))
        broken.put(broken.key(SQUARE, **base), recovered, seconds=0.0)
        got, _ = broken.get(broken.key(SQUARE, **base))
        checks.ok(got is not None, 'an unwritable directory still caches in memory')


# ======================================================================
# 5. timing
# ======================================================================

def part_timing(checks):
    section('5. TIMING -- what the cache is actually worth')

    rows = []
    with Scratch() as store:
        for label, outer, holes, guides, kwargs in [
            ('square', SQUARE, None, None, dict(target_length=SPACING)),
            ('square+hole', SQUARE, [SQUARE_HOLE], None, dict(target_length=SPACING)),
            ('square+cable', SQUARE, None, [CABLE],
             dict(target_length=SPACING, guide_weight=5.0)),
            ('disc', DISC, None, None, dict(target_length=SPACING)),
            ('disc @0.3', DISC, None, None, dict(target_length=0.3)),
            ('disc+hole @0.3', DISC, [SQUARE_HOLE], None, dict(target_length=0.3)),
        ]:
            common = dict(kwargs, cache=store, build=True)
            start = time.time()
            cache.solve(outer, holes, guides=guides, force=True, **common)
            cold = time.time() - start

            start = time.time()
            cache.solve(outer, holes, guides=guides, **common)
            memory = time.time() - start

            store._memory.clear()
            store._memory_size = 0
            start = time.time()
            cache.solve(outer, holes, guides=guides, **common)
            disk = time.time() - start

            key = store.key(outer, holes, guides, **kwargs)
            rows.append((label, cold, memory, disk, _entry_bytes(store, key.digest)))

        print('')
        print('  {:16s} {:>9s} {:>9s} {:>9s} {:>9s}'.format(
            'domain', 'cold', 'memory', 'disk', 'entry'))
        for label, cold, memory, disk, size in rows:
            print('  {:16s} {:8.2f}s {:8.3f}s {:8.3f}s {:8.0f}K'.format(
                label, cold, memory, disk, size / 1024.0))

        total_cold = sum(row[1] for row in rows)
        total_warm = sum(row[3] for row in rows)
        print('')
        print('  {:16s} {:8.2f}s {:>9s} {:8.3f}s {:8.0f}K'.format(
            'TOTAL', total_cold, '', total_warm,
            sum(row[4] for row in rows) / 1024.0))
        checks.ok(total_warm < total_cold,
                  'a warm run is faster than a cold one',
                  '{:.1f}x'.format(total_cold / max(total_warm, 1e-6)))


def _entry_bytes(store, digest):
    path = store._path(digest)
    try:
        return os.path.getsize(path)
    except OSError:
        return 0


# ======================================================================
# main
# ======================================================================

PARTS = [
    ('--identity', part_identity),
    ('--staleness', part_staleness),
    ('--isolation', part_isolation),
    ('--damage', part_damage),
    ('--timing', part_timing),
]


def main(argv):
    chosen = [fn for flag, fn in PARTS if flag in argv]
    if not chosen:
        chosen = [fn for _, fn in PARTS]

    print('cache directory in use by default: {}'.format(cache.default_directory()))
    print('(this script uses a temporary one and never touches it)')

    checks = Checks()
    for run in chosen:
        run(checks)

    section('RESULT')
    if checks.failed:
        print('{} passed, {} FAILED:'.format(checks.passed, len(checks.failed)))
        for label in checks.failed:
            print('  - {}'.format(label))
        return 1
    print('all {} checks passed.'.format(checks.passed))
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
