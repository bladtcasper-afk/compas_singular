"""A cache for the frame-field solve that knows when it is stale.

``FieldDecomposition.from_boundary`` is the expensive step of the whole front
end -- background triangulation, sparse field solve, symmetrisation, separatrix
tracing -- and with ``build=True`` this module also carries the ``_build()``
stage after it, which ``19_field_accuracy.py`` measures as the larger cost of
the two. Recorded solve times run 0.4-2.1 s; a Rhino user stepping
``CMD_coarse_mesh`` -> ``CMD_edit_coarse_mesh`` -> ``CMD_quad_mesh`` on an
unchanged document paid for three identical solves before this existed, and paid
again after every restart.

**Nothing here imports Rhino.** The predecessor
(``compas_topology/rhino_commands/old/ff_common.py::get_decomposition``) cached
into ``scriptcontext.sticky``, so it could not be used from a script, a test, or
the ``CMD_`` family. This one is plain Python plus compas, and the CAD layer
supplies nothing but the boundary points.

Caching this solve at all is only legitimate because it is DETERMINISTIC, and
that is a property the package defends rather than merely hopes for. There is no
RNG anywhere in ``framefield``; the one place randomness would naturally enter --
breaking Qhull's tie on the cocircular corners of a regular seed grid -- is
``background._jitter``, a pure integer hash of the grid indices written for
exactly this reason. ``symmetry.py``, ``arrangement.py`` and ``densify.py`` each
assert byte-identical output in their own docstrings.


The invalidation contract
-------------------------

The cache key is the digest of FOUR sub-digests, and every one of them is
load-bearing. Three of the four are answers to a way the old cache was wrong.

``inputs``
    The outer boundary, the holes and the guides, each **tagged with its role**
    and hashed separately. The old key flattened ``[outer] + inners + guides``
    into one list, so the same polylines redistributed between holes and guides
    hashed identically -- a different problem served the wrong field. Points are
    rounded to :data:`COORDINATE_DIGITS`, because CAD round-trips shed float
    noise that must not thrash the key.

``params``
    Every ``from_boundary`` argument that reaches the solver, listed in
    :data:`SOLVE_PARAMETERS`. The old key carried ``target_length``, ``mode`` and
    ``symmetry`` only, so ``relax``, ``guide_weight``, ``guide_band`` and
    ``field_tau`` all changed the field without changing the key. ``relax`` is
    not a theoretical worry: the ``CMD_`` family passes
    ``relax=resolve_relax(settings, guides)``, which returns a different value
    for the same settings depending on whether the document has guides.

``code``
    sha256 over the CONTENTS of every ``framefield/*.py``. Without this the cache
    serves a field produced by code that no longer exists, silently, which is the
    failure that matters most while the library is being developed. Content and
    not mtime, so a git checkout restoring identical bytes does not invalidate.
    Memoised per process behind an mtime/size stamp, so a loop pays for the read
    once. ``cache.py`` excludes itself -- it does not participate in a solve, so
    editing it cannot change one.

``env``
    :data:`CACHE_VERSION`, the Python minor version, the numpy / scipy / compas
    versions, and ``TOL.precision``. scipy is in here because Qhull is the one
    external decision in the pipeline. ``TOL`` is in here because it is process-
    mutable global state that ``edges_to_curves`` and ``densification`` key
    geometry off.

Storing the four separately is what lets a miss say WHY it missed. The sidecar
beside every entry holds them plus the readable parameter and environment dicts,
so :meth:`SolveCache.get` can answer with ``framefield source changed since these
inputs were last solved`` or ``target_length 0.5 -> 0.3`` rather than re-solving
without comment.


Why every hit is a fresh object
-------------------------------

Both tiers store PICKLED BYTES and a hit is ``pickle.loads``. Two things fall out
of that, and both are the point rather than a side effect.

The tiers cannot drift apart: memory and disk reconstruct by the same code path,
so there is no class of bug that appears only after a restart.

And no two callers ever share an object. A ``FieldDecomposition`` is mutable in
ways that outlive a command -- ``edit_coarse`` sets ``_edited``, replaces
``mesh`` and rewrites ``polylines``; ``quad_mesh`` sets ``dense``;
``decomposition_mesh`` clears and refills ``repair_notes``. Handing one instance
to consecutive commands is what forced ``ff_common.print_warnings`` to
de-duplicate warnings and what forced three separate files to reset
``edit_notes`` defensively on the way in. ``CMD_quad_mesh`` never grew either
workaround, because it has always had a fresh object. It still does.

Pickle is normally fragile across a code change. Here it cannot be: the ``code``
axis means a pickle written by older ``framefield`` is never looked up by newer
``framefield``, because the key does not match. Everything on a decomposition is
plain Python plus one compas ``Mesh`` -- no numpy arrays, no scipy objects, no
file handles, no CAD references -- so nothing resists it. The store sits behind
:class:`SolveCache` so a compas-JSON store can replace it without touching a
caller; two things such a store would have to get right are recorded in
``25_cache.py``.

**Only ever unpickle from a directory you own.** :func:`default_directory` picks
a per-user cache directory for that reason, and nothing here will read an entry
from anywhere else unless you pass the path yourself.


Using it
--------

>>> from compas_singular.framefield import cache              # doctest: +SKIP
>>> d = cache.solve(outer, holes, guides=guides,              # doctest: +SKIP
...                 target_length=0.5, verbose=True)

``solve`` does not touch ``FieldDecomposition.from_boundary`` -- it is a sibling,
not a wrapper, and nothing is monkey-patched. Every existing script and harness
goes on solving cold, which matters: ``24_symmetry.py`` part 3 proves determinism
by solving the same input twice, and a cache underneath it would make that a
tautology.
"""
import hashlib
import json
import os
import pickle
import sys
import tempfile
import time
from collections import OrderedDict
from numbers import Number

from .constraints import as_curve_list


__all__ = [
    'CACHE_VERSION',
    'SolveCache',
    'CACHE',
    'solve',
    'default_directory',
    'source_digest',
]


#: Bump to invalidate every entry ever written, for a reason the four key axes
#: cannot see by themselves -- a change in what is stored rather than in what is
#: computed.
CACHE_VERSION = 1

#: Coordinate rounding for the ``inputs`` digest. Well below anything geometric
#: in the pipeline (``TOL.geometric_key`` works at 3) but coarse enough to absorb
#: the last-bit noise a CAD curve round-trip leaves on a point that has not
#: moved. The old cache used 6; 9 is the same idea with more headroom, since a
#: false MISS only costs a re-solve while a false HIT is wrong.
COORDINATE_DIGITS = 9

#: Everything ``from_boundary`` accepts that can change the solved field. Kept
#: explicit rather than harvested from the signature: a parameter added upstream
#: should fail loudly here (it lands in ``extra`` below and still enters the
#: key) rather than be silently dropped from it.
SOLVE_PARAMETERS = (
    'mode',
    'target_length',
    'orthogonal',
    'guide_weight',
    'guide_band',
    'relax',
    'field_tau',
    'symmetry',
)

#: Pickle protocol. 4 rather than the highest available, so an entry stays
#: readable by any Python this package runs under -- Rhino 8 ships CPython 3.9.
#: The ``env`` axis pins the minor version anyway; this is belt and braces.
PROTOCOL = 4

DEFAULT_MEMORY_BYTES = 256 * 1024 * 1024
DEFAULT_DISK_BYTES = 2 * 1024 * 1024 * 1024
DEFAULT_MAX_AGE_DAYS = 90

#: Environment variable that overrides :func:`default_directory`. Set it to an
#: empty string to disable the disk tier entirely.
DIRECTORY_ENVVAR = 'COMPAS_SINGULAR_CACHE'


# ==============================================================================
# where entries live
# ==============================================================================

def default_directory():
    """The per-user directory entries are written to.

    ``$COMPAS_SINGULAR_CACHE`` wins if set; an empty value disables the disk
    tier. Otherwise ``%LOCALAPPDATA%`` on Windows, ``$XDG_CACHE_HOME`` or
    ``~/.cache`` elsewhere, and the system temp directory if neither resolves.

    Deliberately NOT ``compas_singular.TEMP``. That path is computed relative to
    the source tree, so under a non-editable install it lands inside
    ``site-packages`` -- a directory the user may not own and may not be able to
    write, which is the wrong place to be unpickling from.

    Returns
    -------
    str or None
        ``None`` when the disk tier is switched off.
    """
    override = os.environ.get(DIRECTORY_ENVVAR)
    if override is not None:
        override = override.strip()
        return override or None

    root = None
    if sys.platform == 'win32':
        root = os.environ.get('LOCALAPPDATA')
    if not root:
        root = os.environ.get('XDG_CACHE_HOME')
    if not root:
        home = os.path.expanduser('~')
        if home and home != '~':
            root = os.path.join(home, '.cache')
    if not root:
        root = tempfile.gettempdir()
    return os.path.join(root, 'compas_singular', 'framefield-solve')


# ==============================================================================
# the four key axes
# ==============================================================================

def _digest(obj):
    blob = json.dumps(obj, sort_keys=True, separators=(',', ':'))
    return hashlib.sha256(blob.encode('utf-8')).hexdigest()


def _as_curve(curve):
    """One curve as a plain, rounded list of ``[x, y, z]``.

    Accepts what ``from_boundary`` accepts -- lists, tuples, compas ``Point``s,
    a compas ``Polyline`` -- and drops a repeated closing point, because a loop
    given closed and the same loop given open are the same input to the solver
    (``background._as_open_loop`` drops it too).
    """
    points = []
    for point in curve:
        points.append([round(float(point[i]), COORDINATE_DIGITS) for i in range(3)])
    if len(points) > 1 and points[0] == points[-1]:
        points = points[:-1]
    return points


def _inputs(outer_boundary, inner_boundaries, guides):
    """The geometry, with each loop's ROLE part of what is hashed.

    The role tags are what stop a hole and a guide with identical points from
    colliding. Holes and guides keep their given ORDER, because
    ``Symmetry.detect`` and ``from_curves`` both consume them as sequences.
    """
    return {
        'outer': _as_curve(outer_boundary),
        'inners': [_as_curve(loop) for loop in (inner_boundaries or [])],
        'guides': [_as_curve(curve) for curve in as_curve_list(guides)],
    }


def _plain(value):
    """A solve parameter as something JSON can hold and a human can read."""
    if value is None or isinstance(value, bool) or isinstance(value, str):
        return value
    if isinstance(value, Number):
        return float(value)
    # A Symmetry. Its group IS the identity -- centre, elements and the steps it
    # is enabled for -- and none of that is reconstructible from ``repr``.
    if hasattr(value, 'centre') and hasattr(value, 'steps'):
        return {
            'symmetry': [round(float(c), COORDINATE_DIGITS) for c in value.centre],
            'names': sorted(value.names()),
            'steps': sorted(value.steps),
        }
    if isinstance(value, (list, tuple)):
        return [_plain(item) for item in value]
    if isinstance(value, dict):
        return dict((str(k), _plain(v)) for k, v in value.items())
    return repr(value)


def _parameters(kwargs):
    """The readable parameter dict that goes into the ``params`` digest.

    Only what the caller actually PASSED is recorded, so omitting a parameter
    and passing its default key differently even though they solve the same.
    That is deliberate: pinning the defaults here would mean a second copy of
    ``from_boundary``'s signature to keep in step, and the cost of getting this
    wrong is one extra solve rather than a wrong answer. Anything
    ``from_boundary`` grows that is not yet in :data:`SOLVE_PARAMETERS` enters
    the key too, under ``extra``, for the same reason: an unknown parameter is
    assumed to matter.
    """
    out = {}
    for name in SOLVE_PARAMETERS:
        if name in kwargs:
            out[name] = _plain(kwargs[name])
    extra = dict((k, _plain(v)) for k, v in kwargs.items()
                 if k not in SOLVE_PARAMETERS)
    if extra:
        out['extra'] = extra
    return out


_SOURCE = {'stamp': None, 'digest': None}


def source_digest():
    """sha256 over the contents of every ``framefield/*.py`` but this one.

    Memoised behind an (name, mtime, size) stamp, so a harness in a loop reads
    the sources once rather than once per lookup.
    """
    here = os.path.dirname(os.path.abspath(__file__))
    mine = os.path.basename(__file__)
    names = sorted(name for name in os.listdir(here)
                   if name.endswith('.py') and name != mine)

    stamp = []
    for name in names:
        path = os.path.join(here, name)
        try:
            info = os.stat(path)
        except OSError:                     # deleted between listdir and stat
            continue
        stamp.append((name, info.st_mtime, info.st_size))
    stamp = tuple(stamp)
    if stamp == _SOURCE['stamp']:
        return _SOURCE['digest']

    hasher = hashlib.sha256()
    for name, _, _ in stamp:
        hasher.update(name.encode('utf-8'))
        with open(os.path.join(here, name), 'rb') as handle:
            hasher.update(handle.read())
    _SOURCE['stamp'] = stamp
    _SOURCE['digest'] = hasher.hexdigest()
    return _SOURCE['digest']


def environment():
    """The versions and global tolerances the solve is reproducible under."""
    import compas
    import numpy
    import scipy
    from compas.tolerance import TOL

    return {
        'cache_version': CACHE_VERSION,
        'python': '{}.{}'.format(*sys.version_info[:2]),
        'compas': compas.__version__,
        'numpy': numpy.__version__,
        'scipy': scipy.__version__,
        'tol_precision': TOL.precision,
        'tol_absolute': TOL.absolute,
    }


class CacheKey(object):
    """The four axes, plus the digest that is their conjunction.

    Attributes
    ----------
    digest : str
        The entry's filename stem.
    inputs, params, code, env : str
        The sub-digests, stored in the sidecar so a miss can name its cause.
    parameters, environment : dict
        The readable forms, stored so a miss can name the value that changed.
    """

    def __init__(self, inputs, parameters, env):
        self.parameters = parameters
        self.environment = env
        self.inputs = _digest(inputs)
        self.params = _digest(parameters)
        self.code = source_digest()
        self.env = _digest(env)
        self.digest = _digest([self.inputs, self.params, self.code, self.env])

    def sidecar(self):
        return {
            'inputs': self.inputs,
            'params': self.params,
            'code': self.code,
            'env': self.env,
            'parameters': self.parameters,
            'environment': self.environment,
        }

    def __repr__(self):
        return '<CacheKey {}>'.format(self.digest[:12])


# ==============================================================================
# the cache
# ==============================================================================

class SolveCache(object):
    """A two-tier store of pickled ``FieldDecomposition``s.

    Parameters
    ----------
    directory : str or None, optional
        Where entries are written. ``None`` takes :func:`default_directory`;
        pass ``False`` for a memory-only cache.
    memory_bytes, disk_bytes : int, optional
        LRU budgets. Both count the pickled size, which is the only figure that
        is knowable without reconstructing.
    max_age_days : float, optional
        Entries older than this are pruned on the next write. A stale entry is
        never WRONG -- the key would have to match for it to be served, and a
        matching key means the inputs, the parameters, the source and the
        environment are all unchanged -- so this is housekeeping, not safety.
    """

    def __init__(self, directory=None, memory_bytes=DEFAULT_MEMORY_BYTES,
                 disk_bytes=DEFAULT_DISK_BYTES, max_age_days=DEFAULT_MAX_AGE_DAYS):
        if directory is False:
            self.directory = None
        elif directory is None:
            self.directory = default_directory()
        else:
            self.directory = directory
        self.memory_bytes = memory_bytes
        self.disk_bytes = disk_bytes
        self.max_age_days = max_age_days
        self._memory = OrderedDict()
        self._memory_size = 0

    # ------------------------------------------------------------------
    # keys
    # ------------------------------------------------------------------

    def key(self, outer_boundary, inner_boundaries=None, guides=None, **kwargs):
        """The :class:`CacheKey` for one call to ``from_boundary``."""
        return CacheKey(_inputs(outer_boundary, inner_boundaries, guides),
                        _parameters(kwargs), environment())

    # ------------------------------------------------------------------
    # lookup
    # ------------------------------------------------------------------

    def get(self, key):
        """``(decomposition, reason)``; ``decomposition`` is ``None`` on a miss.

        ``reason`` is always a phrase fit to print -- ``'in memory'`` or
        ``'on disk'`` on a hit, and on a miss the axis that moved.
        """
        blob = self._memory.get(key.digest)
        if blob is not None:
            self._memory.move_to_end(key.digest)
            decomposition = self._loads(blob)
            if decomposition is not None:
                return decomposition, 'in memory'
            self._memory_drop(key.digest)

        path = self._path(key.digest)
        if path and os.path.isfile(path):
            try:
                with open(path, 'rb') as handle:
                    blob = handle.read()
            except OSError:
                blob = None
            decomposition = self._loads(blob) if blob else None
            if decomposition is not None:
                os.utime(path, None)                # keep the LRU honest
                self._memory_put(key.digest, blob)
                return decomposition, 'on disk'
            # Truncated by a killed Rhino, or written by something else. A
            # damaged entry is a miss and never an exception into a command.
            self._discard(key.digest)

        return None, self._reason(key)

    def _loads(self, blob):
        try:
            return pickle.loads(blob)
        except Exception:                           # noqa: BLE001 -- any of it is a miss
            return None

    # ------------------------------------------------------------------
    # why a miss missed
    # ------------------------------------------------------------------

    def _reason(self, key):
        entries = self.entries()
        if not entries:
            return 'nothing cached yet'

        same_inputs = [e for e in entries if e.get('inputs') == key.inputs]
        if same_inputs:
            for entry in same_inputs:
                if entry.get('params') != key.params:
                    continue
                if entry.get('code') != key.code:
                    return 'framefield source changed since these inputs were last solved'
                if entry.get('env') != key.env:
                    return 'environment changed -- {}'.format(
                        _diff(entry.get('environment'), key.environment))
            return 'solver settings changed -- {}'.format(
                _diff(same_inputs[0].get('parameters'), key.parameters))

        for entry in entries:
            if (entry.get('params') == key.params
                    and entry.get('code') == key.code
                    and entry.get('env') == key.env):
                return 'the input geometry changed'
        return 'no entry for these inputs'

    # ------------------------------------------------------------------
    # storing
    # ------------------------------------------------------------------

    def put(self, key, decomposition, seconds=None):
        """Store ``decomposition`` under ``key``. Returns the pickled size."""
        blob = pickle.dumps(decomposition, protocol=PROTOCOL)
        self._memory_put(key.digest, blob)

        path = self._path(key.digest)
        if path:
            meta = key.sidecar()
            meta['created'] = time.time()
            meta['seconds'] = None if seconds is None else round(float(seconds), 3)
            meta['bytes'] = len(blob)
            try:
                _ensure_directory(self.directory)
                _write_atomic(path, blob)
                _write_atomic(self._path(key.digest, '.json'),
                              json.dumps(meta, sort_keys=True, indent=1).encode('utf-8'))
            except OSError:
                # A read-only or full cache directory degrades to memory-only.
                # It must never take the command down with it.
                self._discard(key.digest)
            else:
                self.prune()
        return len(blob)

    # ------------------------------------------------------------------
    # housekeeping
    # ------------------------------------------------------------------

    def entries(self):
        """Every stored entry's sidecar, newest first. Nothing is unpickled."""
        if not self.directory or not os.path.isdir(self.directory):
            return []
        out = []
        for name in os.listdir(self.directory):
            if not name.endswith('.json'):
                continue
            path = os.path.join(self.directory, name)
            try:
                with open(path, 'r') as handle:
                    meta = json.load(handle)
            except (OSError, ValueError):
                continue
            meta['digest'] = name[:-len('.json')]
            out.append(meta)
        out.sort(key=lambda meta: meta.get('created') or 0.0, reverse=True)
        return out

    def prune(self):
        """Drop entries past the age limit, then oldest-first past the budget."""
        if not self.directory or not os.path.isdir(self.directory):
            return 0
        found = []
        for name in os.listdir(self.directory):
            if not name.endswith('.pickle'):
                continue
            path = os.path.join(self.directory, name)
            try:
                info = os.stat(path)
            except OSError:
                continue
            found.append((info.st_mtime, info.st_size, name[:-len('.pickle')]))
        found.sort()

        dropped = 0
        cutoff = time.time() - self.max_age_days * 86400.0
        keep = []
        for mtime, size, digest in found:
            if mtime < cutoff:
                self._discard(digest)
                dropped += 1
            else:
                keep.append((mtime, size, digest))

        total = sum(size for _, size, _ in keep)
        index = 0
        while total > self.disk_bytes and index < len(keep):
            _, size, digest = keep[index]
            self._discard(digest)
            total -= size
            dropped += 1
            index += 1
        return dropped

    def clear(self):
        """Empty both tiers. Returns how many disk entries went."""
        self._memory.clear()
        self._memory_size = 0
        if not self.directory or not os.path.isdir(self.directory):
            return 0
        dropped = 0
        for name in os.listdir(self.directory):
            if name.endswith('.pickle'):
                self._discard(name[:-len('.pickle')])
                dropped += 1
        return dropped

    # ------------------------------------------------------------------
    # internals
    # ------------------------------------------------------------------

    def _path(self, digest, suffix='.pickle'):
        if not self.directory:
            return None
        return os.path.join(self.directory, digest + suffix)

    def _discard(self, digest):
        for suffix in ('.pickle', '.json'):
            path = self._path(digest, suffix)
            if path and os.path.isfile(path):
                try:
                    os.remove(path)
                except OSError:
                    pass
        self._memory_drop(digest)

    def _memory_put(self, digest, blob):
        self._memory_drop(digest)
        if len(blob) > self.memory_bytes:
            return
        self._memory[digest] = blob
        self._memory_size += len(blob)
        while self._memory_size > self.memory_bytes and self._memory:
            _, evicted = self._memory.popitem(last=False)
            self._memory_size -= len(evicted)

    def _memory_drop(self, digest):
        blob = self._memory.pop(digest, None)
        if blob is not None:
            self._memory_size -= len(blob)


def _diff(before, after):
    """The first differing entry of two readable dicts, as ``'name a -> b'``."""
    before = before or {}
    after = after or {}
    for name in sorted(set(before) | set(after)):
        old = before.get(name, '<unset>')
        new = after.get(name, '<unset>')
        if old != new:
            return '{} {} -> {}'.format(name, old, new)
    return 'no visible difference'


def _ensure_directory(path):
    if not os.path.isdir(path):
        os.makedirs(path)


def _write_atomic(path, payload):
    """Write via a sibling temp file and ``os.replace``.

    A Rhino killed mid-write must not leave a half entry that the next session
    would try to unpickle. ``os.replace`` is atomic on Windows and POSIX alike.
    """
    temporary = path + '.{}.tmp'.format(os.getpid())
    with open(temporary, 'wb') as handle:
        handle.write(payload)
    os.replace(temporary, path)


#: The cache :func:`solve` uses when none is passed.
CACHE = SolveCache()


# ==============================================================================
# the entry point
# ==============================================================================

def solve(outer_boundary, inner_boundaries=None, guides=None, cache=None,
          force=False, verbose=False, build=False, **kwargs):
    """A ``FieldDecomposition`` for these inputs, solved or recalled.

    Parameters
    ----------
    outer_boundary : list[[x, y, z]]
    inner_boundaries : list[list[[x, y, z]]], optional
    guides : list[list[[x, y, z]]], optional
    cache : SolveCache or False, optional
        ``None`` uses the module-level :data:`CACHE`; ``False`` disables caching
        for this call, which is what a determinism check wants.
    force : bool, optional
        Re-solve and overwrite the entry. For when you suspect the cache rather
        than the code -- ordinarily the key should be trusted to notice.
    verbose : bool, optional
        Print a line saying what happened and, on a miss, why.
    build : bool, optional
        Run ``_build()`` -- ``build_network`` plus ``planar_arrangement`` --
        before storing, so that stage rides along in the entry. It takes no
        arguments and is a pure function of the solve, so this changes nothing
        about the result. Worth it for a caller that will ask for a layout;
        wasted on one that only wants the field.
    **kwargs
        Passed through to ``FieldDecomposition.from_boundary`` and hashed. See
        :data:`SOLVE_PARAMETERS`.

    Returns
    -------
    FieldDecomposition
        Freshly reconstructed on every hit, so it is yours to mutate.
    """
    from .decomposition import FieldDecomposition

    store = CACHE if cache is None else cache
    key = store.key(outer_boundary, inner_boundaries, guides, **kwargs) if store else None

    if key is not None and not force:
        decomposition, reason = store.get(key)
        if decomposition is not None:
            if verbose:
                print('field: reusing a cached solve ({})'.format(reason))
            return decomposition
        if verbose:
            print('field: cache miss -- {}'.format(reason))

    start = time.time()
    decomposition = FieldDecomposition.from_boundary(
        outer_boundary, inner_boundaries, guides=guides, **kwargs)
    if build:
        decomposition._build()
    seconds = time.time() - start
    if verbose:
        print('field: solved and traced in {:.1f} s'.format(seconds))

    if key is not None:
        store.put(key, decomposition, seconds=seconds)
    return decomposition


# ==============================================================================
# Main
# ==============================================================================

if __name__ == '__main__':
    pass
