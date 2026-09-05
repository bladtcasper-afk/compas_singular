"""**Getting an API key from somewhere that is not a repository.**

A key in the environment is the right default for a terminal, and the wrong one
for Rhino: Rhino reads the environment once when it starts, so a variable set
afterwards never reaches it, and telling a user to restart Rhino to change a key
is a poor trade. So a key may also live in a plain text file whose PATH is
configured, and the path is the only thing that ever appears in code.

**The file belongs outside every repository**, and that is not a formality --
``compas_singular`` pushes to a public fork, and ``compas_topology`` is not
version controlled today but could be tomorrow. A key committed once is a key
that has to be revoked. Somewhere under the user profile is the obvious home::

    C:\\Users\\<you>\\.secrets\\gemini_api_key.txt

The file is read as text and the first line that is not blank and not a comment
is taken. ``KEY=value`` is accepted as well as a bare value, because that is how
half of them arrive by email::

    # Google AI Studio, created 2026-08-28
    GEMINI_API_KEY=AIza...

Nothing here caches. A key read at construction is the key that run uses, and
editing the file changes the next run rather than needing anything restarted.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os


__all__ = ['read_key_file', 'resolve_api_key', 'KeyNotFound']


class KeyNotFound(Exception):
    """No key could be found, with a message saying everywhere that was tried."""


def read_key_file(path):
    """The key in a text file. Raises :class:`KeyNotFound` with a clear reason.

    Blank lines and ``#`` comments are skipped; a ``NAME=value`` line gives up
    its value. Surrounding quotes and whitespace are stripped, because a key
    pasted from a console often arrives wearing them.
    """
    expanded = os.path.expanduser(os.path.expandvars(str(path)))
    if not os.path.isfile(expanded):
        raise KeyNotFound(
            'no key file at {!r}{}. Create it with the key on one line, or '
            'point the path somewhere else.'.format(
                expanded,
                ' (expanded from {!r})'.format(path) if expanded != str(path)
                else ''))

    try:
        with open(expanded) as handle:
            lines = handle.read().splitlines()
    except IOError as exc:
        raise KeyNotFound('could not read {!r}: {}'.format(expanded, exc))

    for line in lines:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        if '=' in line:
            line = line.split('=', 1)[1].strip()
        line = line.strip('"').strip("'").strip()
        if line:
            return line

    raise KeyNotFound(
        '{!r} exists but has no key in it -- every line is blank or a '
        'comment.'.format(expanded))


def resolve_api_key(key=None, key_file=None, env=(), required=True):
    """A key from, in order: the argument, a file, the environment.

    Parameters
    ----------
    key : str, optional
        Wins outright. For a caller that already has one.
    key_file : str, optional
        Path to a text file. See :func:`read_key_file`. A path that is given but
        unusable RAISES rather than falling through to the environment -- being
        handed a stale environment key when the file you edited was not read is
        worse than being told the path is wrong.
    env : sequence[str], optional
        Variable names to try, in order.
    required : bool, optional
        Raise when nothing is found. ``False`` returns ``None`` instead, which
        lets an SDK apply its own resolution.

    Returns
    -------
    str or None
    """
    if key:
        return key
    if key_file:
        return read_key_file(key_file)
    for name in env:
        value = os.environ.get(name)
        if value:
            return value
    if not required:
        return None
    raise KeyNotFound(
        'no API key. Set {}, or pass a key_file pointing at a text file that '
        'holds one.'.format(' or '.join(env) if env else 'the API key variable'))
