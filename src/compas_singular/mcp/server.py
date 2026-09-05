"""**Binding the tools, the session and the library to the protocol's methods.**

:mod:`~compas_singular.mcp.protocol` knows JSON-RPC and nothing about meshes.
This module is the handler it asks: it holds the one session the process is
working on, answers ``tools/*`` from the registry, and answers ``resources/*``
and ``prompts/*`` from the library.

**Importing this module is what registers the tools.** Each ``tools_*`` module
registers into :data:`~compas_singular.mcp.registry.TOOLS` at import time, so the
imports below are load-bearing rather than tidy -- dropping one silently removes
its tools from ``tools/list``. :func:`build` refuses to start on an empty
registry rather than letting that happen quietly.

**One session per process.** An MCP server over stdio serves exactly one client,
so there is nothing to key a session registry by and nothing to isolate. If this
ever grows a socket transport, that assumption is the thing to revisit.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from . import library
from . import registry
from . import tools_io         # noqa: F401  -- registers its tools
from . import tools_library    # noqa: F401  -- registers its tools
from . import tools_mesh       # noqa: F401  -- registers its tools
from . import tools_rhino      # noqa: F401  -- registers its tools
from .describe import describe
from .protocol import Dispatcher
from .session import MeshSession


__all__ = ['Handler', 'build', 'INSTRUCTIONS']


INSTRUCTIONS = """Improve an existing quad mesh from the compas_singular library.

Pull a mesh from the open Rhino document with rhino_pull (or work on one already
loaded), measure it with inspect, improve it with relax / smooth_region /
smooth_boundary_constrained, and push it back with rhino_push.

Read guidance://quality and guidance://smoothing before choosing a smoother --
they carry the measurements that decide which one is right, and the order to try
them in. inspect's share_below is what separates one bad face from a generally
slack mesh.

This version reads and improves. It does not solve a frame field and does not
change topology."""


class Handler(object):
    """What :class:`~compas_singular.mcp.protocol.Dispatcher` asks for content."""

    def __init__(self, session=None):
        self.session = session if session is not None else MeshSession()

    # -- tools ----------------------------------------------------------

    def list_tools(self):
        return registry.schemas()

    def call_tool(self, name, arguments):
        return registry.call(self.session, name, arguments)

    # -- resources ------------------------------------------------------

    def list_resources(self):
        resources = list(library.list_resources())
        # The live journal is not a file, so the library cannot list it.
        resources.insert(0, {
            'uri': 'session://current',
            'name': 'this session',
            'description': 'What has been done so far: state, history, remarks '
                           'and the current reading.',
            'mimeType': 'text/markdown',
        })
        return resources

    def read_resource(self, uri):
        if uri == 'session://current':
            return self.render_session()
        return library.read_resource(uri)

    # -- prompts --------------------------------------------------------

    def list_prompts(self):
        return library.list_prompts()

    def get_prompt(self, name, arguments):
        return library.get_prompt(name, arguments)

    # -- the live journal -----------------------------------------------

    def render_session(self):
        """The session as Markdown, so it survives a context compaction.

        Deliberately the same shape as a stored example, so what a model reads
        about its own run and what it reads about a past one look alike.
        """
        session = self.session
        state = session.state()
        lines = ['# This session', '']
        if not session.loaded:
            lines.append('Nothing is loaded. Call rhino_pull to take a mesh '
                         'out of the open Rhino document.')
            return '\n'.join(lines)

        source = session.source or {}
        lines.append('Mesh: {} faces, {} vertices, from {}.'.format(
            state['faces'], state['vertices'],
            source.get('layer') or source.get('name') or source.get('kind')
            or 'an unknown source'))
        metrics = session.quality()
        told = describe(session.mesh, metrics, library.thresholds())
        lines += ['', '**Now:** {} -- {}'.format(told['verdict'],
                                                 told['reading'])]

        if session.history:
            lines += ['', '## Steps', '']
            marks = {True: 'improved all three',
                     False: 'did NOT improve all three',
                     None: 'no change'}
            for entry in session.history:
                lines.append('{}. `{}` -- {}'.format(
                    entry['step'] + 1, entry['action'],
                    marks.get(entry.get('all_improved'), '')))
                for remark in session.remarks_for(entry['step']):
                    lines.append('   - {}'.format(remark))
        loose = [r for r in session.remarks if r['step'] >= len(session.history)]
        if loose:
            lines += ['', '## Remarks', '']
            lines += ['- {}'.format(r['text']) for r in loose]
        return '\n'.join(lines)


def build(session=None):
    """A dispatcher wired to a handler. The one thing ``__main__`` needs."""
    handler = Handler(session)
    if not handler.list_tools():
        raise RuntimeError(
            'no tools are registered -- the tools_* imports in server.py are '
            'what register them, so one has been dropped')
    return Dispatcher(handler, name='compas_singular.mcp', version='0.1.0',
                      instructions=INSTRUCTIONS)
