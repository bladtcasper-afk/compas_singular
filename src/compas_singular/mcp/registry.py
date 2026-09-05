"""**The tool table: one decorator, one dispatch point, one place schemas live.**

A tool is a plain function taking a
:class:`~compas_singular.mcp.session.MeshSession` first and keyword arguments
after, returning ``{'ok': bool, ...}``. The JSON schema for its arguments sits on
the decorator, immediately above the function it describes, because a schema
kept anywhere else drifts from the code within a release or two.

**Every tool carries annotations, and they are not decoration.** MCP lets a
server declare whether a tool only reads, whether it destroys something, and
whether it reaches outside the process. A client uses those to decide what to
show a person before running it. Three of the tools here touch things the
session does not own -- ``rhino_push`` writes into a live CAD document,
``save_example`` writes into the library on disk -- and a client that cannot tell
those apart from ``inspect`` will either confirm everything or confirm nothing.

**A refusal is a return value, never an exception.** A tool that will not do
something answers ``{'ok': False, 'reason': ...}`` and the reason travels to the
model, which is the only way it can correct itself. Exceptions are reserved for
this server being broken.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function


__all__ = ['ToolSpec', 'TOOLS', 'tool', 'tool_names', 'schemas', 'call']


class ToolSpec(object):
    """One callable tool: the function, what it is for, and its arguments."""

    def __init__(self, name, function, description, schema,
                 read_only=False, destructive=False, idempotent=False,
                 open_world=False, title=None):
        self.name = name
        self.function = function
        self.description = description
        self.schema = schema
        self.read_only = read_only
        #: Removes or overwrites something a person would miss.
        self.destructive = destructive
        self.idempotent = idempotent
        #: Reaches beyond this process -- the Rhino document, or the library.
        self.open_world = open_world
        self.title = title or name

    def annotations(self):
        return {
            'title': self.title,
            'readOnlyHint': self.read_only,
            'destructiveHint': self.destructive,
            'idempotentHint': self.idempotent,
            'openWorldHint': self.open_world,
        }

    def __repr__(self):
        return '<ToolSpec {}>'.format(self.name)


TOOLS = {}


def tool(name, description, properties=None, required=(), read_only=False,
         destructive=False, idempotent=False, open_world=False, title=None):
    """Register a tool. ``properties`` is the JSON-schema body of its arguments."""
    def wrap(function):
        schema = {
            'type': 'object',
            'properties': dict(properties or {}),
            'required': list(required),
            'additionalProperties': False,
        }
        TOOLS[name] = ToolSpec(name, function, description, schema,
                               read_only=read_only, destructive=destructive,
                               idempotent=idempotent, open_world=open_world,
                               title=title)
        return function
    return wrap


def tool_names():
    return sorted(TOOLS)


def schemas():
    """``[{name, description, inputSchema, annotations}, ...]`` for ``tools/list``.

    Note ``inputSchema``, camelCase -- MCP's spelling, which is not the
    ``input_schema`` of Anthropic's own tool API.
    """
    return [{'name': spec.name,
             'description': spec.description,
             'inputSchema': spec.schema,
             'annotations': spec.annotations()}
            for spec in (TOOLS[name] for name in tool_names())]


def call(session, name, arguments=None):
    """Run a tool by name. An unknown name or bad arguments is a refusal.

    Parameters
    ----------
    session : MeshSession
    name : str
    arguments : dict, optional

    Returns
    -------
    dict
    """
    spec = TOOLS.get(name)
    if spec is None:
        return {'ok': False, 'reason': 'no tool named {!r}'.format(name),
                'available': tool_names()}
    arguments = dict(arguments or {})
    try:
        return spec.function(session, **arguments)
    except TypeError as exc:
        # Python names the IMPLEMENTATION in the message, which is not a name
        # the caller has ever seen. Swap the tool's own name in so the
        # correction is obvious, and list what it does accept.
        detail = str(exc).replace(spec.function.__name__, name)
        accepted = sorted(spec.schema['properties'])
        return {'ok': False,
                'reason': 'bad arguments for {}: {}'.format(name, detail),
                'accepted': accepted}
