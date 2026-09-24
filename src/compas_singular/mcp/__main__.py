"""**The entry point: ``python -m compas_singular.mcp``.**

Speaks MCP on stdin and stdout and nothing else. Every diagnostic goes to
stderr, where an MCP client shows it as server log output; anything written to
stdout that is not a protocol frame breaks the session, which is why
:func:`~compas_singular.mcp.protocol.serve` rebinds ``sys.stdout`` before it
reads the first line.

Run it by hand to check the wiring::

    echo {"jsonrpc":"2.0","id":1,"method":"tools/list"} | python -m compas_singular.mcp
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import annotations

import sys
from typing import Sequence

from compas_singular.mcp.protocol import serve
from compas_singular.mcp.server import build


def main(argv: Sequence[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if '--version' in argv:
        sys.stderr.write('compas_singular.mcp 0.1.0\n')
        return 0
    if '--help' in argv or '-h' in argv:
        sys.stderr.write((__doc__ or '') + '\n')
        return 0
    try:
        dispatcher = build()
    except Exception as exc:
        sys.stderr.write('compas_singular.mcp failed to start: {}: {}\n'.format(
            type(exc).__name__, exc))
        return 1
    sys.stderr.write('compas_singular.mcp ready: {} tools\n'.format(
        len(dispatcher.handler.list_tools())))
    serve(dispatcher)
    return 0


if __name__ == '__main__':
    sys.exit(main())
