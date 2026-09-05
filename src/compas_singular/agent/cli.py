"""**Run the agent from a terminal.** ``python -m compas_singular.agent.cli``

The headless entry point. Same session, same tools, same loop the Rhino command
drives -- the only difference is where the boundary comes from and where the
result goes.

    python -m compas_singular.agent.cli --demo plate-with-hole \\
        "even out the elements around the hole"

    python -m compas_singular.agent.cli --boundary site.json --out ./run \\
        --yes "mesh this at about 0.5m and smooth it"

**A dry run costs nothing and calls no model.** ``--dry-run`` swaps in
:class:`~.runner.backend.ScriptedBackend` and executes a fixed plan -- solve,
layout, densify, smooth -- so a boundary file, an install, or a set of
parameters can be checked before any tokens are spent. It is also what to reach
for when something fails and it is not obvious whether the model or the geometry
is at fault.

**Destructive tools are refused unless** ``--yes``. Four of the tools discard
work that cannot be recovered; without the flag they come back to the model as
refusals with an explanation, and the run continues without them. That is a
deliberate default for an unattended terminal.

Boundary files are JSON::

    {"outer":  [[0,0,0], [10,0,0], [10,10,0], [0,10,0]],
     "inners": [[[4,4,0], [6,4,0], [6,6,0], [4,6,0]]],
     "guides": [[[0,5,0], [10,5,0]]]}

``inners`` and ``guides`` are optional, and a closing point that repeats the
first is accepted but not required.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import argparse
import json
import os
import sys
from math import cos
from math import pi
from math import sin


__all__ = ['main', 'load_boundary', 'demo_boundary', 'DEMOS']


def _circle(cx, cy, r, n=32):
    return [[cx + r * cos(2 * pi * i / n), cy + r * sin(2 * pi * i / n), 0.0]
            for i in range(n)]


#: Shapes to try the thing on without drawing anything. Each returns
#: ``(outer, inners, guides)``.
DEMOS = {
    'square': lambda: ([[0.0, 0.0, 0.0], [10.0, 0.0, 0.0],
                        [10.0, 10.0, 0.0], [0.0, 10.0, 0.0]], [], []),
    'plate-with-hole': lambda: ([[0.0, 0.0, 0.0], [10.0, 0.0, 0.0],
                                 [10.0, 10.0, 0.0], [0.0, 10.0, 0.0]],
                                [_circle(5.0, 5.0, 2.0)], []),
    'two-holes': lambda: ([[0.0, 0.0, 0.0], [16.0, 0.0, 0.0],
                           [16.0, 10.0, 0.0], [0.0, 10.0, 0.0]],
                          [_circle(4.5, 5.0, 1.8), _circle(11.5, 5.0, 1.8)], []),
    'cable': lambda: ([[0.0, 0.0, 0.0], [12.0, 0.0, 0.0],
                       [12.0, 8.0, 0.0], [0.0, 8.0, 0.0]], [],
                      [[[0.5, 4.0, 0.0], [11.5, 4.0, 0.0]]]),
    'l-plate': lambda: ([[0.0, 0.0, 0.0], [12.0, 0.0, 0.0], [12.0, 5.0, 0.0],
                         [6.0, 5.0, 0.0], [6.0, 11.0, 0.0], [0.0, 11.0, 0.0]],
                        [], []),
}


def demo_boundary(name):
    if name not in DEMOS:
        raise SystemExit('unknown demo {!r}. Available: {}'.format(
            name, ', '.join(sorted(DEMOS))))
    return DEMOS[name]()


def load_boundary(path):
    """``(outer, inners, guides)`` from a JSON file. See the module docstring."""
    with open(path) as handle:
        data = json.load(handle)
    if 'outer' not in data:
        raise SystemExit('{}: no "outer" key -- a boundary file needs at least '
                         'an outer loop'.format(path))

    def points(seq):
        out = []
        for p in seq:
            p = list(p)
            while len(p) < 3:
                p.append(0.0)
            out.append([float(p[0]), float(p[1]), float(p[2])])
        return out

    return (points(data['outer']),
            [points(loop) for loop in data.get('inners') or []],
            [points(curve) for curve in data.get('guides') or []])


#: What ``--dry-run`` executes. A whole pipeline, so the plan exercises the
#: field, the layout, both density knobs and the smoothing pass.
DRY_PLAN = [
    [{'name': 'check_inputs', 'arguments': {}}],
    [{'name': 'enter_coarse', 'arguments': {}}],
    [{'name': 'inspect', 'arguments': {}}],
    [{'name': 'list_strips', 'arguments': {}}],
    [{'name': 'rebuild_dense', 'arguments': {}}],
    [{'name': 'smooth', 'arguments': {}}],
    'Dry run: solved, laid out, densified and smoothed. No model was called.',
]


def build_parser():
    parser = argparse.ArgumentParser(
        prog='python -m compas_singular.agent.cli',
        description='Drive the compas_singular meshing pipeline with a model.')
    parser.add_argument('instruction', nargs='?', default=None,
                        help='What to do, in words.')

    source = parser.add_argument_group('the domain')
    source.add_argument('--boundary', metavar='FILE',
                        help='JSON file with outer/inners/guides.')
    source.add_argument('--demo', metavar='NAME',
                        help='A built-in shape instead: {}.'.format(
                            ', '.join(sorted(DEMOS))))
    source.add_argument('--spacing', type=float, default=None,
                        help='Background triangulation spacing for the field '
                             'solve. NOT the element size.')
    source.add_argument('--target-length', type=float, default=None,
                        help='Element size for densification. This IS the quad '
                             'size.')
    source.add_argument('--relax', action='store_true',
                        help='Diffusion solver. Right for OPEN guides; it can '
                             'collapse a CLOSED one.')

    run = parser.add_argument_group('the run')
    run.add_argument('--backend', default='anthropic',
                     choices=['anthropic', 'gemini'],
                     help='Which model provider. "gemini" needs google-genai '
                          'and GEMINI_API_KEY.')
    run.add_argument('--model', default=None,
                     help='Defaults to Claude Opus, or Gemini Flash under '
                          '--backend gemini.')
    run.add_argument('--key-file', metavar='PATH', default=None,
                     help='Text file holding the API key, one line. Keep it '
                          'OUTSIDE this repo. Without it the environment is '
                          'used.')
    run.add_argument('--effort', default='xhigh',
                     choices=['low', 'medium', 'high', 'xhigh', 'max'])
    run.add_argument('--max-steps', type=int, default=40)
    run.add_argument('--images', default='transitions',
                     choices=['transitions', 'always', 'never'])
    run.add_argument('--yes', action='store_true',
                     help='Allow the four tools that destroy work. Without '
                          'this they are refused and the run continues.')
    run.add_argument('--ask', action='store_true',
                     help='Prompt on the terminal for each destructive tool.')
    run.add_argument('--dry-run', action='store_true',
                     help='Run a fixed plan with no model. Costs nothing.')
    run.add_argument('--no-cache', action='store_true',
                     help='Turn off prompt caching. For a determinism check.')

    out = parser.add_argument_group('output')
    out.add_argument('--out', metavar='DIR',
                     help='Write mesh JSON, a PNG, an SVG and the transcript.')
    out.add_argument('--quiet', action='store_true')
    return parser


def _ask(name, arguments):
    """Terminal confirmation for one destructive tool."""
    print('\n  {} destroys work that cannot be recovered.'.format(name))
    if arguments:
        print('  arguments: {}'.format(arguments))
    try:
        answer = raw_input('  allow it? [y/N] ')      # noqa: F821  (py2)
    except NameError:
        answer = input('  allow it? [y/N] ')
    return answer.strip().lower().startswith('y')


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)

    if not args.boundary and not args.demo:
        parser.error('give --boundary FILE or --demo NAME')
    if args.boundary and args.demo:
        parser.error('--boundary and --demo are alternatives')
    if not args.instruction and not args.dry_run:
        parser.error('an instruction is required (or use --dry-run)')

    # Imported here rather than at module scope so ``--help`` works instantly
    # and without solving anything.
    from ..framefield import cache as solve_cache
    from .core import MeshEditSession
    from .core.render import render_png, render_svg
    from .runner import ScriptedBackend, run_session

    outer, inners, guides = (load_boundary(args.boundary) if args.boundary
                             else demo_boundary(args.demo))

    def say(message):
        if not args.quiet:
            print(message)

    say('domain: outer {} points, {} hole(s), {} guide(s)'.format(
        len(outer), len(inners), len(guides)))

    solve_kwargs = {}
    if args.spacing is not None:
        solve_kwargs['target_length'] = args.spacing
    if args.relax:
        solve_kwargs['relax'] = True

    # The backend FIRST, before the field solve. A missing SDK or an unreadable
    # key file should cost nothing, and the solve is seconds of work that would
    # be thrown away -- measured on a plate with a hole, the run reached
    # "solving the field..." before noticing the key path was wrong.
    if args.dry_run:
        backend = ScriptedBackend(script=list(DRY_PLAN))
        say('dry run: no model will be called')
    elif args.backend == 'gemini':
        from .runner import GeminiBackend
        from .runner import KeyNotFound
        try:
            backend = GeminiBackend(model=args.model, key_file=args.key_file)
        except (ImportError, KeyNotFound) as exc:
            raise SystemExit(str(exc))
        # --effort and --no-cache are Anthropic-only. Say so rather than
        # letting a flag look as though it took effect.
        for flag, given in (('--effort', args.effort != 'xhigh'),
                            ('--no-cache', args.no_cache)):
            if given:
                say('note: {} has no effect on the gemini backend'.format(flag))
        say('model: {} (gemini)'.format(backend.model))
        if args.images == 'never':
            say('note: --images never means the model cannot see its own work. '
                'It will judge the mesh only from the numbers.')
    else:
        from .runner import AnthropicBackend
        from .runner import KeyNotFound
        try:
            backend = AnthropicBackend(model=args.model, effort=args.effort,
                                       cache=not args.no_cache,
                                       key_file=args.key_file)
        except (ImportError, KeyNotFound) as exc:
            raise SystemExit(str(exc))
        say('model: {}, effort {}'.format(backend.model, backend.effort))

    confirm = True if args.yes else (_ask if args.ask else False)
    if confirm is False:
        # Said on a dry run too. The setting is genuinely in force either way,
        # and a user who dry-runs first should see the same notice they will see
        # live -- a gate that appears only in one of the two modes is worse than
        # a line of noise.
        say('note: the four destructive tools will be refused. Pass --yes to '
            'allow them, or --ask to be prompted.')

    # Only now, with a usable backend in hand, is the solve worth paying for.
    say('solving the field...')
    decomposition = solve_cache.solve(outer, inners or None,
                                      guides=guides or None, build=True,
                                      **solve_kwargs)

    session = MeshEditSession(decomposition, outer=outer, inners=inners,
                              guides=guides, solve_params=solve_kwargs,
                              target_length=args.target_length)

    def on_event(event):
        if args.quiet:
            return
        if event['type'] == 'tool':
            print('  {} {}'.format('.' if event.get('ok') else 'x',
                                   event['message']))
        elif event['type'] == 'refused':
            print('  ! {}'.format(event['message']))
        elif event['type'] == 'text':
            print('\n{}\n'.format(event['message']))

    try:
        result = run_session(session, args.instruction or 'run the pipeline',
                             backend=backend, confirm=confirm,
                             max_steps=args.max_steps, images=args.images,
                             on_event=on_event)
    except Exception as exc:                                # noqa: BLE001
        # An auth failure is the overwhelmingly likely one on a first run, and
        # neither SDK's message names the environment variable to set.
        blob = '{}: {}'.format(type(exc).__name__, exc)
        looks_like_auth = any(marker in blob for marker in (
            'Authentication', 'authentication_error', 'API key', 'api_key',
            'PERMISSION_DENIED', 'UNAUTHENTICATED', '401', '403'))
        if looks_like_auth:
            how = ('Set ANTHROPIC_API_KEY, or run "ant auth login" -- the '
                   'client reads either.' if args.backend == 'anthropic' else
                   'Set GEMINI_API_KEY (or GOOGLE_API_KEY) to a Google AI '
                   'Studio key.')
            raise SystemExit(
                'the API rejected the credentials.\n  {}\n  {}\n'
                '  To check everything else first, re-run with --dry-run, '
                'which calls no model.'.format(how, blob[:160]))
        raise

    say('\n{}'.format(result.summary()))
    for entry in result.calls(ok=False):
        say('  refused: {} -- {}'.format(entry['name'], entry['reason'][:100]))

    quality = (result.digest.get('dense') or {}).get('quality')
    if quality:
        say('quality: min angle {:.1f}, max angle {:.1f}, aspect {:.2f}'.format(
            quality['min_angle'], quality['max_angle'], quality['aspect_max']))
        floor = (result.digest.get('dense') or {}).get('hard_floor') or {}
        if not floor.get('ok', True):
            say('  QUALITY FAILURE: {}'.format(floor.get('reason')))
    for warning in result.digest.get('warnings') or []:
        say('warning: {}'.format(warning))

    if args.out:
        if not os.path.isdir(args.out):
            os.makedirs(args.out)
        if result.mesh is not None:
            # ``to_json``, not ``json.dump(mesh.to_data())``: COMPAS 2 dropped
            # ``to_data`` and the JSON this writes carries the type information
            # ``compas.json_load`` needs to give the mesh back as a mesh.
            result.mesh.to_json(os.path.join(args.out, 'mesh.json'))
            render_png(session, path=os.path.join(args.out, 'mesh.png'))
            render_svg(session, path=os.path.join(args.out, 'mesh.svg'))
        with open(os.path.join(args.out, 'transcript.json'), 'w') as handle:
            json.dump({'instruction': args.instruction,
                       'summary': result.summary(),
                       'stage': result.stage,
                       'stop_reason': result.stop_reason,
                       'usage': result.usage,
                       'calls': result.transcript,
                       'text': result.text}, handle, indent=2, default=str)
        say('wrote {}'.format(args.out))

    if result.usage:
        say('tokens: {}'.format(', '.join(
            '{} {}'.format(k, v) for k, v in sorted(result.usage.items()))))

    return 0 if result.stop_reason in ('end_turn', 'max_steps') else 1


if __name__ == '__main__':
    sys.exit(main())
