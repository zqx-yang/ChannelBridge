# pymol -cq run_pse.py -- complex=<pdb> out=<pse> [tunnel_up=..] [tunnel_dn=..] [up_chains=B,C] [dn_chains=N] [ligands=NAP,BGC] [oic_arm=1.0] [label=..]
#      ... -- --selftest
#
# Thin PyMOL entry point for channelbridge.render.scenes.
#
# Two PyMOL facts this works around, both of which cost a run to discover:
#  * PyMOL's -m is not Python's -m; it runs a FILE, so this has to be a file.
#  * A file PyMOL execfile()s gets __file__ = pymol/__init__.py, NOT the
#    script. sys.argv[0] is the script path, so the package is located from
#    that instead -- no PYTHONPATH needed, and no second install into the
#    read-only pymol environment.
import os
import sys

_here = os.path.dirname(os.path.abspath(sys.argv[0]))
_src = os.path.dirname(os.path.dirname(_here))          # .../src
if os.path.isdir(os.path.join(_src, "channelbridge")) and _src not in sys.path:
    sys.path.insert(0, _src)

from channelbridge.render.scenes import main  # noqa: E402

# Called unconditionally: a file PyMOL execfile()s does NOT get
# __name__ == "__main__", so an `if __name__` guard here would silently do
# nothing and exit 0 -- which it did, once.
_rc = main()
if _rc:
    print("[run_pse] exit %d" % _rc, file=sys.stderr)
