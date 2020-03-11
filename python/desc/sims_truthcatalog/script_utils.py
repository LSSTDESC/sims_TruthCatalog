"""
Utilities for top-level scripts
"""
from datetime import datetime as dt
import sys

__all__ = ['print_callinfo', 'TIME_TO_SECOND_FMT']

TIME_TO_SECOND_FMT = '%Y-%m-%d %H:%M:%S'

def print_callinfo(prog, args):
    """
    Print information about how a script using argparse was called

    Parameters
    ----------
    prog   program name, typically sys.argv[0]
    args   object returned by ArgumentParser.parse_args()
    """

    print('{}   {}  invoked  with arguments'.format(dt.now().strftime(TIME_TO_SECOND_FMT), prog))
    for e in dir(args):
        if not e.startswith('_'):
            nm = 'args.' + e
            print('{}: {}'.format(e, eval(nm)))

    sys.stdout.flush()

