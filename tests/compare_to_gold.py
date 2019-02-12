"""Compares a test result to the gold standard.

Error codes:
----------------------------------------------------------
  G     | Error in the gold file -- bad testname?
  D     | Error in my file -- test did not complete?
  E     | Error in my file size -- test did not complete.
  *     | Unknown error
  F     | Failing test
  .     | Passing test
----------------------------------------------------------

"""

import argparse
import numpy as np
import os

TOL = 1.e-5

def run_comparison(testname):
    fname_gold = os.path.join('..', 'tests_gold', '%s.stdout'%testname)
    try:
        gold = np.loadtxt(fname_gold, skiprows=1)
    except IOError:
        print 'ERROR: cannot find gold file "%s", bad testname?'%fname_gold
        return 'G'

    fname_mine = '%s.stdout'%testname
    try:
        mine = np.loadtxt(fname_mine, skiprows=1)
    except IOError:
        print 'ERROR: cannot find my file "%s", failed test run?'%fname_mine
        return 'D'

    try:
        close = np.allclose(gold, mine)
    except ValueError:
        return 'E'

    if close:
        return '.'
    else:
        return 'F'
                      

if __name__ == "__main__":
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('testnames', nargs='+',
                      help='Name of tests to run')
    args = parser.parse_args()

    results = []
    print 'Comparing tests:',
    for test in args.testnames:
        try:
            results.append(run_comparison(test))
        except Exception as err:
            print err
            results.append('*')
    print ''.join(results)
