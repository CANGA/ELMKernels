"""Compares a test result to the gold standard.

Error codes:
----------------------------------------------------------
  G     | Error in the gold file -- bad testname?
  D     | Error in my file -- test did not complete?
  E     | Error in my file size -- test did not complete.
  N/A   | Unknown error
  FAIL  | Failing test
  PASS  | Passing test
----------------------------------------------------------

"""
from termcolor import colored 
import argparse
import numpy as np
import os

TOL = 1.e-5

def run_comparison(testname, full_message=False):
    # get gold file
    fname_gold = os.path.join('..', 'tests_gold', '%s.soln'%testname)
    try:
        gold = np.loadtxt(fname_gold, skiprows=1)
    except IOError:
        print('ERROR: cannot find gold file "%s", bad testname?'%fname_gold)
        return colored('Error in the gold file ', 'magenta')

    # get my file
    fname_mine = '%s.soln'%testname
    try:
        mine = np.loadtxt(fname_mine, skiprows=1)
    except IOError:
        fname_mine = '%s.stdout'%testname
        try:
            mine = np.loadtxt(fname_mine, skiprows=1)
        except IOError:
            print('ERROR: cannot find my file "%s", failed test run?'%fname_mine)
            return colored('Error in my file ', 'cyan')
        
    # compare
    try:
        close = np.allclose(gold, mine, rtol=1.e-10, atol=1.e-10)
    except ValueError as err:
        if full_message:
            print('')
            print('Test: "%s" FAILED with error: "%r"'%(testname, err))
        return colored('Error in my file size ','yellow')

    if close:
        return colored('PASS ','green')
    else:
        return colored('FAIL ','red')
                      

if __name__ == "__main__":
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('testnames', nargs='+',
                      help='Name of tests to run')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Print error messages')
    args = parser.parse_args()

    results = []
    print('Comparing tests:')
    for test in args.testnames:
        try:
            results.append(run_comparison(test, args.verbose))
        except Exception as err:
            print(err)
            results.append('N/A ')
    print(''.join(results))
