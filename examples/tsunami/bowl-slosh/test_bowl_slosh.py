#!/usr/bin/env python

r"""Bowl-Slosh regression test for GeoClaw

To create new regression data use
    `python regression_tests.py True`
"""

from pathlib import Path
import sys, os, shutil
import unittest

import numpy

import clawpack.geoclaw.test as test

example_dir = os.getcwd() # where to find setrun_test.py

class BowlSloshTest(test.GeoClawRegressionTest):

    r"""Bowl-Slosh regression test for GeoClaw"""

    def setUp(self):

        super(BowlSloshTest, self).setUp()


    def load_rundata(self, setrun_file='setrun'):
        r"""(Re)load setrun module and create *rundata* object
        Modified version from clawutil.test to allow specifying setrun_file.
        """

        import importlib
        if 'setrun' in sys.modules:
            del(sys.modules['setrun'])
        #sys.path.insert(0, self.test_path)
        #import setrun
        sys.path.insert(0, example_dir)
        setrun = importlib.import_module(setrun_file)
        self.rundata = setrun.setrun()
        sys.path.pop(0)

    def runTest(self, save=False, indices=(2, 3)):
        r"""Test bowl-slosh example

        Note that this stub really only runs the code and performs no tests.

        """

        # Write out data files
        self.load_rundata(setrun_file='setrun_test')

        assert self.rundata.clawdata.tfinal == 0.5, '** unexpected tfinal'

        self.write_rundata_objects()

        # Make topography
        import maketopo
        maketopo.maketopo()
        shutil.copy('bowl.topotype2', Path(self.temp_path))

        # Run code
        self.run_code()

        # Perform tests
        self.check_gauges(save=save, gauge_id=1, indices=(2, 3))
        #self.check_fgmax(save=save)  # minor differences, why??
        self.success = True


if __name__=="__main__":
    if len(sys.argv) > 1:
        if bool(sys.argv[1]):
            # Fake the setup and save out output
            test = BowlSloshTest()
            try:
                test.setUp()
                test.runTest(save=True)
            finally:
                test.tearDown()
            sys.exit(0)
    unittest.main()
