"""Test the `curp cal-tc` command
"""

import pytest
from base_test_class import BaseTestClass

class TestCurpCalTc(BaseTestClass):
    def test_curp_cal_tc(self):
        pass
    
    def _assert_output_files(self):
        """Assert that the output files are created and have the correct content
        
        Files to assert:
        - time integrated autocorrelation function values in netcdf4 format
        - overall conductivity
        """
        