"""Test the `curp compute` command
"""

import pytest
from base_test_class import BaseTestClass

class TestCurpCompute(BaseTestClass):
    def test_curp_compute_no_decomp_intra_whole(self):
        pass
    
    def test_curp_compute_no_decomp_inter_residues(self):
        pass
    
    def test_curp_compute_decomp_intra_whole(self):
        pass
    
    def test_curp_compute_decomp_inter_residues(self):
        pass
    
    
    def _assert_output_files(self):
        """Assert that the output files are created and have the correct content
        
        Files to assert:
        - heat flux in netcdf4 format
        """
        # [TODO]
        # 1. assert that the output files are created
        # 2. assert that the output files have the correct content
        pass
        
        # after some assertions, clear the output files
        self.clear_output_files()