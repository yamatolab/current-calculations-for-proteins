"""Base class for all pytest classes
"""

import pytest
import os

BASE_OUTPUT_DIR = "test/pytests/output"
BASE_TRUE_OUTPUT_DIR = "test/pytests/true_output"
BASE_INPUT_DIR = "test/pytests/input"

class BaseTestClass:
    def __init__(self, *args, **kwargs):
        self.base_output_dir = kwargs.get("base_output_dir", BASE_OUTPUT_DIR)
        self.base_true_output_dir = kwargs.get("base_true_output_dir", BASE_TRUE_OUTPUT_DIR)
        self.base_input_dir = kwargs.get("base_input_dir", BASE_INPUT_DIR)
    
    def clear_output_files(self, output_dir: str = None, extensions: list[str] = [".nc", ".log", ".out", ".dat"]):
        """Clear the output files in the output directory
        
        Args:
            output_dir (str): The output directory to clear
            extensions (list[str]): The extensions of the files to clear
        """
        if output_dir is None:
            output_dir = self.base_output_dir
        
        for file in os.listdir(output_dir):
            if any(file.endswith(ext) for ext in extensions):
                os.remove(os.path.join(output_dir, file))
