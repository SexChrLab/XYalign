import os
import pytest
import subprocess
from xyalign import xyalign

# Get current directory
dir = os.path.dirname(__file__)

# Test if "samtools" and "bwa" are available
try:
	a = subprocess.Popen(
		["samtools"], stderr=subprocess.PIPE, stdout=subprocess.PIPE)
	b = subprocess.Popen(
		["bwa"], stderr=subprocess.PIPE, stdout=subprocess.PIPE)
	tools_present = True
except OSError:
	tools_present = False


def teardown_module(function):
	teardown_files = []
	for file_name in teardown_files:
		if os.path.exists(os.path.join(dir, file_name)):
			os.remove(os.path.join(dir, file_name))


def test_xyalign():
	"""
	All this does is ensure xyalign.py is checked for typical Python syntax errors
	"""
	pass
