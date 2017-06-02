import os
import pytest
from xyalign import utils


def teardown_module(function):
	if os.path.exists(os.path.join(".", "temp_test")):
		os.rmdir(os.path.join(".", "temp_test"))


def test_validate_external_prog():
	assert utils.validate_external_prog("echo", "echo") == 0
	with pytest.raises(OSError):
		utils.validate_external_prog(
			"this_program_should_not_exist", "this_program_should_not_exist")


def test_validate_dir_exist():
	assert utils.validate_dir(".", "temp_test") is False
	assert utils.validate_dir(".", "temp_test") is True


# def test_chromosome_bed():
# 	with pytest.raises(RuntimeError):
# 		utils.chromosome_bed("foo", "foo.txt", ["foo"])
