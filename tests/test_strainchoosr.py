#!/usr/bin/env python

import pytest
from strainchoosr.strainchoosr import *


def test_read_weights_file_good():
    weight_dict = read_weights_file('tests/text_files/good_weights_file.txt')
    assert weight_dict['strain1'] == 2.3
    assert weight_dict['strain2'] == 3.6
    assert len(weight_dict) == 2


def test_weights_file_wrong_separator():
    with pytest.raises(RuntimeError):
        read_weights_file('tests/text_files/weights_wrong_separator.txt')


def test_weights_file_bad_second_column():
    with pytest.raises(ValueError):
        read_weights_file('tests/text_files/weights_wrong_second_column.txt')


def test_weights_file_with_blank_lines():
    weight_dict = read_weights_file('tests/text_files/good_weights_file_with_blank_line.txt')
    assert weight_dict['strain1'] == 2.3
    assert weight_dict['strain2'] == 3.6
    assert len(weight_dict) == 2
