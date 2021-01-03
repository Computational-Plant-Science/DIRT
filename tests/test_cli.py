from os.path import join, splitext
from pathlib import Path

import pytest
from click.testing import CliRunner

from dirt import cli

test_dir = Path(__file__).parent
image_files = [
    'cassava1.jpg',
    'cassava2.jpg',
    'hemp1.jpg',
    'hemp2.jpg',
    'maize1.jpg',
    'maize2.jpg',
    'tomato1.jpg',
    'tomato2.jpg',
]


@pytest.mark.parametrize("test_image", image_files)
def test_cli(test_image):
    result = CliRunner().invoke(
        cli,
        [join(test_dir, 'data', test_image), '--output_directory', join(test_dir, 'output')])
    assert result.exit_code == 0
