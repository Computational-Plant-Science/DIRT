from os.path import join, splitext
from pathlib import Path

import pytest
from click.testing import CliRunner

from dirt import cli

test_dir = Path(__file__).parent
output_dir = join(test_dir, 'output')
Path(output_dir).mkdir(parents=True, exist_ok=True)
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
        [join(test_dir, 'data', test_image), '--output_directory', output_dir])
    assert result.exit_code == 0
