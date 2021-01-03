from os.path import join, splitext
from pathlib import Path

import pytest
from click.testing import CliRunner

from dirt import cli

test_dir = Path(__file__).parent
test_files = [
    'cassava1.jpg',
    'cassava2.jpg',
    'hemp1.jpg',
    'hemp2.jpg',
    'maize1.jpg',
    'maize2.jpg',
    'tomato1.jpg',
    'tomato2.jpg',
]


@pytest.mark.parametrize("test_file", test_files)
def test_cli(test_file):
    runner = CliRunner()
    # Path(join(test_dir, 'output')).mkdir(parents=True, exist_ok=True)
    result = runner.invoke(
        cli,
        [join(test_dir, 'data', test_file), '--output_directory', join(test_dir, 'output')])
    assert result.exit_code == 0
