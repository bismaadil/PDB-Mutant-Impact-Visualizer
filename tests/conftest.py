from pathlib import Path

import pytest

from pdb_mutant_impact_visualizer.fetch_structure import load_structure

FIXTURE = Path(__file__).parent / "fixtures" / "tiny_test.pdb"


@pytest.fixture
def fixture_path():
    return FIXTURE


@pytest.fixture
def structure():
    return load_structure(FIXTURE, "TEST")
