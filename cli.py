"""Backward-compatible entry point: ``python cli.py --pdb 1AKI --chain A --mutation S50V``."""

from pdb_mutant_impact_visualizer.commands import run_command as main

if __name__ == "__main__":
    main()
