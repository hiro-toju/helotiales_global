#!/usr/bin/env python3
"""Run SH host-specificity analysis with continent-blocked sampling-unit label shuffling."""

from spatial_blocks_common import run_workflow


if __name__ == "__main__":
    raise SystemExit(run_workflow("sh"))
