#!/usr/bin/env python3
"""Run genus host-specificity analysis with latitude-blocked sampling-unit label shuffling."""

from latitude_blocks_common import run_workflow


if __name__ == "__main__":
    raise SystemExit(run_workflow("genus"))
