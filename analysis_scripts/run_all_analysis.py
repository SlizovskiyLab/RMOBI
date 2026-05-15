#!/usr/bin/env python3
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


ANALYSIS_DIR = Path(__file__).resolve().parent
ROOT_DIR = ANALYSIS_DIR.parent

SCRIPT_ORDER = [
    "1_upset.py",
    "2_transitionRate.py",
    "3_grouped_bar_transition_by_mge.py",
    "4_nested_donut_transition_mge.py",
    "5_stacked_colocalization_status_top10.py",
    "6_pie_charts_resistance_groups_graph2.py",
    "7_heatmaps_graph2.py",
]

RUN_FROM_ROOT: set[str] = set()


def ensure_output_dirs() -> None:
    for rel_dir in [
        "images/upset",
        "images/transition_rates",
        "images/grouped_bar",
        "images/nested_donut",
        "images/stacked",
        "images/piecharts",
        "images/heatmaps",
    ]:
        (ANALYSIS_DIR / rel_dir).mkdir(parents=True, exist_ok=True)


def run_script(script_name: str, python_exe: str) -> int:
    script_path = ANALYSIS_DIR / script_name
    if not script_path.exists():
        print(f"[ERROR] Missing script: {script_path}")
        return 1

    script_cwd = ROOT_DIR if script_name in RUN_FROM_ROOT else ANALYSIS_DIR
    command = [python_exe, f"analysis_scripts/{script_name}"] if script_cwd == ROOT_DIR else [python_exe, script_name]

    print(f"\n[RUNNING] {script_name} (cwd: {script_cwd})")
    result = subprocess.run(command, cwd=str(script_cwd))
    if result.returncode == 0:
        print(f"[OK] {script_name}")
    else:
        print(f"[FAILED] {script_name} (exit code {result.returncode})")
    return result.returncode


def main() -> int:
    parser = argparse.ArgumentParser(description="Run all analysis scripts in sequence.")
    parser.add_argument(
        "--continue-on-error",
        action="store_true",
        help="Continue running remaining scripts after a failure.",
    )
    parser.add_argument(
        "--python",
        default=sys.executable,
        help="Python executable to use (default: current interpreter).",
    )
    args = parser.parse_args()

    ensure_output_dirs()

    failed = []
    for script_name in SCRIPT_ORDER:
        code = run_script(script_name, args.python)
        if code != 0:
            failed.append((script_name, code))
            if not args.continue_on_error:
                break

    print("\n" + "=" * 60)
    if not failed:
        print("All analysis scripts completed successfully.")
        print("Outputs are available in analysis_scripts/images/.")
        return 0

    print("Completed with failures:")
    for script_name, code in failed:
        print(f"- {script_name}: exit code {code}")
    print("Some outputs may still have been generated in analysis_scripts/images/.")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
