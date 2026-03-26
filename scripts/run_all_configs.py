#!/usr/bin/env python3
"""Run the tracker pipeline for every config under `configs/`."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import yaml


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Iterate over tracker configs and execute the tracker script."
    )
    repo_root = Path(__file__).resolve().parent.parent
    tracker_default = repo_root / "image_registration" / "opencv_tracker_v3.py"
    parser.add_argument(
        "--configs-dir",
        type=Path,
        default=repo_root / "configs",
        help="Directory containing YAML configs to run.",
    )
    parser.add_argument(
        "--tracker",
        type=Path,
        default=tracker_default,
        help="Tracker entry point (should accept `--config`).",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show which commands would execute without running the tracker.",
    )
    parser.add_argument(
        "--stop-on-error",
        action="store_true",
        help="Exit immediately if a tracker run fails.",
    )
    return parser.parse_args()


def find_configs(config_dir: Path) -> list[Path]:
    if not config_dir.exists():
        raise FileNotFoundError(f"configs directory not found: {config_dir}")
    return sorted(config_dir.glob("*.yaml"))


def describe_config(config_path: Path) -> tuple[str, str]:
    try:
        payload = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}
    except Exception:
        return ("unknown", "unknown")
    experiment = payload.get("input", {}).get("base_dir") or "unspecified"
    experiment = Path(experiment).name if experiment else "unspecified"
    subfolder = payload.get("input", {}).get("subfolder") or "all"
    return experiment, subfolder


def run_tracker(config_path: Path, tracker: Path, dry_run: bool) -> int:
    cmd = [sys.executable, str(tracker), "--config", str(config_path)]
    experiment, subfolder = describe_config(config_path)
    print(f"\nProcessing experiment '{experiment}', sub-folder '{subfolder}' (config: {config_path.name})")
    print("Command:", " ".join(cmd))
    if dry_run:
        return 0
    result = subprocess.run(cmd)
    return result.returncode


def main() -> None:
    args = parse_args()
    configs = find_configs(args.configs_dir)
    if not configs:
        raise RuntimeError(f"No configs found in {args.configs_dir}")

    for cfg in configs:
        rc = run_tracker(cfg, args.tracker, args.dry_run)
        if rc != 0:
            print(f"Tracker failed for {cfg.name} with exit code {rc}")
            if args.stop_on_error:
                sys.exit(rc)
    print("All configs processed.")


if __name__ == "__main__":
    main()
