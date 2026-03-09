#!/usr/bin/env python3
"""Config-driven OpenCV tracker for image sequences.

Workflow:
1. Load image paths from a glob pattern.
2. Select a tracking ROI on the first frame.
3. Select crop ROI size (or use fixed size from config).
4. Track object across frames, crop, optionally rotate, save the result.
5. Estimate frame rate from filesystem timestamps (mtimes) and stamp rows with `T = frame_index / frame_rate`.
6. Save per-frame metadata and the frame-rate summary to disk.
"""

from __future__ import annotations

import argparse
import csv
import glob
import json
import os
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import cv2
import yaml


BBox = Tuple[int, int, int, int]

ROTATION_MODES = {
    "none": None,
    "cw": cv2.ROTATE_90_CLOCKWISE,
    "ccw": cv2.ROTATE_90_COUNTERCLOCKWISE,
}


DEFAULT_IMAGE_EXTENSIONS = [".tif", ".tiff", ".bmp", ".png", ".jpg", ".jpeg"]


TRACKER_FACTORY_CANDIDATES: Dict[str, Sequence[str]] = {
    "csrt": ["TrackerCSRT_create", "legacy.TrackerCSRT_create"],
    "kcf": ["TrackerKCF_create", "legacy.TrackerKCF_create"],
    "mosse": ["TrackerMOSSE_create", "legacy.TrackerMOSSE_create"],
}


def load_config(config_path: Path) -> Dict[str, Any]:
    with config_path.open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle)

    if not isinstance(config, dict):
        raise ValueError("Config root must be a YAML mapping.")

    required = ["input", "tracker", "crop", "output"]
    missing = [key for key in required if key not in config]
    if missing:
        raise ValueError(f"Missing required config sections: {missing}")

    return config


def resolve_path(path_like: str, config_dir: Path) -> Path:
    path = Path(path_like)
    if path.is_absolute():
        return path
    return (config_dir / path).resolve()


def resolve_pattern(pattern: str, config_dir: Path) -> str:
    if os.path.isabs(pattern):
        return pattern
    return str(config_dir / pattern)


def expand_directory_pattern(pattern: str) -> str:
    if not glob.has_magic(pattern) and os.path.isdir(pattern):
        return os.path.join(pattern, "*")
    return pattern


def create_tracker(tracker_name: str):
    name = tracker_name.strip().lower()
    allowed = list(TRACKER_FACTORY_CANDIDATES.keys())
    if name not in TRACKER_FACTORY_CANDIDATES:
        raise ValueError(
            f"Unsupported tracker '{tracker_name}'. Supported trackers: {allowed}"
        )

    for symbol in TRACKER_FACTORY_CANDIDATES[name]:
        obj: Any = cv2
        ok = True
        for part in symbol.split("."):
            if not hasattr(obj, part):
                ok = False
                break
            obj = getattr(obj, part)
        if ok:
            return obj()

    if hasattr(cv2, "Tracker_create"):
        return cv2.Tracker_create(tracker_name.upper())

    raise RuntimeError(
        f"Tracker '{tracker_name}' is not available in OpenCV {cv2.__version__}."
    )


def select_roi(window_name: str, frame, from_center: bool = False) -> BBox:
    bbox = cv2.selectROI(window_name, frame, fromCenter=from_center, showCrosshair=True)
    x, y, w, h = [int(v) for v in bbox]
    if w <= 0 or h <= 0:
        raise ValueError(f"No ROI selected for '{window_name}'.")
    return x, y, w, h


def clamp_crop(center_x: float, center_y: float, crop_w: int, crop_h: int, img_w: int, img_h: int) -> BBox:
    half_w = crop_w // 2
    half_h = crop_h // 2
    x0 = int(round(center_x)) - half_w
    y0 = int(round(center_y)) - half_h
    x0 = max(0, min(x0, img_w - crop_w))
    y0 = max(0, min(y0, img_h - crop_h))
    x1 = min(img_w, x0 + crop_w)
    y1 = min(img_h, y0 + crop_h)
    return x0, y0, x1 - x0, y1 - y0


def apply_rotation(image, rotation_mode: str):
    mode = rotation_mode.strip().lower()
    if mode not in ROTATION_MODES:
        raise ValueError(
            f"rotation must be one of {list(ROTATION_MODES.keys())}; got '{rotation_mode}'"
        )
    rotate_flag = ROTATION_MODES[mode]
    return cv2.rotate(image, rotate_flag) if rotate_flag is not None else image


def normalize_glob_pattern(pattern: str, base_dir: Path) -> str:
    resolved = resolve_path(pattern, base_dir)
    return str(resolved)


def list_images(patterns: Sequence[str]) -> List[Path]:
    paths: List[Path] = []
    for pattern in patterns:
        matched = glob.glob(pattern, recursive=True)
        paths.extend(Path(p) for p in matched if Path(p).is_file())
    return sorted(paths)


def main() -> None:
    parser = argparse.ArgumentParser(description="Track and crop image sequence using OpenCV.")
    parser.add_argument(
        "--config",
        default="../configs/opencv_tracker_v3.yaml",
        help="Path to YAML configuration file.",
    )
    args = parser.parse_args()

    config_path = Path(args.config).resolve()
    config_dir = config_path.parent
    config = load_config(config_path)

    raw_patterns = config["input"]["image_glob"]
    if isinstance(raw_patterns, str):
        raw_patterns = [raw_patterns]
    image_patterns = [
        expand_directory_pattern(resolve_pattern(pattern, config_dir))
        for pattern in raw_patterns
    ]
    image_paths = list_images(image_patterns)
    if not image_paths:
        raise FileNotFoundError(f"No images found for patterns: {image_patterns}")
    allowed_exts_raw = config["input"].get("extensions", DEFAULT_IMAGE_EXTENSIONS)
    allowed_exts: List[str] = (
        [allowed_exts_raw]
        if isinstance(allowed_exts_raw, str)
        else list(allowed_exts_raw)
    )
    normalized_exts = {
        ext.lower() if ext.startswith(".") else f".{ext.lower()}" for ext in allowed_exts
    }
    filtered_paths = [p for p in image_paths if p.suffix.lower() in normalized_exts]
    if not filtered_paths:
        raise FileNotFoundError(
            f"No images with extensions {sorted(normalized_exts)} found for {image_patterns}"
        )
    image_paths = filtered_paths
    image_patterns_resolved = image_patterns

    first_frame = cv2.imread(str(image_paths[0]), cv2.IMREAD_COLOR)
    if first_frame is None:
        raise RuntimeError(f"Unable to read first image: {image_paths[0]}")

    tracker_name = config["tracker"].get("type", "csrt")
    allowed_override = config["tracker"].get("allowed_types")
    if allowed_override:
        normalized = [t.strip().lower() for t in allowed_override]
        invalid = [
            t for t in normalized if t not in TRACKER_FACTORY_CANDIDATES.keys()
        ]
        if invalid:
            raise ValueError(
                f"Tracker.allowed_types contains unsupported names: {invalid}"
            )
        allowed_list = normalized
    else:
        allowed_list = list(TRACKER_FACTORY_CANDIDATES.keys())

    if tracker_name.strip().lower() not in allowed_list:
        raise ValueError(
            f"{tracker_name} is not enabled. Choose from {allowed_list} or adjust tracker.allowed_types."
        )

    tracker = create_tracker(tracker_name)

    window_track = config["tracker"].get("selection_window", "Select Tracking ROI")
    tracking_bbox = select_roi(window_track, first_frame, from_center=False)
    tracker.init(first_frame, tracking_bbox)

    crop_mode = config["crop"].get("mode", "select_roi_size")
    crop_w = config["crop"].get("width")
    crop_h = config["crop"].get("height")
    if crop_mode == "select_roi_size":
        window_crop = config["crop"].get("selection_window", "Select Crop ROI Size")
        _, _, crop_w, crop_h = select_roi(window_crop, first_frame, from_center=False)
    elif crop_mode == "fixed_size":
        if not crop_w or not crop_h:
            raise ValueError("crop.width and crop.height must be set when mode=fixed_size.")
        crop_w = int(crop_w)
        crop_h = int(crop_h)
    elif crop_mode == "tracker_box":
        crop_w, crop_h = int(tracking_bbox[2]), int(tracking_bbox[3])
    else:
        raise ValueError("crop.mode must be one of: select_roi_size, fixed_size, tracker_box")

    input_dir = image_paths[0].parent
    default_output_root = input_dir.parent / f"{input_dir.name}_cropped"
    configured_root = config["output"].get("root_dir")
    if configured_root:
        output_root = resolve_path(configured_root, config_dir)
    else:
        output_root = default_output_root
    crops_dir = output_root / config["output"].get("crops_subdir", "registered_crops")
    crops_dir.mkdir(parents=True, exist_ok=True)
    metadata_csv = output_root / config["output"].get("metadata_csv", "tracking_metadata.csv")
    summary_json = output_root / config["output"].get("summary_json", "frame_rate_summary.json")

    show_preview = bool(config["output"].get("show_preview", True))
    preview_window = config["output"].get("preview_window", "Tracking Preview")
    keep_source_stem = bool(config["output"].get("keep_source_stem", True))
    save_ext = config["output"].get("save_extension", ".png")
    save_color = bool(config["output"].get("save_color", True))
    rotation_mode = config["output"].get("rotation", "ccw").strip().lower()
    if rotation_mode not in ROTATION_MODES:
        raise ValueError(
            f"Invalid rotation '{rotation_mode}'. Choose from {list(ROTATION_MODES.keys())}."
        )

    rows: List[Dict[str, Any]] = []
    file_times: List[float] = []

    for index, image_path in enumerate(image_paths):
        frame = cv2.imread(str(image_path), cv2.IMREAD_COLOR)
        if frame is None:
            raise RuntimeError(f"Unable to read image: {image_path}")
        img_h, img_w = frame.shape[:2]

        if index == 0:
            success = True
            bbox = tracking_bbox
        else:
            success, updated = tracker.update(frame)
            if not success:
                break
            bbox = tuple(int(v) for v in updated)

        if crop_mode == "tracker_box" and success:
            crop_w = max(1, int(bbox[2]))
            crop_h = max(1, int(bbox[3]))

        center_x = float(bbox[0] + bbox[2] / 2.0)
        center_y = float(bbox[1] + bbox[3] / 2.0)
        crop_x, crop_y, crop_w_eff, crop_h_eff = clamp_crop(
            center_x, center_y, int(crop_w), int(crop_h), img_w, img_h
        )
        crop = frame[crop_y : crop_y + crop_h_eff, crop_x : crop_x + crop_w_eff]
        if not save_color:
            crop = cv2.cvtColor(crop, cv2.COLOR_BGR2GRAY)
        crop = apply_rotation(crop, rotation_mode)

        if keep_source_stem:
            out_name = f"{image_path.stem}_crop{save_ext}"
        else:
            out_name = f"frame_{index:06d}{save_ext}"
        out_path = crops_dir / out_name
        cv2.imwrite(str(out_path), crop)
        file_stat = image_path.stat()
        file_times.append(file_stat.st_mtime)
        saved_h, saved_w = crop.shape[:2]

        rows.append(
            {
                "frame_index": index,
                "image_path": str(image_path),
                "crop_path": str(out_path),
                "tracking_success": int(success),
                "track_x": int(bbox[0]),
                "track_y": int(bbox[1]),
                "track_w": int(bbox[2]),
                "track_h": int(bbox[3]),
                "crop_x": crop_x,
                "crop_y": crop_y,
                "crop_w": crop_w_eff,
                "crop_h": crop_h_eff,
                "file_timestamp_seconds": file_stat.st_mtime,
                "file_size_bytes": file_stat.st_size,
                "saved_width": saved_w,
                "saved_height": saved_h,
                "rotation_applied": rotation_mode,
            }
        )

        if show_preview:
            view = frame.copy()
            tx, ty, tw, th = [int(v) for v in bbox]
            cv2.rectangle(view, (tx, ty), (tx + tw, ty + th), (0, 255, 0), 2)
            cv2.rectangle(
                view,
                (crop_x, crop_y),
                (crop_x + crop_w_eff, crop_y + crop_h_eff),
                (255, 0, 0),
                2,
            )
            cv2.putText(
                view,
                f"{tracker_name.upper()} | frame={index} | ok={int(success)}",
                (10, 25),
                cv2.FONT_HERSHEY_SIMPLEX,
                0.7,
                (0, 255, 255),
                2,
            )
            cv2.imshow(preview_window, view)
            key = cv2.waitKey(1) & 0xFF
            if key in (27, ord("q")):
                break

    cv2.destroyAllWindows()

    frame_rate: Optional[float] = None
    sequence_duration: Optional[float] = None
    if len(file_times) >= 2:
        sequence_duration = file_times[-1] - file_times[0]
        if sequence_duration > 0.0 and len(rows) > 1:
            frame_rate = (len(rows) - 1) / sequence_duration

    delta_t = 1.0 / frame_rate if frame_rate else 0.0
    for idx, row in enumerate(rows):
        row["time_seconds"] = idx * delta_t if frame_rate else 0.0

    if rows:
        fieldnames = list(rows[0].keys())
        with metadata_csv.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            for row in rows:
                writer.writerow(row)

    summary = {
        "input_pattern": image_patterns_resolved,
        "allowed_extensions": sorted(normalized_exts),
        "num_input_images": len(image_paths),
        "num_saved_crops": len(rows),
        "tracker_type": tracker_name,
        "rotation_mode": rotation_mode,
        "frame_rate_hz": frame_rate,
        "frame_duration_seconds": 1.0 / frame_rate if frame_rate else None,
        "sequence_duration_seconds": sequence_duration,
        "metadata_csv": str(metadata_csv),
        "crops_dir": str(crops_dir),
        "output_root": str(output_root),
    }
    with summary_json.open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
