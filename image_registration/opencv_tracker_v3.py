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
import io
import json
import os
import re
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


DEFAULT_CROPS_SUBDIR = "cropped"
DEFAULT_FONT_SCALE = 1.0


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


def choose_start_frame(
    image_paths: List[Path], slider_cfg: Dict[str, Any]
) -> Tuple[int, Optional[Any], bool]:
    window_name = slider_cfg.get("window_name", "Select Start Frame")
    trackbar_name = slider_cfg.get("trackbar_name", "Frame")
    instructions = slider_cfg.get("instructions", "Adjust slider and press 's'.")
    start_index = slider_cfg.get("start_index", 0)
    max_index = len(image_paths) - 1
    start_index = max(0, min(start_index, max_index))
    font_scale = slider_cfg.get("font_scale", DEFAULT_FONT_SCALE)
    font_thickness = max(1, round(font_scale * 2))

    selection = {"index": start_index, "frame": None}

    cv2.namedWindow(window_name, cv2.WINDOW_NORMAL)

    def on_trackbar(pos: int) -> None:
        pos = max(0, min(pos, max_index))
        selection["index"] = pos
        frame = cv2.imread(str(image_paths[pos]), cv2.IMREAD_COLOR)
        if frame is None:
            raise RuntimeError(f"Unable to read image for slider: {image_paths[pos]}")
        display = frame.copy()
        cv2.putText(
            display,
            f"frame {pos + 1}/{max_index + 1}",
            (10, 30),
            cv2.FONT_HERSHEY_SIMPLEX,
            font_scale,
            (0, 255, 255),
            font_thickness,
        )
        cv2.putText(
            display,
            instructions,
            (10, display.shape[0] - 10),
            cv2.FONT_HERSHEY_SIMPLEX,
            font_scale,
            (255, 255, 255),
            font_thickness,
        )
        cv2.imshow(window_name, display)
        selection["frame"] = frame

    cv2.createTrackbar(trackbar_name, window_name, start_index, max_index, on_trackbar)
    on_trackbar(start_index)

    while True:
        key = cv2.waitKey(100) & 0xFF
        if key in (ord("s"), ord("S")):
            break
        if key in (ord("x"), ord("X")):
            cv2.destroyWindow(window_name)
            return selection["index"], selection["frame"], True
        if key in (27, ord("q"), ord("Q")):
            cv2.destroyWindow(window_name)
            raise KeyboardInterrupt("Slider cancelled by user.")

    cv2.destroyWindow(window_name)
    return selection["index"], selection["frame"], False


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


def natural_sort_key(path: Path):
    parts = re.split(r"(\d+)", path.name)
    return [int(part) if part.isdigit() else part.lower() for part in parts]


def collect_image_paths(folder: Path, allowed_exts: Sequence[str]) -> List[Path]:
    filtered = [
        p
        for p in folder.iterdir()
        if p.is_file()
        and not p.name.startswith(".")
        and p.suffix.lower() in allowed_exts
    ]
    return sorted(filtered, key=natural_sort_key)


def format_slider_text(template: str, folder_name: str, default: str) -> str:
    if not template:
        return default
    try:
        return template.format(folder=folder_name)
    except Exception:
        return template


def measure_diameter_from_threshold(roi, cfg: Dict[str, Any]) -> Tuple[Optional[float], Optional[int]]:
    gray = cv2.cvtColor(roi, cv2.COLOR_BGR2GRAY) if roi.ndim == 3 else roi.copy()
    otsu_value, _ = cv2.threshold(gray, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
    otsu_scalar = float(otsu_value)
    default_threshold = cfg.get("default_threshold")
    if default_threshold is not None:
        threshold_value = int(default_threshold)
    else:
        threshold_value = int(round(otsu_scalar))
    min_val = int(cfg.get("min", 0))
    max_val = int(cfg.get("max", 255))
    window_name = cfg.get("window_name", "Threshold tuning")
    instructions = cfg.get(
        "instructions",
        "Adjust threshold and press ENTER/SPACE to accept or 'c' to cancel.",
    )
    threshold_value = max(min(threshold_value, max_val), min_val)
    font_scale = cfg.get("font_scale", DEFAULT_FONT_SCALE)
    font_thickness = max(1, round(font_scale * 2))

    confirmed = False

    def on_change(val):
        nonlocal threshold_value
        threshold_value = max(min_val, min(max_val, val))

    cv2.namedWindow(window_name, cv2.WINDOW_NORMAL)
    cv2.createTrackbar("Threshold", window_name, threshold_value, max_val, on_change)

    mask = None
    while True:
        _, mask = cv2.threshold(gray, threshold_value, 255, cv2.THRESH_BINARY)
        mask_overlay = cv2.cvtColor(mask, cv2.COLOR_GRAY2BGR)
        preview = cv2.addWeighted(
            cv2.cvtColor(gray, cv2.COLOR_GRAY2BGR), 0.6, mask_overlay, 0.4, 0
        )
        contours_data = cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        contours = contours_data[-2] if len(contours_data) == 3 else contours_data[0]
        if contours:
            cv2.drawContours(preview, contours, -1, (0, 0, 255), 1)
        cv2.putText(
            preview,
            f"Threshold {threshold_value}",
            (10, 30),
            cv2.FONT_HERSHEY_SIMPLEX,
            font_scale,
            (0, 255, 255),
            font_thickness,
        )
        cv2.putText(
            preview,
            instructions,
            (10, preview.shape[0] - 10),
            cv2.FONT_HERSHEY_SIMPLEX,
            font_scale,
            (255, 255, 255),
            font_thickness,
        )
        cv2.imshow(window_name, preview)
        key = cv2.waitKey(100) & 0xFF
        if key in (13, 32):
            confirmed = True
            break
        if key in (ord("c"), 27, ord("q")):
            break

    cv2.destroyWindow(window_name)

    if not confirmed:
        return None, None

    contours = cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)[0]
    if not contours:
        return None, threshold_value
    largest = max(contours, key=cv2.contourArea)
    (_, _), radius = cv2.minEnclosingCircle(largest)
    return radius * 2, threshold_value


def analyze_folder_metadata(folder: Path) -> Dict[str, Any]:
    joined = " ".join(folder.parts)
    normalized = re.sub(r"[\W_]+", " ", joined).lower()
    objective = "20x" if "20x" in normalized else "10x"
    pixelperum = 2.4 if objective == "20x" else 1.36
    particle_type = "pyrite" if "pyrite" in normalized else "plume particle"
    return {
        "objective": objective,
        "pixelperum": pixelperum,
        "particle_type": particle_type,
    }


def ensure_experiment_metadata(folder: Path) -> Path:
    json_path = folder / "metadata.json"
    csv_path = folder / "metadata.csv"
    if json_path.exists():
        return json_path
    if csv_path.exists():
        return csv_path
    data = analyze_folder_metadata(folder)
    json_path.write_text(json.dumps(data, indent=2), encoding="utf-8")
    return json_path


def measure_diameter(frame, cfg: Dict[str, Any]) -> Optional[float]:
    window_name = cfg.get("window_name", "Measure Diameter")
    instructions = cfg.get(
        "instructions",
        "Click and drag to draw a circle covering the object; release to finish.",
    )
    center = None
    radius = 0
    finished = False

    def on_mouse(event, x, y, flags, param):
        nonlocal center, radius, finished
        if event == cv2.EVENT_LBUTTONDOWN:
            center = (x, y)
            radius = 0
        elif event == cv2.EVENT_MOUSEMOVE and center is not None:
            dx = x - center[0]
            dy = y - center[1]
            radius = int((dx ** 2 + dy ** 2) ** 0.5)
        elif event == cv2.EVENT_LBUTTONUP and center is not None:
            dx = x - center[0]
            dy = y - center[1]
            radius = int((dx ** 2 + dy ** 2) ** 0.5)
            finished = True

    cv2.namedWindow(window_name, cv2.WINDOW_NORMAL)
    cv2.setMouseCallback(window_name, on_mouse)

    while not finished:
        preview = frame.copy()
        if center and radius > 0:
            cv2.circle(preview, center, radius, (0, 255, 255), 2)
        cv2.putText(
            preview,
            instructions,
            (10, preview.shape[0] - 10),
            cv2.FONT_HERSHEY_SIMPLEX,
            0.8,
            (255, 255, 255),
            2,
        )
        cv2.imshow(window_name, preview)
        key = cv2.waitKey(30) & 0xFF
        if key in (ord("c"), 27):
            finished = True
            radius = 0
    cv2.destroyWindow(window_name)
    return radius * 2 if radius else None


def main() -> None:
    parser = argparse.ArgumentParser(description="Track and crop image sequences per experiment.")
    parser.add_argument(
        "--config",
        default="../configs/opencv_tracker_v3.yaml",
        help="Path to YAML configuration file.",
    )
    args = parser.parse_args()

    config_path = Path(args.config).resolve()
    config_dir = config_path.parent
    config = load_config(config_path)

    input_cfg = config["input"]
    if "base_dir" in input_cfg:
        base_dir = resolve_path(input_cfg["base_dir"], config_dir)
    elif "image_glob" in input_cfg:
        base_dir = resolve_path(input_cfg["image_glob"], config_dir).parent
    else:
        raise ValueError("input.base_dir is required in the config.")

    if not base_dir.exists():
        raise FileNotFoundError(f"Input base directory does not exist: {base_dir}")

    raw_root_cfg = input_cfg.get(
        "raw_root", "/Volumes/Extreme SSD/deep-sea-particles-gm/raw"
    )
    raw_root = resolve_path(raw_root_cfg, config_dir)

    allowed_exts_raw = input_cfg.get("extensions", DEFAULT_IMAGE_EXTENSIONS)
    allowed_exts: List[str] = (
        [allowed_exts_raw]
        if isinstance(allowed_exts_raw, str)
        else list(allowed_exts_raw)
    )
    normalized_exts = {
        ext.lower() if ext.startswith(".") else f".{ext.lower()}" for ext in allowed_exts
    }

    subfolders = [
        entry
        for entry in sorted(base_dir.iterdir())
        if entry.is_dir() and collect_image_paths(entry, normalized_exts)
    ]
    if not subfolders:
        fallback_images = collect_image_paths(base_dir, normalized_exts)
        if fallback_images:
            subfolders = [base_dir]
        else:
            raise FileNotFoundError(f"No images found in {base_dir}")

    slider_cfg = input_cfg.get("slider", {})
    slider_enabled = bool(slider_cfg.get("enabled", True))
    threshold_cfg = input_cfg.get("threshold_slider", {})
    threshold_enabled = bool(threshold_cfg.get("enabled", True))

    processed_root_cfg = config["output"].get(
        "processed_root", "/Volumes/Extreme SSD/deep-sea-particles-gm/processed"
    )
    processed_root = resolve_path(processed_root_cfg, config_dir)
    try:
        relative_path = base_dir.relative_to(raw_root)
    except ValueError:
        relative_path = Path(base_dir.name)
    default_output_root = processed_root / relative_path

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

    configured_root = config["output"].get("root_dir")
    if configured_root:
        output_root = resolve_path(configured_root, config_dir)
    else:
        output_root = default_output_root

    output_root.mkdir(parents=True, exist_ok=True)
    metadata_src = ensure_experiment_metadata(base_dir)
    processed_metadata_dest = output_root / metadata_src.name
    metadata_text = metadata_src.read_text(encoding="utf-8")
    processed_metadata_dest.write_text(metadata_text, encoding="utf-8")
    try:
        exp_metadata = json.loads(metadata_text)
    except json.JSONDecodeError:
        exp_metadata = {}
        reader = csv.DictReader(io.StringIO(metadata_text))
        for row in reader:
            exp_metadata = row
            break

    crops_subdir_value = config["output"].get("crops_subdir")
    crops_subdir_value = (
        crops_subdir_value.strip() if isinstance(crops_subdir_value, str) else ""
    )
    if not crops_subdir_value:
        crops_subdir_value = DEFAULT_CROPS_SUBDIR
    crops_root = output_root / crops_subdir_value
    crops_root.mkdir(parents=True, exist_ok=True)

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
    folder_summaries: List[Dict[str, Any]] = []

    exit_requested = False
    for folder in subfolders:
        image_paths = collect_image_paths(folder, normalized_exts)
        if not image_paths:
            continue

        slider_params = {
            "window_name": format_slider_text(
                slider_cfg.get("window_name", "Select Start Frame"),
                folder.name,
                "Select Start Frame",
            ),
            "trackbar_name": slider_cfg.get("trackbar_name", "Frame"),
            "instructions": format_slider_text(
                slider_cfg.get(
                    "instructions",
                    "Move the slider and press 's' to select, 'x' to skip or 'q'/Esc to stop.",
                ),
                folder.name,
                "Move the slider and press 's' to select, 'x' to skip or 'q'/Esc to stop.",
            ),
            "start_index": int(slider_cfg.get("start_index", 0)),
            "font_scale": float(slider_cfg.get("font_scale", DEFAULT_FONT_SCALE)),
        }

        if slider_enabled:
            try:
                start_index, first_frame, skipped = choose_start_frame(
                    image_paths, slider_params
                )
            except KeyboardInterrupt:
                exit_requested = True
                break
            if skipped:
                continue
        else:
            start_index = max(
                0,
                min(slider_params["start_index"], len(image_paths) - 1),
            )
            first_frame = cv2.imread(str(image_paths[start_index]), cv2.IMREAD_COLOR)
        diameter_px = None
        diameter_um = None
        if first_frame is None:
            raise RuntimeError(f"Unable to read first image: {image_paths[start_index]}")

        tracker = create_tracker(tracker_name)
        window_track = config["tracker"].get("selection_window", "Select Tracking ROI")
        tracking_bbox = select_roi(window_track, first_frame, from_center=False)
        tracker.init(first_frame, tracking_bbox)

        threshold_value = None
        if threshold_enabled:
            x, y, w, h = [int(v) for v in tracking_bbox]
            if w > 0 and h > 0:
                roi = first_frame[y : y + h, x : x + w]
                threshold_settings = {
                    "window_name": format_slider_text(
                        threshold_cfg.get("window_name", "Threshold Tuning"),
                        folder.name,
                        "Threshold Tuning",
                    ),
                    "instructions": format_slider_text(
                threshold_cfg.get(
                    "instructions",
                    "Adjust threshold; press ENTER/SPACE to accept or 'c'/ESC to skip.",
                ),
                folder.name,
                "Adjust threshold; press ENTER/SPACE to accept or 'c'/ESC to skip.",
            ),
            "default_threshold": threshold_cfg.get("default_threshold"),
            "min": threshold_cfg.get("min", 0),
            "max": threshold_cfg.get("max", 255),
            "font_scale": float(threshold_cfg.get("font_scale", DEFAULT_FONT_SCALE)),
        }
                measured_px, measured_threshold = measure_diameter_from_threshold(
                    roi, threshold_settings
                )
                threshold_value = measured_threshold
                if measured_px:
                    diameter_px = measured_px
                    pixel_per_um = exp_metadata.get("pixelperum")
                    if pixel_per_um:
                        diameter_um = diameter_px / float(pixel_per_um)

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

        folder_crops_dir = crops_root / f"{folder.name}_cropped"
        folder_crops_dir.mkdir(parents=True, exist_ok=True)

        folder_rows: List[Dict[str, Any]] = []
        folder_file_times: List[float] = []
        processed_paths = image_paths[start_index:]
        tracking_failed = False
        failure_frame: Optional[int] = None

        for relative_index, image_path in enumerate(processed_paths):
            if relative_index == 0:
                frame = first_frame.copy()
            else:
                frame = cv2.imread(str(image_path), cv2.IMREAD_COLOR)
            if frame is None:
                raise RuntimeError(f"Unable to read image: {image_path}")
            img_h, img_w = frame.shape[:2]

            if relative_index == 0:
                success = True
                bbox = tracking_bbox
            else:
                success, updated = tracker.update(frame)
                if not success:
                    failure_frame = start_index + relative_index
                    tracking_failed = True
                    exit_requested = True
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
            cropped = frame[crop_y : crop_y + crop_h_eff, crop_x : crop_x + crop_w_eff]
            if not save_color:
                cropped = cv2.cvtColor(cropped, cv2.COLOR_BGR2GRAY)
            cropped = apply_rotation(cropped, rotation_mode)

            absolute_index = start_index + relative_index
            if keep_source_stem:
                out_name = f"{image_path.stem}_crop{save_ext}"
            else:
                out_name = f"frame_{absolute_index:06d}{save_ext}"
            out_path = folder_crops_dir / out_name
            cv2.imwrite(str(out_path), cropped)
            file_stat = image_path.stat()
            folder_file_times.append(file_stat.st_mtime)
            saved_h, saved_w = cropped.shape[:2]

            folder_rows.append(
                {
                    "subfolder": folder.name,
                    "frame_number": absolute_index,
                    "frame_index": relative_index,
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
                    "file_size_bytes": file_stat.st_size,
                    "saved_width": saved_w,
                    "saved_height": saved_h,
                    "rotation_applied": rotation_mode,
                    "hydrodynamic_diameter_px": diameter_px,
                    "hydrodynamic_diameter_um": diameter_um,
                    "threshold_value": threshold_value,
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
                    f"{tracker_name.upper()} | folder={folder.name} | frame={absolute_index} | ok={int(success)}",
                    (10, 25),
                    cv2.FONT_HERSHEY_SIMPLEX,
                    0.7,
                    (0, 255, 255),
                    2,
                )
                cv2.imshow(preview_window, view)
                key = cv2.waitKey(1) & 0xFF
                if key in (27, ord("q"), ord("t"), ord("T")):
                    exit_requested = True
                    break

        cv2.destroyAllWindows()
        if tracking_failed:
            print(
                f"Tracker lost the object in folder '{folder.name}' at frame {failure_frame}; "
                "stopping further folders."
            )

        folder_frame_rate: Optional[float] = None
        sequence_duration: Optional[float] = None
        if len(folder_file_times) >= 2:
            sequence_duration = folder_file_times[-1] - folder_file_times[0]
            if sequence_duration > 0.0 and len(folder_rows) > 1:
                folder_frame_rate = (len(folder_rows) - 1) / sequence_duration

        delta_t = 1.0 / folder_frame_rate if folder_frame_rate else 0.0
        for idx, row in enumerate(folder_rows):
            row["time_seconds"] = idx * delta_t if folder_frame_rate else 0.0

        rows.extend(folder_rows)
        folder_summaries.append(
            {
                "subfolder": folder.name,
                "num_frames": len(folder_rows),
                "frame_rate_hz": folder_frame_rate,
                "sequence_duration_seconds": sequence_duration,
                "crops_dir": str(folder_crops_dir),
                "hydrodynamic_diameter_um": diameter_um,
                "threshold_value": threshold_value,
            }
        )
        if exit_requested:
            break
    if rows:
        fieldnames = list(rows[0].keys())
        with metadata_csv.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            for row in rows:
                writer.writerow(row)

    summary = {
        "input_directory": str(base_dir),
        "allowed_extensions": sorted(normalized_exts),
        "total_saved_frames": len(rows),
        "tracker_type": tracker_name,
        "rotation_mode": rotation_mode,
        "num_subfolders": len(folder_summaries),
        "folder_summaries": folder_summaries,
        "metadata_csv": str(metadata_csv),
        "crops_root": str(crops_root),
        "output_root": str(output_root),
    }
    with summary_json.open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)

    print(json.dumps(summary, indent=2))

if __name__ == "__main__":
    main()
