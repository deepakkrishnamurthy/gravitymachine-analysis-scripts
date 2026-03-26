#!/usr/bin/env python3
"""Estimate uniform background velocity from processed PTV data."""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path

import cv2
import numpy as np
import pandas as pd
import trackpy as tp
import yaml


DEFAULT_MAX_FRAMES = 50


def resolve_path(path_like, base_dir: Path) -> Path:
    candidate = Path(path_like)
    return candidate if candidate.is_absolute() else (base_dir / candidate).resolve()


def natural_sort_key(path: Path):
    return [int(part) if part.isdigit() else part.lower() for part in re.split(r'(\d+)', path.name)]


def collect_images(folder: Path, allowed_exts=None):
    allowed_exts = [ext.lower() for ext in (allowed_exts or [])]
    candidates = [p for p in folder.iterdir() if p.is_file() and not p.name.startswith('.')]
    if allowed_exts:
        candidates = [p for p in candidates if p.suffix.lower() in allowed_exts]
    return sorted(candidates, key=natural_sort_key)


def load_experiment_metadata(experiment_dir: Path) -> dict[str, str]:
    for suffix in ('metadata.json', 'metadata.csv'):
        path = experiment_dir / suffix
        if not path.exists():
            continue
        if path.suffix.lower() == '.json':
            return yaml.safe_load(path.read_text(encoding='utf-8')) or {}
        with path.open('r', encoding='utf-8') as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                return row
    return {}


def build_far_field_mask(image: np.ndarray, threshold: float, iterations: int = 5) -> tuple[np.ndarray, np.ndarray]:
    _, mask = cv2.threshold(image, int(threshold), 255, cv2.THRESH_BINARY)
    if not np.any(mask):
        _, mask = cv2.threshold(image, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (11, 11))
    expanded = cv2.dilate(mask, kernel, iterations=iterations)
    far_field = cv2.bitwise_not(expanded)
    return far_field, expanded


def histogram_stretch(image: np.ndarray, lower_pct: float = 0.1, upper_pct: float = 99.9) -> np.ndarray:
    arr = np.asarray(image, dtype=np.float32)
    if arr.ndim != 2:
        arr = cv2.cvtColor(arr, cv2.COLOR_BGR2GRAY)
    low = np.percentile(arr, lower_pct)
    high = np.percentile(arr, upper_pct)
    if high <= low:
        return np.clip(arr, 0, 255).astype(np.uint8)
    stretched = (arr - low) * (255.0 / (high - low))
    return np.clip(stretched, 0, 255).astype(np.uint8)


def filter_tracks_by_mask(tracks: pd.DataFrame, mask: np.ndarray, ycol='y', xcol='x', particle_col='particle', drop_entire_track=True):
    h, w = mask.shape
    yy = np.rint(tracks[ycol].to_numpy()).astype(int)
    xx = np.rint(tracks[xcol].to_numpy()).astype(int)
    yy = np.clip(yy, 0, h - 1)
    xx = np.clip(xx, 0, w - 1)
    inside = mask[yy, xx].astype(bool)
    out = tracks.copy()
    out['inside_mask'] = inside
    if drop_entire_track:
        bad_particles = out.groupby(particle_col)['inside_mask'].any()
        bad_particles = bad_particles[bad_particles].index
        out = out.loc[~out[particle_col].isin(bad_particles)].copy()
    else:
        out = out.loc[~out['inside_mask']].copy()
    return out.drop(columns='inside_mask')


def compute_particle_velocities(tracks: pd.DataFrame) -> pd.DataFrame:
    tracks = tracks.reset_index(drop=True).copy()
    rows = []
    for pid, grp in tracks.groupby('particle'):
        grp = grp.sort_values('frame')
        start = grp.iloc[0]
        end = grp.iloc[-1]
        start_frame = int(start['frame'])
        end_frame = int(end['frame'])
        n_steps = end_frame - start_frame
        if n_steps <= 0:
            continue
        dx = end['x'] - start['x']
        dy = end['y'] - start['y']
        vx = dx / n_steps
        vy = dy / n_steps
        rows.append({
            'particle': pid,
            'start_frame': start_frame,
            'end_frame': end_frame,
            'n_steps': n_steps,
            'x_start': start['x'],
            'y_start': start['y'],
            'x_end': end['x'],
            'y_end': end['y'],
            'dx': dx,
            'dy': dy,
            'vx': vx,
            'vy': vy,
            'speed': np.sqrt(vx**2 + vy**2),
        })
    return pd.DataFrame(rows)


def summarize_particle_velocities(particle_vel: pd.DataFrame) -> dict[str, float]:
    if particle_vel.empty:
        return {
            'n_particles': 0,
            'mean_vx': 0.0,
            'std_vx': 0.0,
            'mean_vy': 0.0,
            'std_vy': 0.0,
            'mean_speed': 0.0,
            'std_speed': 0.0,
        }
    return {
        'n_particles': len(particle_vel),
        'mean_vx': particle_vel['vx'].mean(),
        'std_vx': particle_vel['vx'].std(ddof=1),
        'mean_vy': particle_vel['vy'].mean(),
        'std_vy': particle_vel['vy'].std(ddof=1),
        'mean_speed': particle_vel['speed'].mean(),
        'std_speed': particle_vel['speed'].std(ddof=1),
    }


def draw_tracks_overlay(base_image: np.ndarray, tracks: pd.DataFrame, max_tracks=200) -> np.ndarray:
    overlay = cv2.cvtColor(base_image, cv2.COLOR_GRAY2BGR)
    particles = tracks['particle'].unique()[:max_tracks]
    colors = [tuple(int(v) for v in np.random.randint(0, 255, size=3)) for _ in particles]
    for color, pid in zip(colors, particles):
        grp = tracks[tracks['particle'] == pid]
        pts = np.rint(grp[['x', 'y']].to_numpy()).astype(int)
        for i in range(1, len(pts)):
            pt0 = tuple(np.clip(pts[i - 1], 0, np.array(overlay.shape[1::-1]) - 1).tolist())
            pt1 = tuple(np.clip(pts[i], 0, np.array(overlay.shape[1::-1]) - 1).tolist())
            cv2.line(overlay, pt0, pt1, color, 1)
    return overlay


def process_folder(
    cropped_folder: Path,
    folder_name: str,
    summary_row: dict[str, object],
    allowed_exts: list[str],
    metadata_record: dict[str, str],
    output_root: Path,
    max_frames: int,
    feature_size: int = 11,
    feature_minmass: float = 400,
    link_search_distance: int = 30,
    link_memory: int = 3,
    min_track_length: int = 10,
):
    image_paths = collect_images(cropped_folder, allowed_exts)
    if not image_paths:
        raise RuntimeError(f'No images found in {cropped_folder}')
    image_paths = image_paths[:max_frames]
    grayscale_images = []
    for path in image_paths:
        image = cv2.imread(str(path), cv2.IMREAD_GRAYSCALE)
        if image is None:
            continue
        grayscale_images.append(image)
    if len(grayscale_images) < 2:
        raise RuntimeError('Need at least two valid frames for PTV')
    threshold_value = summary_row.get('threshold_value')
    frame_rate = summary_row.get('frame_rate_hz')
    if threshold_value is None or frame_rate is None:
        raise RuntimeError('Missing threshold or frame rate in folder summary')
    _, expanded_mask = build_far_field_mask(grayscale_images[0], threshold_value, iterations=20)
    expanded_mask = expanded_mask.astype('uint8', copy=False)
    stretch_sequence = [histogram_stretch(img) for img in grayscale_images]
    features = tp.batch(stretch_sequence, diameter=feature_size, minmass=feature_minmass)
    tracks = tp.link(features, search_range=link_search_distance, memory=link_memory)
    tracks = tp.filter_stubs(tracks, min_track_length)
    tracks = filter_tracks_by_mask(tracks, expanded_mask, drop_entire_track=True)
    particle_vel = compute_particle_velocities(tracks)
    summary = summarize_particle_velocities(particle_vel)
    pixel_per_um = float(metadata_record.get('pixelperum', 1.36)) or 1.36
    scale = 1.0 / pixel_per_um if pixel_per_um else 1.0
    mean_speed = summary['mean_speed'] * frame_rate * scale
    std_speed = summary['std_speed'] * frame_rate * scale
    mean_vx = summary['mean_vx'] * frame_rate * scale
    mean_vy = summary['mean_vy'] * frame_rate * scale
    std_vx = summary['std_vx'] * frame_rate * scale
    std_vy = summary['std_vy'] * frame_rate * scale
    overlay_root = output_root / 'diagnostics'
    overlay_root.mkdir(parents=True, exist_ok=True)
    overlay_path = overlay_root / f'{folder_name}_tracks.png'
    overlay_img = draw_tracks_overlay(stretch_sequence[0], tracks)
    cv2.imwrite(str(overlay_path), overlay_img)
    diameter_um = summary_row.get('hydrodynamic_diameter_um')
    if diameter_um is None:
        diameter_um = metadata_record.get('hydrodynamic_diameter_um') or metadata_record.get('object_size_um')
    stats = {
        'subfolder': folder_name,
        'threshold_value': threshold_value,
        'frame_rate_hz': frame_rate,
        'pixelperum': pixel_per_um,
        'mean_speed_um_per_s': mean_speed,
        'std_speed_um_per_s': std_speed,
        'n_particles': summary['n_particles'],
        'mean_vx_um_per_s': mean_vx,
        'std_vx_um_per_s': std_vx,
        'mean_vy_um_per_s': mean_vy,
        'std_vy_um_per_s': std_vy,
        'overlay_path': str(overlay_path),
        'particle_type': metadata_record.get('particle_type', ''),
        'particle_diameter_um': diameter_um,
    }
    return stats


def gather_configs(config_dir: Path, single_config: Path | None) -> list[Path]:
    if single_config:
        return [single_config]
    if not config_dir.exists():
        raise FileNotFoundError(f'Configs directory does not exist: {config_dir}')
    return sorted(config_dir.glob('*.yaml'))


def write_summary(output_path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    df = pd.DataFrame(rows)
    output_path.write_text('', encoding='utf-8')
    df.to_csv(output_path, index=False)


def process_config(config_path: Path, args: argparse.Namespace) -> None:
    config_dir = config_path.parent
    config = yaml.safe_load(config_path.read_text(encoding='utf-8')) or {}
    input_cfg = config.get('input', {})
    base_dir = resolve_path(input_cfg['base_dir'], config_dir)
    raw_root = resolve_path(input_cfg.get('raw_root', '/Volumes/Extreme SSD/deep-sea-particles-gm/raw'), config_dir)
    try:
        relative = base_dir.relative_to(raw_root)
    except ValueError:
        relative = Path(base_dir.name)
    processed_root = resolve_path(input_cfg.get('processed_root', '/Volumes/Extreme SSD/deep-sea-particles-gm/processed'), config_dir)
    processed_experiment_dir = processed_root / relative
    crops_subdir = config.get('output', {}).get('crops_subdir', 'cropped') or 'cropped'
    crops_root = processed_experiment_dir / crops_subdir
    summary_path = processed_experiment_dir / config.get('output', {}).get('summary_json', 'frame_rate_summary.json')
    if not summary_path.exists():
        print(f'Summary missing for {config_path.name}, skipping')
        return
    summary_data = json.loads(summary_path.read_text(encoding='utf-8'))
    folder_summaries = summary_data.get('folder_summaries', [])
    metadata_record = load_experiment_metadata(processed_experiment_dir)
    allowed_exts = [ext.lower() for ext in config.get('input', {}).get('extensions', ['.tif', '.png', '.bmp'])]
    output_dir = processed_experiment_dir / 'ptv_uniform_background'
    output_dir.mkdir(parents=True, exist_ok=True)
    summary_rows = []
    for folder_summary in folder_summaries:
        folder_name = folder_summary.get('subfolder')
        if not folder_name:
            continue
        cropped_folder = crops_root / f'{folder_name}_cropped'
        if not cropped_folder.exists():
            continue
        row = process_folder(
            cropped_folder,
            folder_name,
            folder_summary,
            allowed_exts,
            metadata_record,
            output_dir,
            max_frames=args.max_frames,
        )
        row['experiment'] = processed_experiment_dir.name
        summary_rows.append(row)
    write_summary(output_dir / 'ptv_uniform_background_summary.csv', summary_rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Run uniform background PTV diagnostics for all configs.')
    repo_root = Path(__file__).resolve().parent.parent
    parser.add_argument('--configs-dir', type=Path, default=repo_root / 'configs', help='Directory containing tracker configs.')
    parser.add_argument('--config', type=Path, help='Run only a single config file.')
    parser.add_argument('--max-frames', type=int, default=DEFAULT_MAX_FRAMES, help='Maximum number of frames to use per folder.')
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    configs = gather_configs(args.configs_dir, args.config)
    if not configs:
        raise RuntimeError('No config files found to process.')
    for cfg in configs:
        print(f'Processing config: {cfg}')
        try:
            process_config(cfg, args)
        except Exception as exc:
            print(f'Failed to process {cfg.name}: {exc}')
            continue


if __name__ == '__main__':
    main()
