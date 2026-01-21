from __future__ import annotations

from pathlib import Path
from typing import Optional

import typer

from flagella_estimation.tracking_butt.pipeline import run_tracking_butt


def main(
    config: Path = typer.Option(
        Path("conf/config.yaml"),
        "--config",
        "-c",
        help="Path to config yaml",
    ),
    video_path: Optional[Path] = typer.Option(
        None,
        "--video-path",
        help="Override data.video_path in config (supports mp4/avi etc.)",
    ),
    save_contour: bool = typer.Option(
        False,
        "--save-contour",
        help="Enable contour saving (tracking_butt.save.contour=true)",
    ),
) -> None:
    """CLI entrypoint for tracking + butt estimation with optional overrides."""
    overrides = {}
    if video_path is not None:
        overrides["data"] = {"video_path": str(video_path)}

    run_tracking_butt(
        config_path=config, save_contour_flag=save_contour, overrides=overrides or None
    )


if __name__ == "__main__":
    typer.run(main)
