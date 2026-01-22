from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional

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
    overrides: List[str] = typer.Argument(
        None,
        help="Optional overrides as key=value (e.g., data.video_path=path.avi)",
    ),
) -> None:
    """CLI entrypoint for tracking + butt estimation with optional overrides."""
    override_dict: Dict[str, Any] = _parse_overrides(overrides)
    if video_path is not None:
        override_dict.setdefault("data", {})["video_path"] = str(video_path)

    run_tracking_butt(
        config_path=config,
        save_contour_flag=save_contour,
        overrides=override_dict or None,
    )


def _parse_overrides(items: List[str] | None) -> Dict[str, Any]:
    """Parse flat key=value overrides into nested dict (e.g., data.video_path=...)."""
    if not items:
        return {}
    result: Dict[str, Any] = {}
    for raw in items:
        if "=" not in raw:
            raise typer.BadParameter(f"Override must be key=value, got '{raw}'")
        key, value = raw.split("=", 1)
        parts = key.split(".")
        node = result
        for part in parts[:-1]:
            node = node.setdefault(part, {})
        node[parts[-1]] = value
    return result


if __name__ == "__main__":
    typer.run(main)
