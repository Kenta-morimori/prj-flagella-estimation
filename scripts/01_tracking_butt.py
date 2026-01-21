from __future__ import annotations

from pathlib import Path

import typer

from flagella_estimation.tracking_butt.pipeline import run_tracking_butt


def main(
    config: Path = typer.Option(
        Path("conf/config.yaml"),
        "--config",
        "-c",
        help="Path to config yaml",
    ),
    save_contour: bool = typer.Option(
        False,
        "--save-contour",
        help="Enable contour saving (tracking_butt.save.contour=true)",
    ),
) -> None:
    run_tracking_butt(config_path=config, save_contour_flag=save_contour)


if __name__ == "__main__":
    typer.run(main)
