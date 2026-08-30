from __future__ import annotations

from dataclasses import replace
import shutil
import subprocess

import numpy as np
import pytest

from sim_swim.render.project2d import _camera_center_2d, project_states
from sim_swim.render.render3d import (
    _frame_status_lines,
    _hook_edges,
    _resolve_view_range_um,
    _run_tumble_label,
    _select_frames,
    save_swim_movie,
)
from sim_swim.render.video_writer import open_mp4_writer
from sim_swim.sim.flagella_geometry import FlagellaRig
from sim_swim.sim.core import SimulationState, Simulator
from sim_swim.sim.params import SimulationConfig


def _make_cfg(
    *,
    center_body_in_2d: bool,
    follow_camera_2d: bool,
    enable_switching: bool,
) -> SimulationConfig:
    return SimulationConfig.from_dict(
        {
            "scale": {"b_um": 1.0, "bead_radius_a_over_b": 0.1},
            "body": {
                "prism": {
                    "n_prism": 3,
                    "dz_over_b": 0.5,
                    "radius_over_b": 0.5,
                    "axis": "x",
                },
                "length_total_um": 2.0,
            },
            "flagella": {
                "n_flagella": 1,
                "placement_mode": "uniform",
                "discretization": {"ds_over_b": 0.58},
                "bond_L_over_b": 0.58,
                "length_over_b": 2.32,
                "helix_init": {"radius_over_b": 0.25, "pitch_over_b": 2.5},
            },
            "fluid": {"viscosity_Pa_s": 1.0e-3},
            "motor": {
                "torque_Nm": 1.0e-18,
                "reverse_n_flagella": 1,
                "enable_switching": enable_switching,
            },
            "time": {"duration_s": 0.02, "dt_s": 1.0e-3},
            "render": {
                "center_body_in_2d": center_body_in_2d,
                "follow_camera_2d": follow_camera_2d,
            },
        }
    )


def _state() -> SimulationState:
    return SimulationState(
        t=0.0,
        position_um=(3.0, -2.0, 1.0),
        quaternion=(0.0, 0.0, 0.0, 1.0),
        velocity_um_s=(0.0, 0.0, 0.0),
        omega_rad_s=(0.0, 0.0, 0.0),
        bead_positions_um=np.zeros((1, 3), dtype=float),
        flag_states=(2,),
        reverse_flagella=(0,),
    )


def test_camera_center_2d_uses_body_center_when_enabled() -> None:
    cfg = _make_cfg(
        center_body_in_2d=True,
        follow_camera_2d=False,
        enable_switching=False,
    )
    cam = _camera_center_2d(_state(), cfg)

    assert np.allclose(cam, np.array([3.0, -2.0]))


def test_camera_center_2d_keeps_legacy_behavior_when_disabled() -> None:
    cfg_fixed = _make_cfg(
        center_body_in_2d=False,
        follow_camera_2d=False,
        enable_switching=False,
    )
    cam_fixed = _camera_center_2d(_state(), cfg_fixed)
    assert np.allclose(cam_fixed, np.zeros(2, dtype=float))

    cfg_follow = _make_cfg(
        center_body_in_2d=False,
        follow_camera_2d=True,
        enable_switching=False,
    )
    cam_follow = _camera_center_2d(_state(), cfg_follow)
    assert np.allclose(cam_follow, np.array([3.0, -2.0]))


def test_run_tumble_label_is_always_run_when_switching_disabled() -> None:
    cfg = _make_cfg(
        center_body_in_2d=True,
        follow_camera_2d=False,
        enable_switching=False,
    )

    assert _run_tumble_label(_state(), cfg) == "RUN"


def test_frame_status_lines_include_time_torque_and_camera_mode() -> None:
    cfg = _make_cfg(
        center_body_in_2d=True,
        follow_camera_2d=False,
        enable_switching=False,
    )

    lines = _frame_status_lines(_state(), cfg)

    assert lines[0] == "RUN"
    assert "t = 0.000 τ (0.000000 s, 0 steps)" in lines
    assert "motor_torque_Nm = 1.000e-18" in lines
    assert "follow_camera_3d = True" in lines


def test_frame_status_lines_always_include_tau_seconds_and_steps() -> None:
    base_cfg = _make_cfg(
        center_body_in_2d=True,
        follow_camera_2d=False,
        enable_switching=False,
    )
    cfg = replace(
        base_cfg,
        render=replace(
            base_cfg.render,
            timestamp_3d=False,
            timestamp_fmt="t = {t:.3f} s (τ = {tau_s:.3e} s)",
        ),
    )

    lines = _frame_status_lines(_state(), cfg)

    assert lines[1] == "t = 0.000 τ (0.000000 s, 0 steps)"


def test_hook_edges_expand_triplets_into_two_segments() -> None:
    triplets = np.array([[1, 4, 5], [2, 6, 7]], dtype=int)

    edges = _hook_edges(triplets)

    assert edges.shape == (4, 2)
    assert np.array_equal(edges, np.array([[1, 4], [2, 6], [4, 5], [6, 7]]))


def test_view_range_defaults_to_3_when_no_flagella_are_present() -> None:
    cfg = _make_cfg(
        center_body_in_2d=True,
        follow_camera_2d=False,
        enable_switching=False,
    )
    rig = FlagellaRig(
        body_layer_indices=[np.array([0, 1, 2, 3], dtype=int)],
        body_ring_edges=np.array([[0, 1]], dtype=int),
        body_vertical_edges=np.array([[1, 2]], dtype=int),
        body_spring_edges=np.array([[0, 1]], dtype=int),
        flagella_indices=[],
        hook_triplets=np.array([], dtype=int).reshape(0, 3),
    )

    assert _resolve_view_range_um(cfg, rig) == 3.0


def test_select_frames_uses_3d_fps_when_not_saving_all_steps() -> None:
    states = [replace(_state(), t=0.1 * idx) for idx in range(11)]

    assert _select_frames(states, out_all_steps_3d=True, fps_hint=2.0) == states

    selected = _select_frames(states, out_all_steps_3d=False, fps_hint=2.0)

    assert [state.t for state in selected] == [0.0, 0.5, 1.0]


def test_mp4_writer_requires_shared_h264_encoder(tmp_path, monkeypatch) -> None:
    class DummyWriter:
        def isOpened(self) -> bool:
            return True

        def release(self) -> None:
            return None

    monkeypatch.setattr(
        "sim_swim.render.video_writer._FFmpegVideoWriter",
        lambda *args, **kwargs: DummyWriter(),
    )

    selection = open_mp4_writer(tmp_path / "movie.mp4", fps=25.0, frame_size=(10, 10))

    assert selection.selected_codec == "libx264"
    assert selection.attempted_codecs == ("ffmpeg:libx264",)


def test_mp4_writer_rejects_non_h264_candidates(tmp_path) -> None:
    with pytest.raises(ValueError, match="Only H.264 output"):
        open_mp4_writer(
            tmp_path / "movie.mp4",
            fps=25.0,
            frame_size=(10, 10),
            codec_candidates=("mp4v",),
        )


@pytest.mark.skipif(shutil.which("ffmpeg") is None, reason="FFmpeg is unavailable")
def test_mp4_writer_emits_h264_yuv420p(tmp_path) -> None:
    path = tmp_path / "movie.mp4"
    selection = open_mp4_writer(path, fps=10.0, frame_size=(4, 2))
    selection.writer.write(np.zeros((2, 4, 3), dtype=np.uint8))
    selection.writer.release()

    completed = subprocess.run(
        [
            "ffprobe",
            "-v",
            "error",
            "-select_streams",
            "v:0",
            "-show_entries",
            "stream=codec_name,pix_fmt,profile",
            "-of",
            "default=nokey=1:noprint_wrappers=1",
            str(path),
        ],
        check=True,
        capture_output=True,
        text=True,
    )

    assert completed.stdout.splitlines() == ["h264", "High", "yuv420p"]


@pytest.mark.skipif(shutil.which("ffmpeg") is None, reason="FFmpeg is unavailable")
def test_mp4_writer_pads_odd_dimensions_for_h264(tmp_path) -> None:
    path = tmp_path / "odd-size.mp4"
    selection = open_mp4_writer(path, fps=10.0, frame_size=(3, 3))
    selection.writer.write(np.zeros((3, 3, 3), dtype=np.uint8))
    selection.writer.release()

    completed = subprocess.run(
        [
            "ffprobe",
            "-v",
            "error",
            "-select_streams",
            "v:0",
            "-show_entries",
            "stream=codec_name,pix_fmt,width,height",
            "-of",
            "csv=p=0",
            str(path),
        ],
        check=True,
        capture_output=True,
        text=True,
    )

    assert completed.stdout.strip() == "h264,4,4,yuv420p"


def test_save_swim_movie_emits_render_outputs(tmp_path, monkeypatch) -> None:
    cfg = _make_cfg(
        center_body_in_2d=True,
        follow_camera_2d=False,
        enable_switching=False,
    )
    cfg = replace(cfg, render=replace(cfg.render, save_frames_3d=True))

    beads = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
            [3.0, 1.0, 0.0],
            [3.0, 2.0, 0.0],
        ],
        dtype=float,
    )
    state = SimulationState(
        t=0.0,
        position_um=(0.0, 0.0, 0.0),
        quaternion=(0.0, 0.0, 0.0, 1.0),
        velocity_um_s=(0.0, 0.0, 0.0),
        omega_rad_s=(0.0, 0.0, 0.0),
        bead_positions_um=beads,
        flag_states=(0,),
        reverse_flagella=(0,),
    )
    rig = FlagellaRig(
        body_layer_indices=[np.array([0, 1, 2, 3], dtype=int)],
        body_ring_edges=np.array([[0, 1], [1, 2]], dtype=int),
        body_vertical_edges=np.array([[2, 3]], dtype=int),
        body_spring_edges=np.array([[0, 1], [1, 2], [2, 3]], dtype=int),
        flagella_indices=[np.array([3, 4, 5], dtype=int)],
        hook_triplets=np.array([[1, 4, 5]], dtype=int),
    )

    class DummyWriter:
        def __init__(self) -> None:
            self.frame_shape = None
            self.write_calls = 0
            self.released = False

        def isOpened(self) -> bool:
            return True

        def write(self, frame) -> None:
            self.frame_shape = frame.shape
            self.write_calls += 1

        def release(self) -> None:
            self.released = True

    writer_calls: list[DummyWriter] = []

    def make_writer(*args, **kwargs) -> DummyWriter:
        writer = DummyWriter()
        writer_calls.append(writer)
        return writer

    monkeypatch.setattr(
        "sim_swim.render.video_writer._FFmpegVideoWriter",
        make_writer,
    )

    result = save_swim_movie([state], cfg, rig, tmp_path)

    assert (tmp_path / "swim3d_final.png").exists()
    assert (tmp_path / "frames_3d" / "frame_000000.png").exists()
    assert len(writer_calls) == 1
    assert writer_calls[0].write_calls > 0
    assert result is not None
    assert result.selected_codec == "libx264"
    assert result.frame_count == 1


def test_save_swim_movie_can_overlay_flagella_helix_axes(
    tmp_path,
    monkeypatch,
) -> None:
    cfg = _make_cfg(
        center_body_in_2d=True,
        follow_camera_2d=False,
        enable_switching=False,
    )
    cfg = replace(
        cfg,
        render=replace(cfg.render, show_flagella_helix_axis_3d=True),
    )

    beads = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
            [3.0, 1.0, 0.0],
            [3.4, 2.0, 0.2],
            [2.9, 3.0, -0.2],
        ],
        dtype=float,
    )
    state = SimulationState(
        t=0.0,
        position_um=(0.0, 0.0, 0.0),
        quaternion=(0.0, 0.0, 0.0, 1.0),
        velocity_um_s=(0.0, 0.0, 0.0),
        omega_rad_s=(0.0, 0.0, 0.0),
        bead_positions_um=beads,
        flag_states=(0,),
        reverse_flagella=(0,),
    )
    rig = FlagellaRig(
        body_layer_indices=[np.array([0, 1, 2, 3], dtype=int)],
        body_ring_edges=np.array([[0, 1], [1, 2]], dtype=int),
        body_vertical_edges=np.array([[2, 3]], dtype=int),
        body_spring_edges=np.array([[0, 1], [1, 2], [2, 3]], dtype=int),
        flagella_indices=[np.array([3, 4, 5, 6], dtype=int)],
        hook_triplets=np.array([[1, 4, 5]], dtype=int),
    )

    class DummyWriter:
        def isOpened(self) -> bool:
            return True

        def write(self, frame) -> None:
            self.frame_shape = frame.shape

        def release(self) -> None:
            self.released = True

    monkeypatch.setattr(
        "sim_swim.render.video_writer._FFmpegVideoWriter",
        lambda *args, **kwargs: DummyWriter(),
    )

    save_swim_movie([state], cfg, rig, tmp_path)

    assert (tmp_path / "swim3d_final.png").exists()


@pytest.mark.parametrize(
    "projection_mode_2d", ["bead_projection", "body_capsule_orthographic_v1"]
)
def test_project_states_reports_selected_video_codec(
    tmp_path, monkeypatch, projection_mode_2d: str
) -> None:
    base_cfg = _make_cfg(
        center_body_in_2d=True,
        follow_camera_2d=False,
        enable_switching=False,
    )
    cfg = replace(
        base_cfg,
        render=replace(base_cfg.render, projection_mode_2d=projection_mode_2d),
    )
    beads = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
            [3.0, 1.0, 0.0],
            [3.0, 2.0, 0.0],
        ],
        dtype=float,
    )
    state = SimulationState(
        t=0.0,
        position_um=(0.0, 0.0, 0.0),
        quaternion=(0.0, 0.0, 0.0, 1.0),
        velocity_um_s=(0.0, 0.0, 0.0),
        omega_rad_s=(0.0, 0.0, 0.0),
        bead_positions_um=beads,
        flag_states=(0,),
        reverse_flagella=(0,),
    )
    rig = FlagellaRig(
        body_layer_indices=[np.array([0, 1, 2, 3], dtype=int)],
        body_ring_edges=np.array([[0, 1], [1, 2]], dtype=int),
        body_vertical_edges=np.array([[2, 3]], dtype=int),
        body_spring_edges=np.array([[0, 1], [1, 2], [2, 3]], dtype=int),
        flagella_indices=[np.array([3, 4, 5], dtype=int)],
        hook_triplets=np.array([[1, 4, 5]], dtype=int),
    )

    class DummyWriter:
        def __init__(self) -> None:
            self.write_calls = 0

        def isOpened(self) -> bool:
            return True

        def write(self, frame) -> None:
            self.write_calls += 1

        def release(self) -> None:
            pass

    monkeypatch.setattr(
        "sim_swim.render.video_writer._FFmpegVideoWriter",
        lambda *args, **kwargs: DummyWriter(),
    )

    result = project_states([state], cfg, rig, tmp_path)

    assert result is not None
    assert result.selected_codec == "libx264"
    assert result.frame_count == 1
    assert result.frame_size == (256, 256)


def test_project_states_rejects_unknown_projection_mode(tmp_path) -> None:
    base_cfg = _make_cfg(
        center_body_in_2d=True,
        follow_camera_2d=False,
        enable_switching=False,
    )
    cfg = replace(
        base_cfg,
        render=replace(base_cfg.render, projection_mode_2d="unknown"),
    )
    simulator = Simulator(base_cfg)

    with pytest.raises(ValueError, match="render.projection_mode_2d"):
        project_states([simulator._observe(0.0, None)], cfg, simulator.rig, tmp_path)
