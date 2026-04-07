from __future__ import annotations

import numpy as np

from sim_swim.render.project2d import _camera_center_2d
from sim_swim.render.render3d import (
    _hook_edges,
    _resolve_view_range_um,
    _run_tumble_label,
    save_swim_movie,
)
from sim_swim.sim.flagella_geometry import FlagellaRig
from sim_swim.sim.core import SimulationState
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


def test_save_swim_movie_emits_render_outputs(tmp_path, monkeypatch) -> None:
    cfg = _make_cfg(
        center_body_in_2d=True,
        follow_camera_2d=False,
        enable_switching=False,
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
        def write(self, frame) -> None:
            self.frame_shape = frame.shape

        def release(self) -> None:
            pass

    monkeypatch.setattr(
        "sim_swim.render.render3d.cv2.VideoWriter",
        lambda *args, **kwargs: DummyWriter(),
    )

    save_swim_movie([state], cfg, rig, tmp_path)

    assert (tmp_path / "swim3d_final.png").exists()
    assert not (tmp_path / "swim3d.mp4").exists()
    assert (tmp_path / "frames_3d" / "frame_000000.png").exists()
