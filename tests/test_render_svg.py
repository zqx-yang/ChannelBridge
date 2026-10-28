"""The SVG assembler: silhouette, mouths, outlines, placement, framing.

Every function here decides something that is *silent* in the finished
picture -- a channel painted over solvent, a label leaning on the line it
annotates, an arc drawn the long way round -- so each gets a case built from
geometry whose right answer is known by construction.
"""

from __future__ import annotations

import math
import re
import xml.dom.minidom
from typing import ClassVar

import numpy as np
import pytest

from channelbridge.render.style import DEFAULT_STYLE, RenderStyle, rgb_to_hex
from channelbridge.render.svg import (
    AngleMarker,
    Mouth,
    SvgStyle,
    _arc_flags,
    _Layout,
    _signed_delta,
    angle_marker,
    content_box,
    find_mouths,
    open_runs,
    path_closed,
    path_open,
    place_block,
    svg_document,
    trace_contours,
    tunnel_silhouette,
)

H, W = 200, 300


def _disc(shape, cx, cy, r, value=0.0, base=1.0):
    img = np.full(shape, base, float)
    yy, xx = np.mgrid[0:shape[0], 0:shape[1]]
    img[(xx - cx) ** 2 + (yy - cy) ** 2 <= r * r] = value
    return img


def _bar(shape, x0, x1, y0, y1, value=0.0, base=1.0):
    img = np.full(shape, base, float)
    img[y0:y1, x0:x1] = value
    return img


class TestSilhouette:
    """The clip to the enzymes' footprint -- the rule that red never floats."""

    def test_the_part_of_the_mask_outside_the_footprint_is_dropped(self):
        mask = _bar((H, W), 20, 280, 90, 110)          # a bar right across
        foot = _bar((H, W), 20, 150, 60, 140)          # only the left half is protein
        sil = tunnel_silhouette(mask, foot, smooth=0.5, dilate=1)
        assert sil[:, 20:140].any(), "the part inside the enzyme must survive"
        assert not sil[:, 200:].any(), "nothing may be drawn over solvent"

    def test_a_mask_entirely_in_solvent_leaves_nothing(self):
        mask = _bar((H, W), 200, 280, 90, 110)
        foot = _bar((H, W), 20, 150, 60, 140)
        assert not tunnel_silhouette(mask, foot, smooth=0.5, dilate=1).any()

    def test_an_opaque_chain_in_front_hides_the_channel(self):
        mask = _bar((H, W), 20, 280, 90, 110)
        foot = _bar((H, W), 10, 290, 60, 140)
        occ = _bar((H, W), 100, 160, 60, 140, value=0.0)
        plain = tunnel_silhouette(mask, foot, smooth=0.5, dilate=1)
        clipped = tunnel_silhouette(mask, foot, smooth=0.5, dilate=1, occlusion=occ)
        assert plain[:, 110:150].any()
        assert not clipped[:, 110:150].any(), "the overlay must not cross an opaque chain"
        assert clipped[:, 30:80].any(), "the rest of the channel stays"

    def test_renders_of_different_sizes_are_refused(self):
        with pytest.raises(ValueError, match="share a camera"):
            tunnel_silhouette(np.ones((10, 10)), np.ones((12, 12)))

    def test_an_occlusion_render_of_the_wrong_size_is_refused(self):
        with pytest.raises(ValueError, match="occlusion"):
            tunnel_silhouette(np.ones((10, 10)), np.ones((10, 10)), occlusion=np.ones((8, 8)))


class TestFindMouths:
    """Upstream and downstream are told apart by grey level, by construction."""

    def setup_method(self):
        self.mouth = np.ones((H, W))
        self.mouth = np.minimum(self.mouth, _disc((H, W), 60, 100, 6, 0.0))    # upstream, black
        self.mouth = np.minimum(self.mouth, _disc((H, W), 240, 100, 6, 0.5))   # downstream, grey
        self.sil = _bar((H, W), 60, 240, 96, 104, value=True, base=False).astype(bool)

    def test_the_black_disc_is_upstream_and_the_grey_one_downstream(self):
        m = find_mouths(self.mouth, self.sil)
        assert set(m) == {"up", "dn"}
        assert m["up"].centre[0] == pytest.approx(60, abs=1.5)
        assert m["dn"].centre[0] == pytest.approx(240, abs=1.5)

    def test_forward_points_out_of_the_protein(self):
        m = find_mouths(self.mouth, self.sil)
        # the channel lies to the RIGHT of the upstream exit, so out is left
        assert m["up"].forward[0] < -0.8
        assert m["dn"].forward[0] > 0.8

    def test_the_radius_is_the_disc_radius(self):
        m = find_mouths(self.mouth, self.sil)
        assert m["up"].radius == pytest.approx(6.0, abs=0.7)

    def test_a_missing_exit_is_an_error_not_a_half_figure(self):
        only_up = np.minimum(np.ones((H, W)), _disc((H, W), 60, 100, 6, 0.0))
        with pytest.raises(ValueError, match="no dn exit"):
            find_mouths(only_up, self.sil)

    def test_specks_are_ignored(self):
        noisy = np.minimum(self.mouth, _disc((H, W), 150, 30, 1, 0.0))
        m = find_mouths(noisy, self.sil, min_pixels=20)
        assert len(m) == 2

    def test_the_largest_disc_of_a_class_wins(self):
        # a assembly render can show a second copy of an exit; the figure
        # annotates one channel, so the bigger (nearer) disc is the subject
        two = np.minimum(self.mouth, _disc((H, W), 150, 160, 3, 0.0))
        m = find_mouths(two, self.sil)
        assert m["up"].centre[1] == pytest.approx(100, abs=2), "the big black disc, not the speck"


class TestTraceContours:
    def test_one_outline_per_channel(self):
        sil = np.zeros((H, W), bool)
        sil[90:110, 20:120] = True
        sil[90:110, 180:280] = True
        polys = trace_contours(sil, smooth=0.6, min_area=50)
        assert len(polys) == 2

    def test_specks_below_min_area_are_dropped(self):
        sil = np.zeros((H, W), bool)
        sil[90:110, 20:120] = True
        sil[10:13, 10:13] = True
        assert len(trace_contours(sil, smooth=0.6, min_area=50)) == 1

    def test_contours_are_closed_loops_in_image_coordinates(self):
        sil = np.zeros((H, W), bool)
        sil[80:120, 100:200] = True
        (p,) = trace_contours(sil, smooth=0.6, min_area=50)
        assert np.allclose(p[0], p[-1], atol=1.5), "a traced silhouette is a closed loop"
        assert 70 < p[:, 1].min() < 90 and 110 < p[:, 1].max() < 130


class TestOpenRuns:
    """The outline is cut open at the mouth: a channel is not a closed blob."""

    @staticmethod
    def rect(x0, x1, y0, y1, n=40):
        t = np.linspace(0, 1, n)
        top = np.stack([x0 + (x1 - x0) * t, np.full(n, y0)], 1)
        right = np.stack([np.full(n, x1), y0 + (y1 - y0) * t], 1)
        bot = np.stack([x1 + (x0 - x1) * t, np.full(n, y1)], 1)
        left = np.stack([np.full(n, x0), y1 + (y0 - y1) * t], 1)
        return np.vstack([top, right, bot, left])

    def test_a_contour_far_from_every_mouth_stays_closed(self):
        p = self.rect(10, 60, 10, 40)
        m = Mouth("up", np.array([500.0, 500.0]), 6.0, np.array([1.0, 0.0]))
        (run,) = open_runs(p, [m])
        assert np.allclose(run[0], run[-1]), "closed: the loop is explicitly re-joined"

    def test_the_outline_opens_where_it_crosses_a_mouth(self):
        p = self.rect(10, 110, 80, 120)
        m = Mouth("up", np.array([110.0, 100.0]), 10.0, np.array([1.0, 0.0]))
        runs = open_runs(p, [m])
        assert runs, "something must survive"
        pts = np.vstack(runs)
        d = np.linalg.norm(pts - m.centre, axis=1)
        assert d.min() > 1.0, "no point is left sitting on the exit disc"
        assert not np.allclose(runs[0][0], runs[0][-1]), "the run is open"

    def test_a_run_straddling_the_start_of_the_contour_is_not_split(self):
        # the erased stretch is placed at the contour's first point: without
        # the roll, the single surviving run would come back as two
        p = self.rect(10, 110, 80, 120)
        m = Mouth("up", np.array([10.0, 80.0]), 12.0, np.array([-1.0, -1.0]) / math.sqrt(2))
        runs = open_runs(p, [m])
        assert len(runs) == 1, "one opening, one run"

    def test_only_the_outward_side_of_the_disc_is_erased(self):
        p = self.rect(10, 110, 80, 120)
        far = Mouth("up", np.array([60.0, 100.0]), 12.0, np.array([0.0, -1.0]))
        runs = open_runs(p, [far])
        pts = np.vstack(runs)
        assert (pts[:, 1] > 110).any(), "the bottom edge, behind the mouth, must survive"

    def test_a_contour_entirely_inside_a_mouth_disappears(self):
        p = self.rect(98, 102, 98, 102)
        m = Mouth("up", np.array([100.0, 100.0]), 30.0, np.array([1.0, 0.0]))
        assert open_runs(p, [m]) == []


class TestPaths:
    def test_closed_path_ends_with_z(self):
        d = path_closed(np.array([[0.0, 0], [10, 0], [10, 10]]), step=1)
        assert d.startswith("M0.0,0.0") and d.endswith(" Z")

    def test_open_path_does_not(self):
        d = path_open(np.array([[0.0, 0], [10, 0]]), step=1)
        assert "Z" not in d

    def test_the_step_thins_the_point_list(self):
        p = np.zeros((20, 2))
        assert path_open(p, step=4).count("L") == 4


class TestPlaceBlock:
    """Blocks sit clear of the line, on the side asked for, over empty space."""

    def setup_method(self):
        self.alpha = np.zeros((400, 600), bool)
        self.alpha[150:250, 100:500] = True          # protein across the middle
        self.up = np.array([120.0, 200.0])
        self.dn = np.array([480.0, 200.0])

    def test_above_is_above_and_below_is_below(self):
        a = place_block(self.alpha, self.up, self.dn, 120, 40, True, clearance=10)
        b = place_block(self.alpha, self.up, self.dn, 120, 40, False, clearance=10)
        assert a.y + 40 < 200 - 10
        assert b.y > 200 + 10

    def test_it_prefers_empty_space(self):
        p = place_block(self.alpha, self.up, self.dn, 120, 40, True, clearance=10)
        assert p.overlap < 0.2, "a block should not be parked on the protein"

    def test_a_vertical_line_still_has_two_sides(self):
        # the earlier slope test divides by zero here; the normal-based one does not
        up, dn = np.array([300.0, 60.0]), np.array([300.0, 340.0])
        a = place_block(self.alpha, up, dn, 100, 30, True, clearance=10)
        b = place_block(self.alpha, up, dn, 100, 30, False, clearance=10)
        assert a.x != b.x, "the two blocks must land on opposite sides, not on top"
        assert (a.x + 100 < 300) != (b.x + 100 < 300)

    def test_the_two_exits_at_one_point_is_an_error(self):
        with pytest.raises(ValueError, match="same point"):
            place_block(self.alpha, self.up, self.up, 10, 10, True)

    def test_an_impossible_clearance_still_returns_a_position(self):
        # better a label in a poor place than a figure with no numbers
        p = place_block(self.alpha, self.up, self.dn, 120, 40, True, clearance=10_000)
        assert p.overlap == 1.0
        assert np.isfinite([p.x, p.y]).all()


class TestAngleMarker:
    def test_the_arc_spans_the_reported_angle(self):
        up, dn = np.array([100.0, 100.0]), np.array([300.0, 100.0])
        m = angle_marker(up, dn, 150.0, radius=50, arm=80)
        a = m.arc_start - m.vertex
        b = m.arc_end - m.vertex
        cos = float(a @ b / (np.linalg.norm(a) * np.linalg.norm(b)))
        assert math.degrees(math.acos(max(-1.0, min(1.0, cos)))) == pytest.approx(150.0, abs=0.1)

    def test_the_arm_points_along_the_upstream_axis_and_the_arrow_ends_at_the_exit(self):
        up, dn = np.array([100.0, 100.0]), np.array([300.0, 100.0])
        m = angle_marker(up, dn, 150.0, radius=50, arm=80)
        assert np.allclose(m.arrowhead[0], up), "the arrow points AT P2"
        assert np.linalg.norm(m.arm_start - up) == pytest.approx(80.0)

    def test_a_measured_axis_overrides_the_analytic_rule(self):
        up, dn = np.array([100.0, 100.0]), np.array([300.0, 100.0])
        m = angle_marker(up, dn, 150.0, radius=50, arm=80, axis_px=[0.0, -1.0])
        assert m.exact
        assert m.arm_start[1] == pytest.approx(20.0), "straight up, as given"
        assert m.arm_start[0] == pytest.approx(100.0)

    def test_the_analytic_marker_says_it_is_derived(self):
        m = angle_marker([0.0, 0], [10.0, 0], 90.0, 5, 5)
        assert not m.exact

    def test_a_zero_axis_is_refused(self):
        with pytest.raises(ValueError, match="zero length"):
            angle_marker([0.0, 0], [10.0, 0], 90.0, 5, 5, axis_px=[0, 0])

    def test_the_arc_is_drawn_the_short_way_round(self):
        up, dn = np.array([100.0, 100.0]), np.array([300.0, 100.0])
        for oic in (30.0, 90.0, 150.0, 179.0):
            m = angle_marker(up, dn, oic, radius=50, arm=80)
            large, _ = _arc_flags(m)
            assert large == 0, "O_IC below 180 deg never needs the large arc"

    def test_an_obtuse_turn_the_other_way_sets_the_sweep_flag(self):
        up, dn = np.array([100.0, 100.0]), np.array([300.0, 100.0])
        a = _arc_flags(angle_marker(up, dn, 150.0, 50, 80))
        b = _arc_flags(angle_marker(up, dn, -150.0, 50, 80))
        assert a[1] != b[1], "the two turning directions must not draw the same arc"


class TestSignedDelta:
    def test_wraps_to_the_short_way(self):
        assert _signed_delta(0.0, math.radians(350)) == pytest.approx(math.radians(-10))
        assert _signed_delta(0.0, math.radians(10)) == pytest.approx(math.radians(10))

    def test_half_turn_is_positive_pi(self):
        assert _signed_delta(0.0, math.pi) == pytest.approx(math.pi)


class TestContentBox:
    def test_it_crops_to_the_drawn_content(self):
        alpha = np.zeros((400, 600), bool)
        alpha[100:200, 200:300] = True
        sil = np.zeros_like(alpha)
        x, y, w, h = content_box(alpha, sil, margin=0.0)
        assert (x, y) == (200.0, 100.0)
        assert w == pytest.approx(99.0) and h == pytest.approx(99.0)

    def test_extra_points_pull_the_frame_out(self):
        alpha = np.zeros((400, 600), bool)
        alpha[100:200, 200:300] = True
        x, _, w, _ = content_box(alpha, np.zeros_like(alpha), [(20.0, 150.0)], margin=0.0)
        assert x == 20.0 and w == pytest.approx(279.0)

    def test_the_margin_is_a_fraction_of_the_larger_side(self):
        alpha = np.zeros((400, 600), bool)
        alpha[100:200, 100:400] = True
        _, _, w, _ = content_box(alpha, np.zeros_like(alpha), margin=0.1)
        assert w == pytest.approx(299.0 + 2 * 0.1 * 299.0, rel=0.02)

    def test_the_frame_never_leaves_the_render(self):
        alpha = np.ones((40, 60), bool)
        x, y, w, h = content_box(alpha, np.zeros_like(alpha), margin=0.5)
        assert x >= 0 and y >= 0 and x + w <= 60 and y + h <= 40

    def test_an_empty_render_gives_the_whole_frame_rather_than_crashing(self):
        empty = np.zeros((40, 60), bool)
        assert content_box(empty, empty, margin=0.0) == (0.0, 0.0, 60.0, 40.0)


class TestSvgStyle:
    def test_defaults_match_the_published_figure(self):
        s = SvgStyle()
        assert s.fill_opacity == 0.30
        assert s.labels == ["D_IC", "O_IC", "DO"]

    def test_overrides_from_a_pair_file(self):
        s = SvgStyle.from_dict({"fill_opacity": 0.5, "font_size": 60, "labels": "DO"})
        assert s.fill_opacity == 0.5 and s.font_size == 60.0
        assert s.labels == ["DO"]

    def test_an_unknown_key_is_an_error_not_a_silent_default(self):
        with pytest.raises(ValueError, match="unknown key"):
            SvgStyle.from_dict({"fill_opactiy": 0.5})

    def test_an_unknown_label_is_an_error(self):
        with pytest.raises(ValueError, match="unknown label"):
            SvgStyle.from_dict({"labels": ["D_IC", "RMSD"]})

    def test_out_of_range_opacity_is_refused(self):
        with pytest.raises(ValueError, match=r"\[0, 1\]"):
            SvgStyle(fill_opacity=1.4)

    def test_a_nonpositive_size_is_refused(self):
        with pytest.raises(ValueError, match="positive"):
            SvgStyle(font_size=0.0)

    def test_both_spellings_of_colour_at_once_are_refused(self):
        with pytest.raises(ValueError, match="only one"):
            SvgStyle.from_dict({"label_color": "black", "label_colour": "black"})

    def test_the_british_spelling_works(self):
        assert SvgStyle.from_dict({"label_colour": "crimson"}).label_color == "crimson"

    def test_it_rides_along_in_the_render_section(self):
        st = RenderStyle.from_dict({"svg": {"fill_opacity": 0.4}})
        assert SvgStyle.from_dict(st.svg).fill_opacity == 0.4


class TestDocument:
    """The written SVG: well-formed, ids present, colours from the style."""

    def _layout(self, **kw):
        alpha = np.zeros((300, 400), bool)
        alpha[100:200, 50:350] = True
        sil = np.zeros_like(alpha)
        sil[140:160, 80:320] = True
        poly = np.array([[80.0, 140], [320, 140], [320, 160], [80, 160], [80, 140]])
        mouths = {
            "up": Mouth("up", np.array([80.0, 150.0]), 8.0, np.array([-1.0, 0.0])),
            "dn": Mouth("dn", np.array([320.0, 150.0]), 8.0, np.array([1.0, 0.0])),
        }
        marker = angle_marker(mouths["up"].centre, mouths["dn"].centre, 165.0, 40, 60)
        base = dict(
            width=400, height=300, scale=0.2, mouths=mouths, polys=[poly],
            silhouette=sil, alpha=alpha, marker=marker,
            blocks={"above": __import__("channelbridge.render.svg", fromlist=["Placement"])
                    .Placement(100.0, 60.0, 0.0),
                    "below": __import__("channelbridge.render.svg", fromlist=["Placement"])
                    .Placement(140.0, 220.0, 0.0)},
            box=(0.0, 0.0, 400.0, 300.0),
        )
        base.update(kw)
        return _Layout(**base)

    VALUES: ClassVar[dict] = {"D_IC": 17.0, "O_IC": 165.0, "DO": 0.850}

    def doc(self, style=None, svg_style=None, **kw):
        return svg_document(self._layout(**kw), "AAAA", self.VALUES,
                            style or DEFAULT_STYLE, svg_style or SvgStyle())

    def test_it_is_well_formed_xml(self):
        xml.dom.minidom.parseString(self.doc())

    def test_every_layer_has_an_id(self):
        d = self.doc()
        for want in ('id="background"', 'id="proteins"', 'id="tunnels"',
                     'id="tunnel_fill_0"', 'id="channel_line"', 'id="angle"',
                     'id="angle_arc"', 'id="label_DIC"', 'id="label_OIC"', 'id="label_DO"'):
            assert want in d, want

    def test_the_background_is_transparent_by_default(self):
        assert 'id="background"' in self.doc() and 'fill-opacity="0"' in self.doc()

    def test_the_numbers_are_real_text_with_real_subscripts(self):
        d = self.doc()
        assert "17.0" in d and "165.0" in d and "0.850" in d
        assert "<tspan>D</tspan>" in d and ">IC</tspan>" in d

    def test_the_tunnel_takes_its_colour_from_the_palette(self):
        st = RenderStyle.from_dict({"colors": {"tunnel": "#1133ff"}})
        assert "#1133ff" in self.doc(style=st)

    def test_the_fill_opacity_is_the_styles(self):
        assert 'fill-opacity="0.45"' in self.doc(svg_style=SvgStyle(fill_opacity=0.45))

    def test_dropping_the_labels_leaves_the_geometry(self):
        d = self.doc(svg_style=SvgStyle(labels=[]), blocks={})
        assert 'id="channel_line"' in d
        assert 'id="label_DO"' not in d
        xml.dom.minidom.parseString(d)

    def test_without_a_marker_nothing_angular_is_drawn(self):
        d = self.doc(marker=None)
        assert 'id="angle"' not in d and 'id="channel_line"' in d
        xml.dom.minidom.parseString(d)

    def test_the_viewbox_is_the_layouts_box(self):
        d = self.doc(box=(10.0, 20.0, 100.0, 80.0))
        assert 'viewBox="10.0 20.0 100.0 80.0"' in d

    def test_the_raster_keeps_the_renders_own_pixel_frame(self):
        # the crop lives in the viewBox, so the overlay's coordinates stay the
        # render's coordinates and can be checked against the PNG
        d = self.doc(box=(10.0, 20.0, 100.0, 80.0))
        assert '<image id="proteins" x="0" y="0" width="400" height="300"' in d

    def test_a_title_is_escaped(self):
        d = svg_document(self._layout(), "AAAA", self.VALUES, DEFAULT_STYLE, SvgStyle(),
                         title="GDH & FNR <geom001>")
        assert "&amp;" in d and "&lt;geom001&gt;" in d
        xml.dom.minidom.parseString(d)


class TestMarkerDataclass:
    def test_it_reports_whether_the_arm_was_measured(self):
        m = AngleMarker(np.zeros(2), np.zeros(2), np.zeros((3, 2)), np.zeros(2),
                        np.zeros(2), np.zeros(2), 1.0, True)
        assert m.exact


class TestTheShippedExample:
    """The example pair file must actually load: a broken example misleads."""

    def test_the_render_section_produces_both_styles(self):
        import os

        from channelbridge.config import load_pair_file
        from channelbridge.render.hotspots import ClassPalette

        path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                            "examples", "gdh_fnr", "pair.yaml")
        if not os.path.exists(path):
            pytest.skip("example not present")
        style = load_pair_file(path).render
        assert SvgStyle.from_dict(style.svg).fill_opacity == 0.30
        assert ClassPalette.from_dict(style.hotspots).classes["CHANNEL"] == "magenta"
        assert style.transparency_for("inspect") == 0.55


class TestFocusOutline:
    """The two protomers the geometry belongs to, picked out of an assembly.

    In the assembly render every copy of both enzymes is drawn and the pose's
    own pair is distinguished only by being a shade darker -- behind two
    translucent neighbours that is not much. The outline says it outright.
    It is traced from the `footprint` render, which holds exactly those two
    chains, so no extra scene is needed.
    """

    def _frames(self, w=200, h=160):
        import numpy as np

        footprint = np.ones((h, w))          # white = nothing
        footprint[40:120, 40:120] = 0.0      # the focus pair, flat black
        return footprint

    def test_nothing_is_drawn_on_a_pair_render(self):
        # there the footprint IS the whole picture, so an outline would just
        # trace the outside of the frame
        import numpy as np

        from channelbridge.render.svg import focus_outline

        fp = self._frames()
        alpha = fp < 0.5                      # only the pair is drawn
        assert focus_outline(fp, alpha, min_area=10.0) == []
        assert isinstance(np.asarray(alpha), np.ndarray)

    def test_it_is_drawn_when_other_protomers_are_present(self):
        from channelbridge.render.svg import focus_outline

        fp = self._frames()
        alpha = fp < 0.5
        alpha[10:30, 10:30] = True            # a neighbouring copy
        alpha[130:150, 130:150] = True        # and another
        polys = focus_outline(fp, alpha, min_area=10.0)
        assert polys, "an assembly render should get an outline"

    def test_mismatched_renders_are_refused(self):
        import numpy as np
        import pytest

        from channelbridge.render.svg import focus_outline

        with pytest.raises(ValueError, match="different renders"):
            focus_outline(self._frames(), np.zeros((10, 10), bool))

    def test_a_zero_width_switches_it_off(self):
        from channelbridge.render.svg import SvgStyle

        assert SvgStyle.from_dict({"focus_stroke_width": 0}).focus_stroke_width == 0.0


class TestMouthRing:
    """An exit that does not reach the silhouette's rim gets a drawn opening.

    On the worked example's pose the downstream channel leaves at 58 degrees to the view
    plane -- it points at the camera -- so its exit projects 31 px INSIDE the
    traced outline. There is no gap to cut, and without a ring the D_IC line
    appears to stop in the middle of the channel for no reason. The upstream
    channel leaves at 11 degrees and sits exactly on the rim, where the gap
    is the cue and a ring would be clutter.
    """

    def _mouth_render(self, w=200, h=160, up=(60, 80), dn=(140, 80),
                      r_up=6, r_dn=6):
        import numpy as np

        img = np.ones((h, w))
        yy, xx = np.mgrid[0:h, 0:w]
        img[(xx - up[0]) ** 2 + (yy - up[1]) ** 2 <= r_up ** 2] = 0.0    # black
        img[(xx - dn[0]) ** 2 + (yy - dn[1]) ** 2 <= r_dn ** 2] = 0.5    # grey
        return img

    def test_the_projected_ellipse_is_measured_from_the_render(self):
        import numpy as np

        from channelbridge.render.svg import find_mouths

        sil = np.zeros((160, 200), bool)
        sil[60:100, 100:180] = True          # silhouette around the dn disc
        mouths = find_mouths(self._mouth_render(), sil)
        dn = mouths["dn"]
        assert dn.semi_major > 0 and dn.semi_minor > 0
        # a round disc projects round: the axes should be close
        assert abs(dn.semi_major - dn.semi_minor) < 0.25 * dn.semi_major

    def test_inset_is_zero_on_the_rim_and_positive_inside(self):
        import numpy as np

        from channelbridge.render.svg import find_mouths

        sil = np.zeros((160, 200), bool)
        sil[60:100, 100:180] = True          # dn sits inside this
        mouths = find_mouths(self._mouth_render(), sil)
        assert mouths["dn"].inset_px > 0, "inside the silhouette"
        assert mouths["up"].inset_px == 0, "nowhere near it"

    def test_a_zero_width_switches_it_off(self):
        from channelbridge.render.svg import SvgStyle

        assert SvgStyle.from_dict({"mouth_ring_width": 0}).mouth_ring_width == 0.0


class TestTheArmMustBeMeasuredOnARolledFrame:
    """`angle_marker` without `axis_px` assumes the P1-P2-P3 plane is face-on.

    That holds in the angle frame and nowhere else. On the worked example's pose rolled
    160 degrees about the arm, the fallback drew the arm at -151 deg where
    the session has it at -180 -- 29 deg of pure invention, while still
    showing the correct 165 deg BETWEEN arm and line, which is what made it
    look right.
    """

    def test_the_fallback_says_it_is_a_fallback(self):
        from channelbridge.render.svg import angle_marker

        m = angle_marker((100.0, 100.0), (200.0, 120.0), 165.0, 110.0, 150.0)
        assert not m.exact, "an assumed arm must be flagged as assumed"

    def test_a_given_axis_is_used_verbatim(self):
        import numpy as np

        from channelbridge.render.svg import angle_marker

        m = angle_marker((100.0, 100.0), (200.0, 120.0), 165.0, 110.0, 150.0,
                         axis_px=(-1.0, 0.0))
        assert m.exact
        arm = np.asarray(m.arm_start, float) - np.asarray(m.vertex, float)
        # horizontal, pointing the way we said
        assert abs(float(arm[1])) < 1e-6, "a horizontal arm must draw horizontal"
        assert float(arm[0]) < 0.0

    def test_a_zero_axis_is_refused_not_normalised_to_nothing(self):
        import pytest

        from channelbridge.render.svg import angle_marker

        with pytest.raises(ValueError, match="zero length"):
            angle_marker((100.0, 100.0), (200.0, 120.0), 165.0, 110.0, 150.0,
                         axis_px=(0.0, 0.0))


class TestTheChannelAndTheMeasurementReadTheirOwnRoles:
    """`tunnel` colours the channel, `annotation` the D_IC/O_IC marks.

    They hold the same magenta by design, so these assertions are about
    which ROLE each element consults, not about the values. That matters:
    while the two briefly differed, this caught a call site that had mixed
    them up for as long as they were both red -- the PNG preview drew the
    channel OUTLINE with the annotation colour. A shared value hides that
    class of bug, so the roles stay separate and are checked separately.
    """

    def _doc(self):
        return TestDocument().doc()

    def test_the_channel_is_the_tunnel_colour(self):
        d = self._doc()
        tunnel = rgb_to_hex(DEFAULT_STYLE.palette.rgb("tunnel"))
        fill = re.search(r'id="tunnel_fill_0"[^/]*fill="([^"]+)"', d)
        assert fill and fill.group(1) == tunnel
        # this layout's two mouths may consume the whole outline, so the
        # outline runs are optional -- but any that IS drawn is the channel
        outlines = re.findall(r'id="tunnel_outline_\d+_\d+"[^/]*stroke="([^"]+)"', d)
        assert all(c == tunnel for c in outlines)

    def test_the_measurement_is_the_annotation_colour(self):
        d = self._doc()
        annot = rgb_to_hex(DEFAULT_STYLE.palette.rgb("annotation"))
        for pat in (r'id="channel_line"[^/]*stroke="([^"]+)"',
                    r'id="exit_upstream"[^/]*fill="([^"]+)"',
                    r'id="exit_downstream"[^/]*fill="([^"]+)"'):
            m = re.search(pat, d)
            assert m, pat
            assert m.group(1) == annot, pat

    def test_an_override_of_one_does_not_move_the_other(self):
        """The proof that these are two roles and not one."""
        from channelbridge.render.style import Palette

        style = RenderStyle(palette=Palette().with_overrides(
            {"annotation": (0.0, 0.0, 0.0)}))
        d = TestDocument().doc(style=style)
        assert re.search(r'id="channel_line"[^/]*stroke="#000000"', d)
        tunnel = rgb_to_hex(style.palette.rgb("tunnel"))
        fill = re.search(r'id="tunnel_fill_0"[^/]*fill="([^"]+)"', d)
        assert fill and fill.group(1) == tunnel != "#000000"
