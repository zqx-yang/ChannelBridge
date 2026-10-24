"""The figure selection language and its dependency rules."""

from __future__ import annotations

import pytest

from channelbridge.render.select import (
    VARIANTS,
    parse_spec,
    parse_variant,
    rank_poses,
    resolve,
)

# deliberately out of order, with a tie, so ranking is actually exercised
POSES = [
    ("geom003", 0.70),
    ("geom001", 0.85),
    ("geom004", 0.70),
    ("geom002", 0.80),
    ("geom005", 0.60),
]


@pytest.fixture
def ranked():
    return rank_poses(POSES)


class TestRanking:
    def test_ranks_by_do_descending(self, ranked):
        assert [p.pose_id for p in ranked] == [
            "geom001", "geom002", "geom003", "geom004", "geom005"]
        assert [p.rank for p in ranked] == [1, 2, 3, 4, 5]

    def test_ties_break_on_id_so_rank_is_stable(self):
        a = rank_poses([("z", 0.5), ("a", 0.5)])
        b = rank_poses([("a", 0.5), ("z", 0.5)])
        assert [p.pose_id for p in a] == [p.pose_id for p in b] == ["a", "z"]

    def test_accepts_pose_objects(self, ranked):
        assert rank_poses(ranked) == ranked

    def test_duplicate_ids_are_rejected(self):
        with pytest.raises(ValueError, match="duplicate pose id"):
            rank_poses([("geom001", 0.8), ("geom001", 0.7)])


class TestSelectionLanguage:
    def test_none_and_empty(self, ranked):
        for spec in (None, "", "none", "NONE", " "):
            assert parse_spec(spec, ranked) == []

    def test_all(self, ranked):
        assert len(parse_spec("all", ranked)) == 5

    def test_a_bare_count_is_top_n(self, ranked):
        assert parse_spec("3", ranked) == ["geom001", "geom002", "geom003"]
        assert parse_spec(3, ranked) == ["geom001", "geom002", "geom003"]

    def test_topn_spellings(self, ranked):
        for spec in ("top3", "TOP3", "top:3", "top 3"):
            assert parse_spec(spec, ranked) == ["geom001", "geom002", "geom003"]

    def test_a_single_rank(self, ranked):
        assert parse_spec("rank3", ranked) == ["geom003"]
        assert parse_spec("#3", ranked) == ["geom003"]

    def test_a_rank_range(self, ranked):
        assert parse_spec("rank2-4", ranked) == ["geom002", "geom003", "geom004"]
        assert parse_spec("#2-4", ranked) == ["geom002", "geom003", "geom004"]

    def test_by_pose_id(self, ranked):
        assert parse_spec("geom004", ranked) == ["geom004"]

    def test_terms_combine(self, ranked):
        assert parse_spec("top2,geom005,rank3", ranked) == [
            "geom001", "geom002", "geom003", "geom005"]

    def test_result_is_in_rank_order_not_typed_order(self, ranked):
        # outputs should be produced best-first whatever the user wrote
        assert parse_spec("geom005,geom001", ranked) == ["geom001", "geom005"]

    def test_duplicates_collapse(self, ranked):
        assert parse_spec("top3,geom001,rank2", ranked) == [
            "geom001", "geom002", "geom003"]

    def test_whitespace_is_tolerated(self, ranked):
        assert parse_spec("  top2 , geom005  ", ranked) == [
            "geom001", "geom002", "geom005"]

    def test_a_count_beyond_the_end_is_not_an_error(self, ranked):
        assert len(parse_spec("99", ranked)) == 5

    def test_a_rank_beyond_the_end_is_an_error(self, ranked):
        # asking for rank 9 of 5 is a mistake worth catching, unlike top9
        with pytest.raises(ValueError, match="only 5 pose"):
            parse_spec("rank9", ranked)

    def test_a_rank_range_beyond_the_end_clips_from_a_valid_start(self, ranked):
        assert parse_spec("rank4-99", ranked) == ["geom004", "geom005"]

    def test_rank_zero_is_an_error(self, ranked):
        with pytest.raises(ValueError, match="start at 1"):
            parse_spec("rank0", ranked)

    def test_an_inverted_range_is_an_error(self, ranked):
        with pytest.raises(ValueError, match="rank range"):
            parse_spec("rank4-2", ranked)

    def test_an_unknown_id_is_an_error_and_suggests(self, ranked):
        with pytest.raises(ValueError, match="did you mean"):
            parse_spec("geom00", ranked)

    def test_gibberish_explains_the_language(self, ranked):
        with pytest.raises(ValueError, match="Use a count"):
            parse_spec("best ones please", ranked)

    def test_the_error_names_which_option_was_wrong(self, ranked):
        with pytest.raises(ValueError, match=r"--svg"):
            parse_spec("nonsense", ranked, "--svg")


class TestVariants:
    def test_the_four_values(self):
        for v in VARIANTS:
            assert parse_variant(v) == v

    def test_none_by_default(self):
        assert parse_variant(None) == "none"
        assert parse_variant("") == "none"

    def test_aliases(self):
        assert parse_variant("mono") == "monomer"
        assert parse_variant("multi") == "multimer"
        assert parse_variant("complex") == "monomer"
        assert parse_variant("occupancy") == "multimer"
        assert parse_variant("all") == "both"

    def test_case_insensitive(self):
        assert parse_variant("Monomer") == "monomer"

    def test_an_unknown_value_lists_the_options(self):
        with pytest.raises(ValueError, match="monomer"):
            parse_variant("solid")


class TestDependencyAndPlan:
    def test_an_svg_pulls_in_its_session(self):
        plan = resolve(POSES, svg="geom004", render="monomer")
        assert plan.pse == ["geom004"]
        assert plan.pse_added_for_svg == ["geom004"]
        assert "added because an SVG needs one" in plan.summary()

    def test_an_existing_session_is_not_rebuilt_for_an_svg(self):
        plan = resolve(POSES, svg="geom004", render="monomer",
                       existing_pse=["geom004"])
        assert plan.pse == []
        assert plan.pse_reused == ["geom004"]
        assert "already on disk, not rebuilt" in plan.summary()

    def test_rebuild_forces_it_anyway(self):
        plan = resolve(POSES, svg="geom004", render="monomer",
                       existing_pse=["geom004"], rebuild=True)
        assert plan.pse == ["geom004"]

    def test_a_session_already_requested_is_not_double_counted(self):
        plan = resolve(POSES, pse="geom004", svg="geom004", render="monomer")
        assert plan.pse == ["geom004"]
        assert plan.pse_added_for_svg == []

    def test_the_three_outputs_are_independent(self):
        plan = resolve(POSES, pse="top3", svg="rank1", pdb="all",
                       render="monomer")
        assert plan.pse == ["geom001", "geom002", "geom003"]
        assert plan.svg == ["geom001"]
        assert len(plan.pdb) == 5

    def test_pse_list_stays_in_rank_order_after_dependency_merge(self):
        plan = resolve(POSES, pse="geom005", svg="geom001", render="monomer")
        assert plan.pse == ["geom001", "geom005"]

    def test_an_svg_with_no_drawable_scene_is_flagged(self):
        # the SVG is traced from render / cross_section; with both off there
        # is nothing to trace
        plan = resolve(POSES, svg="rank1", render="none", cross_section="none")
        assert any("set at least one" in n for n in plan.notes)

    def test_no_flag_when_a_scene_is_on(self):
        plan = resolve(POSES, svg="rank1", cross_section="monomer")
        assert not any("set at least one" in n for n in plan.notes)

    def test_overwriting_an_existing_session_is_announced(self):
        plan = resolve(POSES, pse="geom001", existing_pse=["geom001"])
        assert any("will be overwritten" in n for n in plan.notes)

    def test_nothing_selected_is_reported_helpfully(self):
        plan = resolve(POSES)
        assert plan.nothing_to_do
        assert "nothing selected" in plan.summary()

    def test_summary_shows_ranks(self):
        plan = resolve(POSES, pse="geom003")
        assert "geom003 (rank 3)" in plan.summary()

    def test_summary_shows_the_scene_variants(self):
        plan = resolve(POSES, pse="rank1", render="both", cross_section="monomer")
        s = plan.summary()
        assert "render=both" in s and "cross_section=monomer" in s
        assert "classes, occupancy" in s, "the always-present scenes should be named"

    def test_rank_of_resolves(self):
        plan = resolve(POSES, pse="all")
        assert plan.rank_of("geom001") == 1
        assert plan.rank_of("geom005") == 5
        with pytest.raises(KeyError):
            plan.rank_of("nope")

    def test_a_bad_variant_fails_before_any_work(self):
        with pytest.raises(ValueError, match="render"):
            resolve(POSES, pse="rank1", render="wireframe")
