from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "scripts" / "basic_qxycell_step_by_step.py"


def test_basic_step_by_step_script_covers_full_pipeline():
    source = SCRIPT.read_text(encoding="utf-8")
    compile(source, str(SCRIPT), "exec")

    required_calls = [
        "qxy.check(",
        "qxy.run(",
        "qxy.threshold(",
        "qxy.add_metadata(",
        "qxy.remove_ignore(",
        "qxy.celltype(",
        "qxy.qc(",
        "qxy.plot_stacked_bar(",
        "qxy.plot_spatial(",
        "qxy.plot_annotation_polygons(",
        "qxy.plot_marker_positivity_heatmap(",
        "qxy.plot_marker_intensity_heatmap(",
        "qxy.plot_cell_boundaries(",
        "qxy.cn_knn(",
        "qxy.cn_kmeans(",
        "qxy.cn_name(",
        "qxy.plot_cn_heatmap(",
        "qxy.save(",
    ]
    positions = [source.index(call) for call in required_calls]
    assert positions == sorted(positions)
    assert "qxy.workflow(" not in source
    assert 'if __name__ == "__main__":' in source
