from protopo.protopo import ProTopo


def test_add_alpha_records_index_map():
    pt = ProTopo()
    pt.add_alpha(10, 13, label="H", to="→", scale=1.0, add_label=False)

    # 检查 index_map 中是否记录了三个位置
    for i in range(10, 13):
        assert i in pt.index_map
        assert pt.index_map[i]["type"] == "alpha"
        assert pt.index_map[i]["to"] == "→"


def test_add_beta_records_index_map():
    pt = ProTopo()
    pt.add_beta(20, 23, label="S", to="↓", scale=1.0, add_label=False)

    for i in range(20, 23):
        assert i in pt.index_map
        assert pt.index_map[i]["type"] == "beta"
        assert pt.index_map[i]["to"] == "↓"


def test_add_linker_steps_and_directions():
    pt = ProTopo()
    pt.add_linker(1, 5, to="→↓", steps=(0.5, 0.5), scale=1.0)

    # 检查 index_map 是否分段记录
    assert 1 in pt.index_map
    assert 2 in pt.index_map
    assert 3 in pt.index_map
    assert 4 in pt.index_map
    # 总共 4 residues，应分两段

    count_right = sum(1 for v in pt.index_map.values() if v["to"] == "→")
    count_down = sum(1 for v in pt.index_map.values() if v["to"] == "↓")

    assert count_right + count_down == 4
    assert count_right in (2, 3)  # 因为四个残基无法正好均分
    assert count_down in (1, 2)


def test_marker_addition_does_not_fail():
    pt = ProTopo()
    pt.add_beta(100, 105)
    pt.add_triangles([101, 104])  # 应该不报错
