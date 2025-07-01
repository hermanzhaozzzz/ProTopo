from protopo import ProTopo


def draw_egfp(figsize=(8, 6)):
    pt = ProTopo(figsize=figsize)
    pt.add_n_term()  # a N term text marker
    pt.add_linker(1, 3, to="←", scale=0.5)  # [1, 3) -> 1, 2
    pt.add_alpha(3, 9, label="H1", to="←", scale=0.5)
    pt.add_linker(9, 12, to="←↑", steps=(0.1, 0.9))
    pt.add_beta(12, 23, label="S1", to="↑")
    pt.add_linker(23, 25, to="↑→↓", steps=(0.2, 0.6, 0.2), scale=1.2)
    pt.add_beta(25, 36, label="S2", to="↓")
    pt.add_linker(36, 40, to="↓→↑", steps=(0.2, 0.6, 0.2))
    pt.add_beta(40, 49, label="S3", to="↑")
    pt.add_linker(49, 56, to="↑→↓", steps=(0.2, 0.6, 0.2), scale=0.5)
    pt.add_alpha(56, 61, label="H2", to="↓")
    pt.add_star(61, 68, to="↓", shape="★", color="gold", size=48, scale=0.3)
    pt.add_alpha(68, 72, label="H2", to="↓")
    pt.add_linker(72, 75, to="↓→", steps=(0.4, 0.6))
    pt.add_alpha(75, 82, label="H3", to="→", scale=0.8)
    pt.add_linker(82, 83, to="→")
    pt.add_alpha(83, 87, label="H4", to="→")
    pt.add_linker(87, 93, to="→↑", steps=(0.4, 0.6))
    pt.add_beta(93, 101, label="S4", to="↑")
    pt.add_linker(101, 104, to="↑→↓", steps=(0.2, 0.6, 0.2))
    pt.add_beta(104, 115, label="S5", to="↓")
    pt.add_linker(115, 118, to="↓→↑", steps=(0.2, 0.6, 0.2))
    pt.add_beta(118, 128, label="S6", to="↑")
    pt.add_linker(128, 148, to="↑←↓", steps=(0.2, 0.6, 0.2), scale=1)
    pt.add_beta(148, 155, label="S7", to="↓")
    pt.add_linker(155, 160, to="↓→↑", steps=(0.2, 0.6, 0.2))
    pt.add_beta(160, 171, label="S8", to="↑")
    pt.add_linker(171, 175, to="↑→↓", steps=(0.2, 0.6, 0.2))
    pt.add_beta(175, 188, label="S9", to="↓", scale=0.8)
    pt.add_linker(188, 199, to="↓←↑", steps=(0.2, 0.6, 0.2), scale=1.1)
    pt.add_beta(199, 208, label="S10", to="↑", scale=1.2)
    pt.add_linker(208, 216, to="↑←↓", steps=(0.2, 0.6, 0.2), scale=0.5)
    pt.add_beta(216, 227, label="S11", to="↓", scale=0.7)
    pt.add_linker(227, 239, to="↓", scale=0.3)
    pt.add_c_term()  # end with K239
    return pt


if __name__ == "__main__":
    egfp = draw_egfp()
    egfp.add_triangles(idx=[6, 211])
    egfp.show()
    # egfp.save("../egfp.svg")
