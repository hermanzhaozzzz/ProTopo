import math

import matplotlib.patches as patches
import matplotlib.pyplot as plt


class ProTopo:
    def __init__(self):
        self.fig, self.ax = plt.subplots(figsize=(12, 8))
        self.x, self.y = 0, 0
        self.path = [(self.x, self.y)]
        self.markers = []  # 存储已放置 marker 以避免重叠
        self.index_map = {}  # 记录每个 index 的位置、方向和类型
        self.ax.set_aspect("equal")
        self.ax.axis("off")

    def _direction_vector(self, direction):
        return {"↑": (0, 1), "↓": (0, -1), "→": (1, 0), "←": (-1, 0)}[direction]

    def _perpendicular_direction(self, direction):
        return {"↑": "→", "↓": "→", "→": "↓", "←": "↓"}[direction]

    def _record_index_map(self, start, end, dx, dy, to, structure_type):
        for i in range(start, end):
            self.index_map[i] = {
                "x": self.x + dx * (i - start),
                "y": self.y + dy * (i - start),
                "to": to,
                "type": structure_type,
            }

    def _move(self, direction, length):
        dx, dy = self._direction_vector(direction)
        self.x += dx * length
        self.y += dy * length
        self.path.append((self.x, self.y))

    def _draw_line(self, p1, p2, color="black"):
        self.ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color=color, linewidth=1)

    def add_triangles(
        self,
        idx,
        scale=0.1,
        offset=1,
        facecolor="gray",
        edgecolor="black",
        linewidth=1,
    ):
        for i in idx:
            info = self.index_map[i]
            direction = info["to"]
            info = self.index_map.get(i, {})

            if not info:
                print(f"No element found for index {i}!")
                continue

            x, y = info["x"], info["y"]

            if info["type"] != "linker":
                bias = 0
            else:
                bias = 0.4

            if direction == "↑":
                triangle = plt.Polygon(
                    [
                        [x - offset + bias, y - scale],
                        [x - offset + bias, y + scale],
                        [x - 0.4 * offset + bias, y],
                    ],
                    facecolor=facecolor,
                    edgecolor=edgecolor,
                    linewidth=linewidth,
                )
            elif direction == "↓":
                triangle = plt.Polygon(
                    [
                        [x + offset - bias, y - scale],
                        [x + offset - bias, y + scale],
                        [x + 0.4 * offset - bias, y],
                    ],
                    facecolor=facecolor,
                    edgecolor=edgecolor,
                    linewidth=linewidth,
                )
            elif direction == "←":
                triangle = plt.Polygon(
                    [
                        [x - scale, y + offset - bias],
                        [x + scale, y + offset - bias],
                        [x, y + 0.4 * offset - bias],
                    ],
                    facecolor=facecolor,
                    edgecolor=edgecolor,
                    linewidth=linewidth,
                )
            elif direction == "→":
                triangle = plt.Polygon(
                    [
                        [x - scale, y - offset + bias],
                        [x + scale, y - offset + bias],
                        [x, y - 0.4 * offset + bias],
                    ],
                    facecolor=facecolor,
                    edgecolor=edgecolor,
                    linewidth=linewidth,
                )
            else:
                raise ValueError(f"Invalid direction: {direction}")
            self.ax.add_patch(triangle)

    def _add_label_start(self, text, dx, dy):
        offset = 0.3
        x = self.x + dx * offset
        y = self.y + dy * offset
        self.ax.text(x, y, str(text), ha="center", va="center", fontsize=10)

    def _add_label_end(self, text, dx, dy, length):
        offset = 0.3
        x = self.x + dx * (length - offset)
        y = self.y + dy * (length - offset)
        self.ax.text(x, y, str(text), ha="center", va="center", fontsize=10)

    def _add_label_center(self, text, dx, dy, length):
        cx = self.x + dx * length / 2
        cy = self.y + dy * length / 2
        self.ax.text(cx, cy, text, ha="center", va="center", fontsize=12)

    def add_n_term(self):
        self.ax.text(
            self.x,
            self.y,
            "N",
            ha="center",
            va="center",
            fontsize=12,
            fontweight="bold",
        )

    def add_c_term(self):
        self.ax.text(
            self.x,
            self.y,
            "C",
            ha="center",
            va="center",
            fontsize=12,
            fontweight="bold",
        )

    def add_linker(self, start, end, to="→", steps=(), color="black", scale=1.0):
        total_residues = end - start
        assert total_residues >= 1, "Total residues must be >= 1."
        to = list(to)
        if not steps:
            assert len(to) == 1, "When steps is empty, to should be a single direction."
            steps = [1.0]
        else:
            assert len(steps) == len(to), "steps and to must be of the same length."
            assert len(to) in (2, 3), "Only 2- or 3-turn linkers are supported."
            assert abs(sum(steps) - 1.0) < 1e-6, "steps must sum to 1.0"

        i = start
        for _to, step_ratio in zip(to, steps):
            step_length = math.ceil(step_ratio * total_residues)
            dx, dy = self._direction_vector(_to)
            self._record_index_map(
                i,
                i + step_length,
                dx * scale,
                dy * scale,
                _to,
                "linker",
            )
            p0 = (self.x, self.y)
            self._move(_to, step_length * scale)
            p1 = (self.x, self.y)
            self._draw_line(p0, p1, color=color)
            i += step_length

    def add_alpha(
        self, start, end, label="", to="→", color="brown", scale=1.0, add_label=True
    ):
        length = (end - start) * scale
        dx, dy = self._direction_vector(to)
        self._record_index_map(start, end, dx * scale, dy * scale, to, "alpha")

        is_horizontal = to in "→←"
        width = 1.2 if is_horizontal else 0.8
        height = 0.8 if is_horizontal else 1.2

        x0 = self.x - (0 if is_horizontal else width / 2)
        y0 = self.y - (height / 2 if is_horizontal else 0)
        w = length * dx if is_horizontal else width
        h = height if is_horizontal else length * dy

        rect = patches.Rectangle(
            (x0, y0), w, h, linewidth=1, edgecolor="black", facecolor=color
        )
        self.ax.add_patch(rect)

        self._add_label_center(label, dx, dy, length)
        if add_label:
            self._add_label_start(start, dx, dy)
            self._add_label_end(end, dx, dy, length)
        self._move(to, length)

    def add_beta(
        self, start, end, label="", to="↑", color="green", scale=1.0, add_label=True
    ):
        length = (end - start) * scale
        dx, dy = self._direction_vector(to)
        self._record_index_map(start, end, dx * scale, dy * scale, to, "beta")

        arrow = patches.FancyArrow(
            self.x,
            self.y,
            dx * length,
            dy * length,
            width=0.8,
            length_includes_head=True,
            head_width=1.4,
            head_length=0.5,
            color=color,
        )
        self.ax.add_patch(arrow)

        self._add_label_center(label, dx, dy, length)
        if add_label:
            self._add_label_start(start, dx, dy)
            self._add_label_end(end, dx, dy, length)
        self._move(to, length)

    def add_star(self, start, end, shape="★", to="→", color="gold", size=24, scale=1.0):
        length = (end - start) * scale
        dx, dy = self._direction_vector(to)
        self._record_index_map(start, end, dx * scale, dy * scale, to, "star")

        tx = self.x + dx * length / 2
        ty = self.y + dy * length / 2
        self.ax.text(
            tx, ty, shape, fontsize=size, color=color, ha="center", va="center"
        )

        self._move(to, length)

    def show(self):
        self.ax.plot(*zip(*self.path), linestyle="dotted", color="gray")
        plt.tight_layout()
        plt.show()

    def save(self, path="protopo_plot.pdf", **kwargs):
        """
        Save the current topology plot to a file.

        Parameters:
            path (str): Output file path. Extension determines format (e.g. .pdf, .png).
        """
        self.ax.plot(*zip(*self.path), linestyle="dotted", color="gray")
        self.fig.savefig(path, **kwargs)
        print(f"Saved topology figure to {path}")
