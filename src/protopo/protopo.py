import math

import matplotlib.patches as patches
import matplotlib.pyplot as plt


class ProTopo:
    """A class for creating and visualizing protein topology diagrams.
    This class provides methods to add various protein structures like linkers,
    alpha helices, beta strands, and stars to a plot. It supports customization
    of colors, scales, and labels, and can save the resulting plot to a file.

    Example usage:
    >>> from protopo import ProTopo
    >>> pt = ProTopo()
    >>> pt.add_n_term()
    >>> pt.add_alpha(1, 10, label="H1", to="→")
    >>> pt.add_linker(10, 13, to="↓→", steps=(0.5, 0.5))
    >>> pt.add_beta(13, 22, label="S1", to="→")
    >>> pt.add_triangles([15])
    >>> pt.add_c_term()
    >>> pt.show()

    Methods:
    - add_n_term: Adds an N-terminal label to the current position.
    - add_c_term: Adds a C-terminal label to the current position.
    - add_linker: Adds a linker structure to the plot.
    - add_alpha: Adds an alpha helix structure to the plot.
    - add_beta: Adds a beta strand structure to the plot.
    - add_star: Adds a star structure to the plot.
    - add_triangles: Adds triangles to the plot based on specified indices.

    Attributes:
        fig (matplotlib.figure.Figure): The figure object for the plot.
        ax (matplotlib.axes.Axes): The axes object for the plot.
        x (float): Current x-coordinate in the plot.
        y (float): Current y-coordinate in the plot.
        path (list): List of tuples representing the path taken in the plot.
        index_map (dict): A mapping of indices to their positions, directions, and types.
    """

    def __init__(self, figsize: tuple = (12, 8)):
        """
        Initializes the ProTopo class with a specified figure size.
        Parameters:
            figsize (tuple): Size of the figure in inches (width, height).
        Example:
            >>> pt = ProTopo(figsize=(10, 6))
        Initializes a ProTopo instance with a figure size of 10x6 inches.
        """
        self.fig, self.ax = plt.subplots(figsize=figsize)
        self.x, self.y = 0, 0
        self.path = [(self.x, self.y)]
        self.index_map = {}  # 记录每个 index 的位置、方向和类型
        self.ax.set_aspect("equal")
        self.ax.axis("off")

    def _direction_vector(self, direction):
        """Converts a direction character to a vector (dx, dy).
        Parameters:
            direction (str): A single character representing the direction ('↑', '↓', '→', '←').
        Returns:
            tuple: A tuple (dx, dy) representing the direction vector.
        Example:
            >>> self._direction_vector("↑")
            (0, 1)
        """
        return {"↑": (0, 1), "↓": (0, -1), "→": (1, 0), "←": (-1, 0)}[direction]

    def _perpendicular_direction(self, direction):
        """Returns the perpendicular direction to the given direction.
        Parameters:
            direction (str): A single character representing the direction ('↑', '↓', '→', '←').
        Returns:
            str: The perpendicular direction character.
        Example:
            >>> self._perpendicular_direction("↑")
            "→"
        """
        return {"↑": "→", "↓": "→", "→": "↓", "←": "↓"}[direction]

    def _record_index_map(self, start, end, dx, dy, to, structure_type):
        """Records the position, direction, and type of a structure in the index map.
        Parameters:
            start (int): The starting index of the structure.
            end (int): The ending index of the structure.
            dx (float): The change in x-coordinate for the structure.
            dy (float): The change in y-coordinate for the structure.
            to (str): The direction of the structure ('↑', '↓', '→', '←').
            structure_type (str): The type of structure ('linker', 'alpha', 'beta').
        Example:
            >>> self._record_index_map(0, 3, 1.0, 0.0, "→", "alpha")
            >>> print(self.index_map)
            {
                0: {'x': 0.0, 'y': 0.0, 'to': '→', 'type': 'alpha'},
                1: {'x': 1.0, 'y': 0.0, 'to': '↑', 'type': 'linker'},
                2: {'x': 2.0, 'y': 0.0, 'to': '→', 'type': 'beta'},
            }
        """
        for i in range(start, end):
            self.index_map[i] = {
                "x": self.x + dx * (i - start),
                "y": self.y + dy * (i - start),
                "to": to,
                "type": structure_type if i != 1 else "N_term",
            }

    def _move(self, direction, length):
        """Moves the current position in the specified direction by a given length.
        Parameters:
            direction (str): A single character representing the direction ('↑', '↓', '→', '←').
            length (float): The distance to move in the specified direction.
        Example:
            >>> self._move("→", 5)
            >>> self.x, self.y
            (5.0, 0.0)
        """
        dx, dy = self._direction_vector(direction)
        self.x += dx * length
        self.y += dy * length
        self.path.append((self.x, self.y))

    def _draw_line(self, p1, p2, color="black"):
        """Draws a line between two points in the plot.
        Parameters:
            p1 (tuple): The starting point (x, y).
            p2 (tuple): The ending point (x, y).
            color (str): The color of the line.
        Example:
            >>> self._draw_line((0, 0), (1, 1), color="red")
            >>> self.ax.lines  # This will show the line added
        """
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
        """Adds triangles to the plot based on the specified indices and properties.
        Parameters:
            idx (list): A list of indices where triangles should be added.
            scale (float): The scale factor for the triangle size.
            offset (float): The offset for triangle positioning.
            facecolor (str): The fill color of the triangles.
            edgecolor (str): The edge color of the triangles.
            linewidth (float): The width of the triangle edges.
        Example:
            >>> pt.add_linker(0, 5, to="→")
            >>> pt.add_triangles([0, 1, 2])
            Adds triangles at indices 0, 1, and 2 of the linker.
        """
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
        """Adds a label at the start of a structure.
        Parameters:
            text (str): The label text to add.
            dx (float): The x-direction component for positioning.
            dy (float): The y-direction component for positioning.
        Example:
            >>> self._add_label_start("Start", 1, 0)
            Adds the label "Start" at the beginning of the current structure.
        """
        offset = 0.3
        x = self.x + dx * offset
        y = self.y + dy * offset
        self.ax.text(x, y, str(text), ha="center", va="center", fontsize=10)

    def _add_label_end(self, text, dx, dy, length):
        """Adds a label at the end of a structure.
        Parameters:
            text (str): The label text to add.
            dx (float): The x-direction component for positioning.
            dy (float): The y-direction component for positioning.
            length (float): The total length of the structure.
        Example:
            >>> self._add_label_end("End", 1, 0, 5)
            Adds the label "End" at the end of the current structure.
        """
        offset = 0.3
        x = self.x + dx * (length - offset)
        y = self.y + dy * (length - offset)
        self.ax.text(x, y, str(text), ha="center", va="center", fontsize=10)

    def _add_label_center(self, text, dx, dy, length):
        """Adds a label at the center of a structure.
        Parameters:
            text (str): The label text to add.
            dx (float): The x-direction component for positioning.
            dy (float): The y-direction component for positioning.
            length (float): The total length of the structure.
        Example:
            >>> self._add_label_center("Center", 1, 0, 5)
            Adds the label "Center" at the center of the current structure.
        """
        cx = self.x + dx * length / 2
        cy = self.y + dy * length / 2
        self.ax.text(cx, cy, text, ha="center", va="center", fontsize=12)

    def add_n_term(self):
        """Adds an N-terminal label to the current position in the plot.
        Example:
            >>> pt.add_n_term()
            Adds the label "N" at the current position.
        """
        assert self.x == 0
        assert self.y == 0
        self.ax.text(
            0,
            0,
            "N",
            ha="center",
            va="center",
            fontsize=12,
            fontweight="bold",
        )
        self._record_index_map(1, 2, 0, 0, "", "N_term")

    def add_c_term(self):
        """Adds a C-terminal label to the current position in the plot.
        Example:
            >>> pt.add_c_term()
            Adds the label "C" at the current position.
        """
        self.ax.text(
            self.x,
            self.y,
            "C",
            ha="center",
            va="center",
            fontsize=12,
            fontweight="bold",
        )
        idx = max(self.index_map.keys()) + 1
        self._record_index_map(
            idx, idx + 1, 0, 0, self.index_map[idx - 1]["to"], "C_term"
        )

    def add_linker(self, start, end, to="→", steps=(), color="black", scale=1.0):
        """Adds a linker structure to the plot.
        Parameters:
            start (int): The starting index of the linker.
            end (int): The ending index of the linker.
            to (str or list): The direction(s) of the linker ('↑', '↓', '→', '←').
            steps (list or tuple): The relative lengths of each segment of the linker.
            color (str): The color of the linker.
            scale (float): The scale factor for the linker length.
        Example:
            >>> pt.add_linker(4, 10, to="→", steps=[1.0])
            Adds a linker from index 4 to 10, moving right with a single segment.

            >>> pt.add_linker(4, 10, to="→↓", steps=[0.5, 0.5])
            Adds a two-segment linker from index 4 to 10, first moving right and
            then down, with each segment taking up half the total length.

            >>> pt.add_linker(4, 10, to="↑→↓", steps=[0.2, 0.6, 0.2])
            Adds a three-segment linker from index 4 to 10, first moving up and
            then right, and finally down, with each segment taking up a proportion
            of the total length.
        """
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
        for _to, step_ratio in zip(to, steps, strict=False):
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
        """Adds an alpha helix structure to the plot.
        Parameters:
            start (int): The starting index of the helix.
            end (int): The ending index of the helix.
            label (str): The label for the helix.
            to (str): The direction of the helix ('↑', '↓', '→', '←').
            color (str): The color of the helix.
            scale (float): The scale factor for the helix length.
            add_label (bool): Whether to add start and end labels.
        Example:
            >>> pt.add_alpha(1, 10, label="H1", to="→")
            Adds an alpha helix from index 1 to 10, moving right, with the label "H1".
        """
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
        """Adds a beta strand structure to the plot.
        Parameters:
            start (int): The starting index of the beta strand.
            end (int): The ending index of the beta strand.
            label (str): The label for the beta strand.
            to (str): The direction of the beta strand ('↑', '↓', '→', '←').
            color (str): The color of the beta strand.
            scale (float): The scale factor for the beta strand length.
            add_label (bool): Whether to add start and end labels.
        Example:
            >>> pt.add_beta(13, 22, label="S1", to="→")
            Adds a beta strand from index 13 to 22, moving right, with the label "S1".
        """
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
        """Adds a star structure to the plot.
        Parameters:
            start (int): The starting index of the star.
            end (int): The ending index of the star.
            shape (str): The star shape to use in a text style (default is "★").
            to (str): The direction of the star ('↑', '↓', '→', '←').
            color (str): The color of the star.
            size (int): The font size of the star.
            scale (float): The scale factor for the star length.
        Example:
            >>> pt.add_star(5, 10, shape="★", to="→", color="gold", size=24)
            Adds a star from index 5 to 10, moving right, with the shape "★",
            color "gold", and font size 24.
        """
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
        """Displays the current topology plot.
        Example:
            >>> pt.show()
            Displays the current plot with all added structures.
        """
        self.ax.plot(*zip(*self.path, strict=False), linestyle="dotted", color="gray")
        plt.tight_layout()
        plt.show()

    def close(self):
        """Closes the current figure and releases memory.

        Example:
            >>> pt.close()  # releases figure from memory
        """
        plt.close(self.fig)

    def save(self, path="protopo_plot.pdf", **kwargs):
        """
        Save the current topology plot to a file.

        Parameters:
            path (str): Output file path. Extension determines format (e.g. .pdf, .png).
            **kwargs: Additional keyword arguments passed to `matplotlib.pyplot.savefig`.
        """
        self.ax.plot(*zip(*self.path, strict=False), linestyle="dotted", color="gray")
        self.fig.savefig(path, **kwargs)
