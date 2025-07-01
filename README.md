# ProTopo

**ProTopo** is a lightweight Python package for drawing 2D protein topology diagrams — including α-helices, β-strands, linkers, N/C-terminals, and customizable markers.

It helps visualize and annotate simplified protein secondary structure layouts for publications, teaching, or design explanation.

[GitHub repo](https://github.com/hermanzhaozzzz/ProTopo) | [Issue page](https://github.com/hermanzhaozzzz/ProTopo/issues) | [Author](https://github.com/hermanzhaozzzz)

---

## ✨ Features

- Simple, declarative interface (`add_alpha`, `add_beta`, `add_linker`)
- Tracks residue indices and directions automatically
- Custom triangle/marker support (for cuts, domains, etc.)
- Outputs high-quality vector graphics via Matplotlib
- Built-in EGFP topology example

---

## 📦 Installation

```bash
pip install protopo
# or if developing locally
poetry install
```

## 🚀 Quick Example

```python
from protopo import ProTopo

pt = ProTopo()
pt.add_n_term()
pt.add_alpha(1, 10, label="H1", to="→")
pt.add_linker(10, 13, to="↓→", steps=(0.5, 0.5))
pt.add_beta(13, 22, label="S1", to="→")
pt.add_triangles([15])
pt.add_c_term()
pt.show()
```

![Qukck Example](resources/plot_demo_quickstart.png)

## 🧪 Run EGFP Example

This renders the schematic topology diagram of Enhanced GFP, following literature-reported secondary structure layout.

```python
from protopo.demo import draw_egfp
egfp = draw_egfp()
egfp.add_triangles(idx=[6, 211])
egfp.show()
# egfp.save("egfp.pdf")
```

then you will see

![EGFP demo](resources/plot_demo_egfp.png)

## 🧪 Testing

```bash
poetry run pytest
```

## License

MIT License © [Huanan Herman Zhao](https://github.com/hermanzhaozzzz)

## Contribution

Welcome to contribute codes and ideas
