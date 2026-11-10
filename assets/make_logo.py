"""Generate the course logo (assets/logo.svg).

Circular badge (circle = training material) in the MedicagoXeniumSpatialViewer palette: a tissue
section with segmented cells and transcripts, circled by the four weeks of
the course. Pure Python, fixed seed.

Usage: python3 assets/make_logo.py
"""
import math
import random
from pathlib import Path

OUT = Path(__file__).parent
BG, EDGE, TISSUE, TEDGE, SOFT, WHITE = "#1A1A2E", "#F2A541", "#23233D", "#3B3B60", "#A9ABC8", "#FFFFFF"
GENES = ["#4C9BE8", "#F28C3A", "#45B87A", "#B57BE0"]
HE_BG, HE_LINE, HE_NUC = "#F3C6DD", "#9C4F86", "#5B2A63"
FONT = 'font-family="Helvetica, Arial, sans-serif"'
S = 240  # badge size


def clip_halfplane(poly, a, b, c):
    out = []
    for i, p in enumerate(poly):
        q = poly[(i + 1) % len(poly)]
        fp, fq = a * p[0] + b * p[1] - c, a * q[0] + b * q[1] - c
        if fp <= 0:
            out.append(p)
        if fp * fq < 0:
            t = fp / (fp - fq)
            out.append((p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1])))
    return out


def voronoi(seeds, box):
    x0, y0, x1, y1 = box
    res = []
    for i, (px, py) in enumerate(seeds):
        poly = [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]
        for j, (qx, qy) in enumerate(seeds):
            if i != j:
                poly = clip_halfplane(poly, 2 * (qx - px), 2 * (qy - py), qx * qx + qy * qy - px * px - py * py)
        res.append(poly)
    return res


def in_poly(x, y, poly):
    inside = False
    for i, (ax, ay) in enumerate(poly):
        bx, by = poly[i - 1]
        if (ay > y) != (by > y) and x < (bx - ax) * (y - ay) / (by - ay) + ax:
            inside = not inside
    return inside


def pts(poly):
    return " ".join(f"{x:.1f},{y:.1f}" for x, y in poly)


class Badge:
    def __init__(self, seed):
        self.rng = random.Random(seed)
        self.b, self.defs = [], []

    def add(self, s):
        self.b.append(s)

    def g(self):
        r = self.rng
        return (r.random() + r.random() + r.random() - 1.5) / 1.5

    def frame(self):
        self.add(f'<rect x="6" y="6" width="{S - 12}" height="{S - 12}" rx="42" fill="{BG}" stroke="{EDGE}" stroke-width="7"/>')

    def name(self):
        self.add(f'<text x="{S / 2}" y="196" text-anchor="middle" {FONT} font-size="25" font-weight="bold" fill="{WHITE}">Relabel</text>')
        self.add(f'<text x="{S / 2}" y="218" text-anchor="middle" {FONT} font-size="13.5" letter-spacing="2.4" fill="{EDGE}">+ RESEGMENT</text>')

    def svg(self):
        return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {S} {S}" width="{S}" height="{S}" role="img" '
                f'aria-label="Spatial training course logo">\n<defs>{"".join(self.defs)}</defs>\n'
                + "\n".join(self.b) + "\n</svg>\n")

    def cells_in_circle(self, cx, cy, r, n, min_d):
        seeds = []
        while len(seeds) < n:
            a, d = self.rng.random() * 6.283, math.sqrt(self.rng.random()) * (r + 6)
            p = (cx + d * math.cos(a), cy + d * math.sin(a))
            if all(math.hypot(p[0] - q[0], p[1] - q[1]) > min_d for q in seeds):
                seeds.append(p)
        return seeds, voronoi(seeds, (cx - r - 10, cy - r - 10, cx + r + 10, cy + r + 10))

    def dom(self, x, y, cx, cy, r):
        d = int(((x - cx) / r * 0.8 + (y - cy) / r * 0.4 + 1.1) / 2.2 * 4 + self.g() * 0.4)
        return max(0, min(3, d))


def name(b):
    b.add(f'<text x="{S / 2}" y="178" text-anchor="middle" {FONT} font-size="24" font-weight="bold" fill="{WHITE}">Spatial</text>')
    b.add(f'<text x="{S / 2}" y="199" text-anchor="middle" {FONT} font-size="13" letter-spacing="2.6" fill="{EDGE}">TRAINING</text>')


def tissue(b, cx, cy, r, n_cells=14, dots=160, cid="tis"):
    b.defs.append(f'<clipPath id="{cid}"><circle cx="{cx}" cy="{cy}" r="{r}"/></clipPath>')
    b.add(f'<circle cx="{cx}" cy="{cy}" r="{r}" fill="{TISSUE}" stroke="{TEDGE}" stroke-width="1.5"/>')
    seeds, polys = b.cells_in_circle(cx, cy, r, n_cells, r * 0.34)
    g = [f'<g clip-path="url(#{cid})">']
    for (sx, sy), poly in zip(seeds, polys):
        g.append(f'<polygon points="{pts(poly)}" fill="none" stroke="#7C7CA8" stroke-width="1"/>')
    for _ in range(dots):
        a, d = b.rng.random() * 6.283, math.sqrt(b.rng.random()) * r
        x, y = cx + d * math.cos(a), cy + d * math.sin(a)
        g.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="1.4" fill="{GENES[b.dom(x, y, cx, cy, r)]}"/>')
    b.add("".join(g) + "</g>")


def logo():
    """Circle = training/course (hexagon = app/package, rounded square = pipeline)."""
    b = Badge(3)
    b.add(f'<circle cx="{S / 2}" cy="{S / 2}" r="{S / 2 - 6}" fill="{BG}" stroke="{EDGE}" stroke-width="7"/>')
    tissue(b, 120, 94, 34, n_cells=9, dots=100)
    pts4 = [(76, 52), (164, 52), (164, 136), (76, 136)]
    b.add(f'<path d="M{pts4[0][0]},{pts4[0][1]} L{pts4[1][0]},{pts4[1][1]} L{pts4[2][0]},{pts4[2][1]} L{pts4[3][0]},{pts4[3][1]}" '
          f'fill="none" stroke="{SOFT}" stroke-width="2" stroke-dasharray="5 4"/>')
    for i, (x, y) in enumerate(pts4, 1):
        b.add(f'<circle cx="{x}" cy="{y}" r="13" fill="{EDGE if i == 4 else BG}" stroke="{EDGE}" stroke-width="3"/>')
        b.add(f'<text x="{x}" y="{y + 5.5}" text-anchor="middle" {FONT} font-size="14.5" font-weight="bold" '
              f'fill="{BG if i == 4 else WHITE}">{i}</text>')
    name(b)
    return b


if __name__ == "__main__":
    out = OUT / "logo.svg"
    out.write_text(logo().svg())
    print("wrote", out)
