
# Momentum-resolved excitation couplings in graphite {#sec:graphite}

Ultrafast electron diffuse scattering is an experimental technique pioneered over the past few years [@Chase2016; @Harb2016; @Waldecker2017]. The author also participated in the early development of the technique with work by M. Stern, L. P. René de Cotret *et al* [@Stern2018]. This early work on graphite was an experimental tour-de-force, but the observations remained qualitative rather than quantitative.

This chapter will detail how to extract the quantitative information encoded in UEDS experiments. 

TODO: important differentiator for UEDS vs X-ray is the amount of redundancy. The redundancy leads to the ability to invert the energy-integrated nature of UEDS.

## Single-crystal graphite

Single-crystal graphite is a two-dimensional material made of a repeating layers of carbon atoms arranged in a hexagonal lattice. Its lattice vectors $\left\{ \vec{a}_i \right\}$ are:
$$
\begin{pmatrix}
    \vec{a}_1 \\
    \vec{a}_2 \\
    \vec{a}_3
\end{pmatrix}
=
\begin{pmatrix}
    \sqrt{3}a  & -a & 0 \\
    0          & 2a & 0 \\
    0          &  0 & c 
\end{pmatrix}
\begin{pmatrix}
    \vec{e}_1 \\
    \vec{e}_2 \\
    \vec{e}_3
\end{pmatrix}
$$
where $a=\SI{1.232}{\angstrom}$, $c=\SI{6.711}{\angstrom}$, and $\left\{ \vec{e}_i \right\}$ are understood to be the usual Euclidean vectors. Carbon atoms $\left\{ \vec{c}_i \right\}$ are positioned as follows:
$$
\begin{pmatrix}
    \vec{c}_1 \\
    \vec{c}_2 \\
    \vec{c}_3 \\
    \vec{c}_4 
\end{pmatrix}
=
\begin{pmatrix}
    0            & 0            & 0            \\
    \sfrac{1}{3} & \sfrac{2}{3} & 0            \\
    0            & 0            & \sfrac{1}{2} \\
    \sfrac{2}{3} & \sfrac{1}{3} & \sfrac{1}{2} \\
\end{pmatrix}
\begin{pmatrix}
    \vec{a}_1 \\
    \vec{a}_2 \\
    \vec{a}_3
\end{pmatrix}
$$

Note carbon atoms $\vec{c}_1$ and $\vec{c}_2$ form a sheet of graphene, while $\vec{c}_3$ and $\vec{c}_4$ form a separate separate sheet of graphene rotated by $\tfrac{\pi}{3}$ (\ang{60}). The stacking of both sublattices along the $\vec{a}_3$ axis confers graphite increased discrete rotational symmetry (6-fold) about the stacking axis, compared to the 3-fold rotational symmetry of graphene. More specifically, the point group for graphite is $6/mmm$, while the Hermann-Mauguin symbol for the space group is $P6_3/mmc$. 

The real-space 6-fold rotational symmetry is mirrored in reciprocal space. The geometry of the Brillouin zone is shown in @fig:graphite-bz, including the high-symmetry points  $\vec{\Gamma}$ (zone-center), $\vec{M}$, and $\vec{K}$.

```{.python #fig:graphite-bz .matplotlib caption="In-plane section of the Brillouin zone of graphite, spanned by reciprocal lattice vectors $\vec{b}_1$ and $\vec{b}_2$. The three classes of high-symmetry points are shown. $\vec{M}$ and $\vec{K}$ have symmetry equivalents on every edge and vertex respectively."}
from math import radians, sqrt
from matplotlib.patches import Circle, RegularPolygon

def named_arrow(ax, x, y, dx, dy, text, tkwds=dict(), **kwargs):
    """ Draw an arrow annotated with some text """
    ax.arrow(x=x, y=y, dx=dx, dy=dy, **kwargs)
    ax.text(x=x + dx / 2, y=y + dy / 2, s=text, **tkwds)

fig, ax = plt.subplots(1, 1, figsize=(2, 2))
ax.set_xlim([-1.1, 1.1])
ax.set_ylim([-1.1, 1.1])
ax.set_aspect("equal")
ax.axis("off")

# reciprocal lattice vectors
arrow_kwds = dict(
    x=0, y=0, length_includes_head=True, width=0.001, head_width=0.07, fc="k"
)

named_arrow(
    ax,
    dx=-1 / 2,
    dy=-sqrt(3) / 2,
    text=r"$\mathbf{b}_1$",
    tkwds=dict(ha="right", va="bottom"),
    **arrow_kwds
)
named_arrow(
    ax,
    dx=1,
    dy=0,
    text=r"$\mathbf{b}_2$",
    tkwds=dict(va="bottom", ha="center"),
    **arrow_kwds
)

GAMMA = np.array((0, 0))
K = np.array((-1, 0))
M = np.array((0, sqrt(3) / 2))

R60deg = np.array([[1 / 2, -sqrt(3) / 2], [sqrt(3) / 2, 1 / 2]])

hexagon = RegularPolygon(
    xy=(0, 0),
    numVertices=6,
    orientation=radians(30),
    radius=1,
    facecolor="None",
    fill=False,
    edgecolor="k",
)
ax.add_patch(hexagon)
offset = np.array([-0.06, 0.06])
for point, label, sym_color in zip(
    [GAMMA, M, K],
    [r"$\mathbf{\Gamma}$", r"$\mathbf{M}$", r"$\mathbf{K}$"],
    ["k", "firebrick", "royalblue"],
):
    point_ = point
    for n in range(1, 7):
        point_ = R60deg @ point_
        ax.add_patch(
            Circle(xy=point_, radius=0.05, edgecolor="None", facecolor=sym_color)
        )

    ax.annotate(
        xy=point,
        ha="right",
        text=label,
        xytext=np.array(point) + offset,
        color=sym_color
    )
```

### Phonon landscape

With four atoms per unit cell, the structure of graphite supports twelve distinct phonon modes. The layered nature of graphite results in a clear distinction between in-plane and out-of-plane phonon modes. These modes are longitudinal modes LA, LO1 - LO3, in-plane transverse modes TA, TO1 - TO3, and out-of-plane traverse modes ZA, ZO1 - ZO3. The phonon dispersion relation for in-plane modes was calculated as described further below in @sec:computational-details, and is shown in @fig:graphite-static-dispersion. The out-of-plane transverse modes are not taken into account as the experiments are not sensitive to them. 

@fig:graphite-static-dispersion reveals why graphite is the *perfect* benchmark material for UEDS. All of the lattice modes are so stiff (high-energy) that, at room temperature, thermally-occupied modes are concentrated at $\vec{\Gamma}$. Therefore, the contrast between static diffuse intensity and diffuse intensity after photoexcitation is bound to be maximal.

```{.python #fig:graphite-static-dispersion .matplotlib caption="Phonon dispersion relation of graphite for in-plane modes LA, TA, and two-fold degenerate modes LO and TO. The path in reciprocal space is shown in the center. A horizontal dashed line at \SI{25}{\milli\electronvolt} indicates the average energy stored in the phonon modes at room temperature (\SI{300}{\kelvin})."}
from math import sqrt
from scipy.constants import physical_constants
import matplotlib.ticker as ticker
from pathlib import Path

HBAR = physical_constants["Planck constant over 2 pi in eV s"][0]
HZ_TO_EV = HBAR * 2 * np.pi

modes = ["LA", "LO2", "TA", "TO2"]

fig, ax = plt.subplots(1, 1, figsize=(FIGURE_WIDTH, 3))
ax_thz = ax.twinx()

for mode in modes:
    f = np.load(Path("data") / "graphite" / "static-dispersion" / f"{mode}.npy")
    ax.scatter(x=range(np.size(f)), y=f * HZ_TO_EV * 1000, s=2)  # meV
    ax_thz.scatter(x=range(np.size(f)), y=f * 1e-12, s=0)  # THz,

labels = [r"$\mathbf{\Gamma}$", r"$\mathbf{M}$", r"$\mathbf{K}$", r"$\mathbf{\Gamma}$"]
nsteps = np.size(f) / 3
for i in range(1, len(labels)):
    ax.axvline(x=i * nsteps, color="k", linestyle="--", linewidth=1)

ax.axhline(y=25, color="k", linestyle="--", linewidth=1)  # room temperature

# Draw BZ
# draw_hexagon(ax, center=(0,0), radius=100, color='k', transform=ax.transAxes)
hex_radius = 0.4
w, h = fig.get_size_inches()
path_ax = fig.add_axes([0.45, 0.10, 0.15, (w / h) * 0.15])
path_ax.axis("off")
path_ax.set_xlim([-0.5, 0.5])
path_ax.set_ylim([-0.5, 0.5])
draw_hexagon(path_ax, center=(0, 0), radius=hex_radius, color="k", facecolor='w')
for (x, y), label, ha in zip(
    [
        (0, 0),
        (0, hex_radius * sqrt(3) / 2),
        (hex_radius/2, hex_radius * sqrt(3) / 2),
    ],
    [r"$\mathbf{\Gamma}$", r"$\mathbf{M}$", r"$\mathbf{K}$"],
    ["right", "right", "left"],
):
    path_ax.scatter(x, y, s=5, c="k", zorder=np.inf)
    path_ax.annotate(label, (x, y), va="bottom", ha=ha)
    path_ax.arrow(x=0, y=0, dx=x, dy=y, linewidth=1, color='k', length_includes_head=True, )

# Annotate the modes
for label, pt in zip(
    ["TA", "LA", "TO", "LO"],
    [(0.14, 0.15), (0.07, 0.35), (0.1, 0.8), (0.25, 0.94)],
):
    ax.annotate(label, pt, color="k", textcoords=ax.transAxes)

ax.xaxis.set_minor_locator(ticker.NullLocator())
ax.xaxis.set_major_locator(
    ticker.FixedLocator([i * nsteps for i in range(len(labels))])
)

ax.xaxis.set_major_formatter(plt.FixedFormatter(labels))
ax.set_xlim([0, np.size(f)])
ax.set_ylabel("Energy [meV]")
ax.set_ylim([0, 210])

ax_thz.set_ylabel("Frequency [THz]")
ax_thz.yaxis.set_label_position("right")
```

## The first hundred femtoseconds viewed by trARPES

Electron scattering measurements are not able to isolate to the dynamics of the electronic system. Although UED and UEDS are sensitive to the dynamics of electrons whizzing about in a material, understanding the electronic contribution to measured signals requires some prior knowledge. To understand the ultrafast electron diffuse scattering results presented further, it is important to understand the effect of photoexcitation during first \SI{100}{\femto\second} as observed by time- and angle-resolved photoemission spectroscopy (trARPES).

## Experimental and computational methods

### Sample preparation

Single-crystal flakes of natural graphite \SIrange{10}{90}{\nano\meter} thick were prepared using a mechanical exfoliation procedure analogous to the work by Novoselov *et al* [@Novoselov2004], briefly described here. Thick flakes were embedded in Crystalbond glue on a \SI{3}{\milli\meter} copper TEM grid (200 lines per inch). The embedded flakes are then exfoliated using ordinary adhesive tape. The procedure was repeated until the embedded flakes were translucent when observed under an optical microscope. The glue is then delicately washed away with a solvent. The choice of the solvent is dependent on the glue used; in the present case, acetone was used. Sample thickness has been measured directly using atomic force microscopy characterization (see @fig:graphite-afm). Once a suitable sample has been identified, an aperture was made using aluminum foil to isolate a sample region with uniform thickness. The resulting sample used in this work covered 500 × 500 \si{\square\micro\meter}, with a thickness of \SI{70}{\nano\meter}.

![Characterization of the thickness of potential graphite samples. **Left** optical microscope image of exfoliated graphite. Translucent regions A and B are highlighted as potential samples. **Right** Atomic force microscope measurement shows a sample thickness of \SI{88}{\nano\meter}. Modified from Chatelain [@Chatelain2014].](diagrams/graphite_exfoliation.pdf){#fig:graphite-afm}

### Data acquisition 

The UEDS experiments presented in this chapter made use of the experimental setup presented in @sec:experimental_setup. Ultrashort laser pulses of light  were shone at $t=t_0$ on a thin single-crystal specimen of graphite, oriented in the $\left[001\right]$ direction. Compressed electron bunches with 10^7^ electrons per bunch were transmitted through the sample at $t=t_0 + \tau$. The time-delay $\tau$ was scanned from \SI{-40}{\pico\second} to \SI{680}{\pico\second}.

The interrogated film were pumped with a pump spot of 1 × 1 \si{\square\milli\meter} full-width at half-maximum (FWHM), ensuring nearly uniform illumination of the probed volume. The film was pumped at a fluence of \SI{12}{\milli\joule\per\square\centi\meter}, resulting in an absorbed energy density of \SI{8}{\joule\per\cubic\meter}. The scattering patterns are collected with a Gatan Ultrascan 895 camera: a 2.54 × 2.54 \si{\square\cm} phosphor screen fiber coupled to a 2048 px × 2048 px charge-coupled detector (CCD) placed \SI{25}{\centi\meter} away from the sample. The experiment herein consists of time delays in the range of $\SIrange{-40}{680}{\pico\second}$. Per-pixel scattering intensity fluctuations over laboratory time reveals a transient dynamic range of 1 : 10^8^, allowing the acquisition of diffraction patterns and diffuse scattering patterns simultaneously. A static diffraction pattern is shown in @fig:graphite-static a).

Due to the flatness of the Ewald sphere for \SI{90}{\kilo\electronvolt} electrons, many symmetry-related reflections are visible within each pattern. The information contained in a set of symmetry-equivalent reflections is redundant due to the point-group symmetry of the scattering crystal. As long as the point-group symmetry is not broken by photoexcitation its consequences, it is possible to harness the redundancy to enhance the signal-to-noise ratio of a UEDS dataset. In the case of graphite, no observable symmetry-breaking phenomena is brought on by photoexcitation at \SI{1.55}{\electronvolt} when looking at the raw data. Moreover, trARPES experiments [@Stange2015] do not show the opening of a gap in the electronic band structure, albeit at much lower photoexcitation densities, which would be indicative of point-group symmetry breaking. The point-group of graphite is $6/mmm$, which encompasses 6-fold discrete rotational symmetry in the $\vec{a} \times \vec{b}$ plane. Therefore, specifically in the case of graphite oriented in the $\left[ 001 \right]$ direction, we can safely enhance the diffuse signals by a factor of $\sqrt{6}$ by the use of a six-fold discrete azimuthal average:
$$
    I(\vec{q}, \tau) \to \frac{1}{6} \sum_{n=1}^6 I( \vec{R}(\tfrac{\pi n}{3}) \cdot \vec{q}, \tau)
$$
where $\vec{R}(\theta)$ is the in-plane rotation matrix:
$$
    \vec{R}(\theta) = \begin{pmatrix}
                      \cos{\theta} & -\sin{\theta} & 0\\
                      \sin{\theta} & \cos{\theta}  & 0\\
                      0            & 0             & 1
                     \end{pmatrix}
$$
Throughout the rest of this chapter, "scattering intensity" will imply discrete azimuthal average unless otherwise noted. An example of six-fold averaged diffraction pattern is shown in @fig:graphite-static b).

```{.python #fig:graphite-static .matplotlib caption="Static diffraction pattern of graphite. Brillouin zones are shown around each reflection to guide the eye. **a)** static, unprocessed diffraction pattern. **b)** Six-fold discrete azimuthal average of the diffraction pattern in a) results in $\sqrt{6}$ increase in signal-to-noise ratio."}
from matplotlib.ticker import FixedLocator, FixedFormatter
from pathlib import Path
from iris import DiffractionDataset
from skimage.transform import rotate
from skued import nfold, detector_scattvectors, combine_masks
from math import floor
import itertools as it

xc, yc = GRAPHITE_CENTER

xx, yy = np.meshgrid(np.arange(0, 2048), np.arange(0, 2048))
rr = np.sqrt(np.square(xx - xc) + np.square(yy - yc))

beamblock = np.ones((2048, 2048), dtype=np.bool)
beamblock[0:1260, 900:1130] = False

artifact_mask = np.ones((2048, 2048), dtype=np.bool)
artifact_mask[1084::, 437:482] = False
artifact_mask[0:932, 1296:1324] = False

mask = combine_masks(beamblock, artifact_mask)

with DiffractionDataset(Path("data") / "graphite" / "graphite_time_corrected_iris5.hdf5") as source:
    b4t0 = source.diff_eq()

b4t0_symmetrized = np.array(b4t0, copy=True)
b4t0_symmetrized = nfold(b4t0, mod=6, center=GRAPHITE_CENTER, mask=mask)
b4t0_symmetrized[rr < 125] = 0

b4t0[:] = rotate(b4t0, angle=GRAPHITE_ANGLE, center=GRAPHITE_CENTER, mode="reflect")
b4t0_symmetrized[:] = rotate(b4t0_symmetrized, angle=GRAPHITE_ANGLE, center=GRAPHITE_CENTER, mode="reflect")

qx, qy, _ = detector_scattvectors(
    keV=90, 
    camera_length=GRAPHITE_CAMERA_LENGTH,
    shape=(2048, 2048),
    pixel_size=14e-6, 
    center=(yc, xc),
)

# Determine the smallest center -> side distance, and crop around that
side_length = floor(min([xc, abs(xc - 2048), yc, abs(yc - 2048)]))
xs, ys = (
    slice(yc - side_length, yc + side_length),
    slice(xc - side_length, xc + side_length),
)
b4t0 = b4t0[xs, ys]
b4t0_symmetrized = b4t0_symmetrized[xs, ys]
qx = qx[ys, xs]
qy = qy[ys, xs]

fig = plt.figure(figsize = (FIGURE_WIDTH, FIGURE_WIDTH/2))
grid = ImageGrid(fig, 111, nrows_ncols=(1,2), cbar_location="top")

for ax, im, label in zip(grid, [b4t0, b4t0_symmetrized], ["a)", "b)"]):
    m = ax.imshow(
        im,
        vmin=0,
        vmax=200,
        cmap="magma",
        extent=[qx.min(), qx.max(), qy.min(), qy.max()],
    )
    draw_hexagon_field(
        ax,
        radius=1.7,
        crystal=Crystal.from_pwscf(Path("data") / "graphite" / "output.out"),
        reflections=it.product(
            range(-4,4), range(-4,4), [0]
        ),
        color="w",
        linewidth=0.5,
        linestyle=":",
    )
    tag_axis(ax, text=label)

for ax in grid:
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)

cbar = grid[0].cax.colorbar(
    m,
    ticks=FixedLocator(locs=[0, 200]),
    format=FixedFormatter(["0", "1"]),
)
cbar.ax.xaxis.set_label_position("top")
cbar.ax.set_xlabel("Scattering intensity [a.u.]")
```

### Computational details {#sec:computational-details}

## Diffuse intensity dynamics

The change in scattering intensity for representative time-delays is shown in

```{.python #fig:graphite-ueds .matplotlib caption="Change in scattering intensity $\Delta I(\vec{q}, t=\tau) = I(\vec{q}, \tau) - I(\vec{q}, -\infty)$ of photoexcited graphite for a few representative time-delays $\tau$. Hexagonal Brillouin zones are shown on half of the reflections to guide the eye. Scattering patterns show diffuse dynamics in the range of $|\vec{q}| < \SI{12}{\per\angstrom}$. Negative going features (blue) are exclusively due to the transient Debye-Waller effect on the Bragg peaks. All positive changes (red) are dynamics of the diffuse scattering intensity."}
from iris import DiffractionDataset
from math import floor
from pathlib import Path
from matplotlib.ticker import FixedLocator, FixedFormatter
from skimage.transform import rotate
from skued import detector_scattvectors
import itertools as it

xc, yc = GRAPHITE_CENTER
qx, qy, _ = detector_scattvectors(
    keV=90, 
    camera_length=GRAPHITE_CAMERA_LENGTH,
    shape=(2048, 2048),
    pixel_size=14e-6, 
    center=(yc, xc),
)

# Determine the smallest center -> side distance, and crop around that
side_length = floor(min([xc, abs(xc - 2048), yc, abs(yc - 2048)]))
xs = slice(yc - side_length, yc + side_length)
ys = slice(xc - side_length, xc + side_length)

qx = qx[ys, xs]
qy = qy[ys, xs]
qq = np.sqrt(qx ** 2 + qy ** 2)

fig = plt.figure(figsize = (FIGURE_WIDTH, FIGURE_WIDTH))
grid = ImageGrid(fig, 111, nrows_ncols=(2,2), cbar_location="top")

with DiffractionDataset(Path("data") / "graphite" / "TDS_delta.hdf5") as dset:
    for time, ax, letter in zip([0.5, 1.5, 5, 100], grid, "abcd"):

        image = dset.diff_data(time)
        image[:] = rotate(image, angle=GRAPHITE_ANGLE, center=GRAPHITE_CENTER, mode="reflect")
        image = image[xs, ys]
        image[qq < 1.5] = 0
        m = ax.imshow(
            image,
            cmap="seismic",
            vmin=-0.4,
            vmax=0.4,
            extent=[qx.min(), qx.max(), qy.min(), qy.max()],
        )

        draw_hexagon_field(
            ax,
            radius=1.7,
            crystal=Crystal.from_pwscf(Path("data") / "graphite" / "output.out"),
            reflections=it.product(
                [-4, -3, -2, -1, 0], [-4, -3, -2, -1, 0, 1, 2, 3], [0]
            ),
            color="k",
            linestyle=":",
        )

        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)

        tag = f"{letter}) {1e3 * time:0.0f} fs" if time < 1 else f"{letter}) {time} ps"
        tag_axis(ax, text=tag)

cbar = ax.cax.colorbar(
    m,
    ticks=FixedLocator(locs=[-0.4, 0, 0.4]),
    format=FixedFormatter(["-1", "0", "1"]),
)
cbar.ax.xaxis.set_label_position("top")
cbar.ax.set_xlabel("Scattering intensity change [a.u.]")
```

## The one-phonon structure factor

```{.python .matplotlib caption="Calculated one-phonon structure factors $|F_{1j}(\vec{q}, t=0)|^2$ of in-plane transverse modes at \SI{300}{\kelvin}, for scattering vectors $\vec{q}$ equivalents to the detector area shown in @fig:graphite-ueds."}
from pathlib import Path
from matplotlib.ticker import FixedFormatter, FixedLocator
from crystals import Crystal

INPUT = Path("data") / "graphite"

# Mode ordering of graphite according to the file
# Gra-C_XDM_mode_grid_new2.json
MODE_ORDERING = {
    "LA": 0,
    "TA": 1,
    "ZA": 2,
    "LO1": 3,
    "LO2": 4,
    "LO3": 5,
    "TO1": 6,
    "TO2": 7,
    "TO3": 8,
    "ZO1": 9,
    "ZO2": 10,
    "ZO3": 11,
}
MODES = sorted(MODE_ORDERING.keys())
IN_PLANE_MODES = sorted(set(MODE_ORDERING.keys()) - {"ZA", "ZO1", "ZO2", "ZO3"})

in_plane_refls = filter(
    lambda tup: tup[2] == 0, Crystal.from_database("C").bounded_reflections(12)
)

reflections = tuple(in_plane_refls)

fig = plt.figure(figsize=(FIGURE_WIDTH, 1.1 * FIGURE_WIDTH))
grid = ImageGrid(
    fig,
    111,
    nrows_ncols=(2, 2),
    cbar_location="top",
)

qx = np.load(INPUT / "oneph" / "qx.npy")
qy = np.load(INPUT / "oneph" / "qy.npy")
bragg_peaks = np.load(INPUT / "oneph" / "bragg_peaks.npy")
cryst = Crystal.from_pwscf(INPUT / "output.out")

# Only Longitudinal modes here
modes = filter(lambda s: s.startswith("L"), IN_PLANE_MODES)
for mode, ax in zip(modes, grid):
    image = np.load(INPUT / "oneph" / f"{mode}.npy")

    # Image is scaled so maximum is always 1
    m = ax.imshow(
        image / image.max(),
        extent=[qx.min(), qx.max(), qy.min(), qy.max()],
        cmap="inferno",
        vmin=0,
        vmax=1,
    )
    ax.scatter(x=bragg_peaks[:, 0], y=bragg_peaks[:, 1], s=1, c="w")

    # Bragg peaks might extend beyond the rest of the image
    ax.set_xlim([qx.min(), qx.max()])
    ax.set_ylim([qy.min(), qy.max()])
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)

    tag_axis(ax, text=f"{mode}")

    draw_hexagon_field(
        ax=ax,
        radius=1.7,
        crystal=cryst,
        color=(0.7, 0.7, 0.7, 1),  # light gray
        linestyle=":",
        reflections=reflections,
    )

cbar = ax.cax.colorbar(
    m,
    ticks=FixedLocator([0, m.get_array().max()]),
    format=FixedFormatter(["0", "1"]),
)
cbar.ax.xaxis.set_label_position("top")
cbar.ax.set_xlabel(r"$|F_{1j}(\mathbf{q}, t_0)|^2$ [a.u.]")

set_height_auto(fig)
```

```{.python .matplotlib caption="Calculated one-phonon structure factors $|F_{1j}(\vec{q}, t=0)|^2$ of in-plane transverse modes at \SI{300}{\kelvin}, for scattering vectors $\vec{q}$ equivalents to the detector area shown in @fig:graphite-ueds."}
from pathlib import Path
from matplotlib.ticker import FixedFormatter, FixedLocator
from crystals import Crystal

INPUT = Path("data") / "graphite"

# Mode ordering of graphite according to the file
# Gra-C_XDM_mode_grid_new2.json
MODE_ORDERING = {
    "LA": 0,
    "TA": 1,
    "ZA": 2,
    "LO1": 3,
    "LO2": 4,
    "LO3": 5,
    "TO1": 6,
    "TO2": 7,
    "TO3": 8,
    "ZO1": 9,
    "ZO2": 10,
    "ZO3": 11,
}
MODES = sorted(MODE_ORDERING.keys())
IN_PLANE_MODES = sorted(set(MODE_ORDERING.keys()) - {"ZA", "ZO1", "ZO2", "ZO3"})

in_plane_refls = filter(
    lambda tup: tup[2] == 0, Crystal.from_database("C").bounded_reflections(12)
)

reflections = tuple(in_plane_refls)

fig = plt.figure(figsize=(FIGURE_WIDTH, 1.1 * FIGURE_WIDTH))
grid = ImageGrid(
    fig,
    111,
    nrows_ncols=(2, 2),
    cbar_location="top",
)

qx = np.load(INPUT / "oneph" / "qx.npy")
qy = np.load(INPUT / "oneph" / "qy.npy")
bragg_peaks = np.load(INPUT / "oneph" / "bragg_peaks.npy")
cryst = Crystal.from_pwscf(INPUT / "output.out")

# Only transverse modes here
modes = filter(lambda s: s.startswith("T"), IN_PLANE_MODES)
for mode, ax in zip(modes, grid):
    image = np.load(INPUT / "oneph" / f"{mode}.npy")

    # Image is scaled so maximum is always 1
    m = ax.imshow(
        image / image.max(),
        extent=[qx.min(), qx.max(), qy.min(), qy.max()],
        cmap="inferno",
        vmin=0,
        vmax=1,
    )
    ax.scatter(x=bragg_peaks[:, 0], y=bragg_peaks[:, 1], s=1, c="w")

    # Bragg peaks might extend beyond the rest of the image
    ax.set_xlim([qx.min(), qx.max()])
    ax.set_ylim([qy.min(), qy.max()])
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)

    tag_axis(ax, text=f"{mode}")

    draw_hexagon_field(
        ax=ax,
        radius=1.7,
        crystal=cryst,
        color=(0.7, 0.7, 0.7, 1),  # light gray
        linestyle=":",
        reflections=reflections,
    )

cbar = ax.cax.colorbar(
    m,
    ticks=FixedLocator([0, m.get_array().max()]),
    format=FixedFormatter(["0", "1"]),
)
cbar.ax.xaxis.set_label_position("top")
cbar.ax.set_xlabel(r"$|F_{1j}(\mathbf{q}, t_0)|^2$ [a.u.]")

set_height_auto(fig)
```

\FloatBarrier
## References {.unnumbered}
\printbibliography[heading=none]
