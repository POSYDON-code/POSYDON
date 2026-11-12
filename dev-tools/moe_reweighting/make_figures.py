"""Figures for the Moe & Di Stefano (2017) reweighting in norm_pop.

For each case, a simulation ("original") is reweighted to a target population
with ``norm_pop.calculate_model_weights``, and compared to a population drawn
directly from the target ("new sample"). Uncertainties come from resampling:
both populations are drawn ``N_REAL`` times with independent seeds; lines are
the median over the realisations and bands the 16-84 percentiles.

All distributions are per unit initial stellar mass (Msun^-1), so the
normalisation of the weights (mean-mass correction) is tested as well.

Run from the repository root:

    python dev-tools/moe_reweighting/make_figures.py

It writes ``figures/*.png`` and ``figures/summary.md`` next to this script.
"""

import os
import warnings

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from posydon.popsyn import norm_pop
from posydon.popsyn.independent_sample import generate_independent_samples
from posydon.popsyn.initial_distributions import (
    get_binary_fraction_distribution,
)

warnings.simplefilter("ignore")

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "figures")
N_SYS = 50_000      # systems per realisation
N_REAL = 30         # independent realisations

# categorical slots 1-3 of the dataviz reference palette (validated);
# each series also has its own line style
COLORS = {"original": "#2a78d6", "new": "#eb6834", "reweighted": "#1baf7a"}
STYLES = {"original": dict(ls="-", lw=1.5),
          "new": dict(ls="-", lw=2.0),
          "reweighted": dict(ls="--", lw=2.0)}
LABELS = {"original": "1. original sample",
          "new": "2. new sample (target)",
          "reweighted": "3. original, reweighted"}
INK = "#3d3d3a"
MUTED = "#8a8a85"
GRID = "#e6e5e0"

BASE = {
    "number_of_binaries": N_SYS,
    "primary_mass_scheme": "Kroupa2001",
    "primary_mass_min": 7.0,
    "primary_mass_max": 40.0,
    "binary_fraction_scheme": "const",
    "binary_fraction_const": 0.7,
    "secondary_mass_scheme": "flat_mass_ratio",
    "secondary_mass_min": 0.35,
    "secondary_mass_max": 140.0,
    "orbital_scheme": "period",
    "orbital_period_scheme": "Sana+12_period_extended",
    "orbital_period_min": 0.35,
    "orbital_period_max": 6e3,
    "eccentricity_scheme": "zero",
}


def moe(**overrides):
    kwargs = {k: v for k, v in BASE.items()
              if k not in ("binary_fraction_scheme", "binary_fraction_const",
                           "secondary_mass_min", "secondary_mass_max")}
    kwargs.update(secondary_mass_scheme="Moe+17-PsandQs",
                  orbital_period_scheme="Moe+17-PsandQs",
                  orbital_period_min=10**0.15, orbital_period_max=5e3)
    kwargs.update(overrides)
    return kwargs


def const(**overrides):
    kwargs = dict(BASE)
    kwargs.update(overrides)
    return kwargs


CASES = [
    dict(name="sana_to_moe",
         title="Sana+12, flat q, f_b = 0.7  →  Moe & Di Stefano (2017)",
         sim=const(),
         pop=moe()),
    dict(name="moe_to_sana",
         title="Moe & Di Stefano (2017)  →  Sana+12, flat q, f_b = 0.5",
         sim=moe(),
         # q >= 4 / M1 >= 0.1 and P inside the Moe range
         pop=const(binary_fraction_const=0.5, secondary_mass_min=4.0,
                   orbital_period_min=10**0.15, orbital_period_max=5e3)),
    dict(name="moe_to_moe",
         title="Moe & Di Stefano (2017): P = 1.4-5000 d  →  10-3000 d",
         sim=moe(),
         pop=moe(orbital_period_min=10.0, orbital_period_max=3e3)),
]

BINS = {
    "logM1": np.linspace(np.log10(7.0), np.log10(40.0), 21),
    "logP": np.linspace(-0.5, 3.8, 44),
    "q": np.linspace(0.0, 1.0, 26),
}


def draw(kwargs, seed):
    """Primary mass, q, log10 P (NaN for single stars) of one realisation."""
    kw = dict(kwargs, RNG=np.random.default_rng(seed))
    P, _, m1, m2 = generate_independent_samples(**kw)
    binary = np.isfinite(m2)
    data = {"S1_mass_i": m1,
            "S2_mass_i": np.where(binary, m2, 0.0),
            "orbital_period_i": P,
            "state_i": np.where(binary, "detached", "initially_single_star")}
    return data, binary


def hist(x, w, bins):
    """Density per unit x and per Msun formed."""
    h, _ = np.histogram(x, bins=bins, weights=w)
    return h / np.diff(bins)


def ratio(num, den):
    with np.errstate(invalid="ignore", divide="ignore"):
        return num / den


def statistics(data, binary, w):
    """Histograms and scalar summaries of one weighted realisation."""
    m1 = data["S1_mass_i"]
    lm = np.log10(m1)
    q = data["S2_mass_i"] / m1
    lp = np.log10(np.where(binary, data["orbital_period_i"], 1.0))
    wb = np.where(binary, w, 0.0)
    out = {
        "M1": hist(lm, w, BINS["logM1"]),
        "fbin": ratio(np.histogram(lm, BINS["logM1"], weights=wb)[0],
                      np.histogram(lm, BINS["logM1"], weights=w)[0]),
        "logP": hist(lp, wb, BINS["logP"]),
        "q": hist(q, wb, BINS["q"]),
        "q_vs_logP": ratio(np.histogram(lp, BINS["logP"], weights=wb * q)[0],
                           np.histogram(lp, BINS["logP"], weights=wb)[0]),
        "n_per_Msun": w.sum(),
        "f_bin": wb.sum() / w.sum(),
        "mean_logP": np.average(lp, weights=wb),
        "mean_q": np.average(q, weights=wb),
    }
    return out


def run_case(case):
    stats = {k: [] for k in COLORS}
    for r in range(N_REAL):
        sim, sim_bin = draw(case["sim"], seed=1000 + r)
        new, new_bin = draw(case["pop"], seed=5000 + r)
        M_sim = sim["S1_mass_i"].sum() + sim["S2_mass_i"].sum()
        M_new = new["S1_mass_i"].sum() + new["S2_mass_i"].sum()
        weights = norm_pop.calculate_model_weights(sim, M_sim, case["sim"],
                                                   case["pop"])
        stats["original"].append(statistics(
            sim, sim_bin, np.full(sim_bin.size, 1.0 / M_sim)))
        stats["new"].append(statistics(
            new, new_bin, np.full(new_bin.size, 1.0 / M_new)))
        stats["reweighted"].append(statistics(sim, sim_bin, weights))
    return {k: {q: np.array([s[q] for s in v]) for q in v[0]}
            for k, v in stats.items()}


def band(ax, x, arr, key, step=True):
    lo, med, hi = np.nanpercentile(arr, [16, 50, 84], axis=0)
    kw = dict(step="mid") if step else {}
    ax.fill_between(x, lo, hi, color=COLORS[key], alpha=0.22, lw=0, **kw)
    if step:
        ax.step(x, med, where="mid", color=COLORS[key], label=LABELS[key],
                **STYLES[key])
    else:
        ax.plot(x, med, color=COLORS[key], label=LABELS[key], **STYLES[key])


def style(ax, xlabel, ylabel):
    ax.set_xlabel(xlabel, color=INK)
    ax.set_ylabel(ylabel, color=INK)
    ax.tick_params(colors=MUTED, labelcolor=INK)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)
    ax.grid(color=GRID, lw=0.8)
    ax.set_axisbelow(True)


def plot_case(case, res):
    centre = {k: 0.5 * (b[1:] + b[:-1]) for k, b in BINS.items()}
    fig, axes = plt.subplots(2, 3, figsize=(15, 8.5))
    panels = [
        ("M1", "logM1", r"$\log_{10}(M_1/M_\odot)$",
         r"$dN/d\log_{10}M_1$ per $M_\odot$ formed", True),
        ("fbin", "logM1", r"$\log_{10}(M_1/M_\odot)$",
         "binary fraction", True),
        ("q", "q", r"$q = M_2/M_1$ (binaries)",
         r"$dN/dq$ per $M_\odot$ formed", True),
        ("logP", "logP", r"$\log_{10}(P/\mathrm{d})$ (binaries)",
         r"$dN/d\log_{10}P$ per $M_\odot$ formed", True),
        ("q_vs_logP", "logP", r"$\log_{10}(P/\mathrm{d})$ (binaries)",
         r"mean $q$ at given $P$", False),
    ]
    for ax, (key, bkey, xl, yl, step) in zip(axes.flat, panels):
        for series in COLORS:
            band(ax, centre[bkey], res[series][key], series, step=step)
        style(ax, xl, yl)
    # exact target binary fraction
    m1 = np.geomspace(7.0, 40.0, 400)
    exact, = axes[0, 1].plot(
        np.log10(m1),
        get_binary_fraction_distribution(case["pop"]).binary_fraction(m1),
        color=INK, lw=1.0, ls=":", label="target, exact")
    axes[0, 1].set_ylim(0, 1.05)
    axes[0, 1].legend(handles=[exact], frameon=False, fontsize=9,
                      loc="lower right")

    # relative residual of the reweighted sample to the new sample in log P
    ax = axes[1, 2]
    new_med = np.nanmedian(res["new"]["logP"], axis=0)
    for series in ("new", "reweighted"):
        rel = ratio(res[series]["logP"], new_med) - 1.0
        band(ax, centre["logP"], rel, series)
    ax.axhline(0, color=MUTED, lw=1)
    ax.set_ylim(-0.4, 0.4)
    style(ax, r"$\log_{10}(P/\mathrm{d})$ (binaries)",
          r"$dN/d\log_{10}P$ relative to new sample $-1$")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, frameon=False,
               fontsize=11, bbox_to_anchor=(0.5, 0.935))
    fig.suptitle(case["title"], color=INK, fontsize=14, y=0.995)
    fig.text(0.5, 0.945, f"median and 16–84% band over {N_REAL} independent "
             f"realisations of {N_SYS:,} systems each", ha="center",
             color=MUTED, fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.89))
    path = os.path.join(OUT, f"{case['name']}.png")
    fig.savefig(path, dpi=120, facecolor="white")
    plt.close(fig)
    return path


def fmt(arr, digits):
    lo, med, hi = np.percentile(arr, [16, 50, 84])
    return f"{med:.{digits}f} ± {0.5 * (hi - lo):.{digits}f}"


def summary_rows(case, res):
    rows = [f"**{case['title']}**", "",
            "| | systems per Msun (×10⁻²) | binary fraction | mean log₁₀ P "
            "(binaries) | mean q (binaries) |",
            "|---|---|---|---|---|"]
    for series in COLORS:
        s = res[series]
        rows.append(f"| {LABELS[series]} | {fmt(100 * s['n_per_Msun'], 3)} | "
                    f"{fmt(s['f_bin'], 3)} | {fmt(s['mean_logP'], 3)} | "
                    f"{fmt(s['mean_q'], 3)} |")
    return rows + [""]


def main():
    os.makedirs(OUT, exist_ok=True)
    lines = []
    for case in CASES:
        res = run_case(case)
        print("wrote", plot_case(case, res))
        lines += summary_rows(case, res)
    with open(os.path.join(OUT, "summary.md"), "w") as f:
        f.write("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
