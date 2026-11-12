# Moe & Di Stefano (2017): from the paper's equations to the code

This document puts each equation of Moe & Di Stefano (2017, ApJS 230, 15; arXiv:1606.05347)
next to the code that implements it in
[`posydon/popsyn/Moes_distributions.py`](../posydon/popsyn/Moes_distributions.py).
The paper's fits are module-level functions, and the class `Moe_17_PsandQs` combines
them into a sampler and a probability density function (PDF).

Sect. 8 lists how this differs from the previous, grid-based implementation (POSYDON
≤ v2.4), and which of its problems are fixed.

Equation numbers are the paper's. $q = M_2/M_1$; $\log$ is base 10, with $P$ in days and
$M_1$ in $M_\odot$ (units are left out of the equations). In the code, `m` is $M_1$ and
`lp` is $\log P$.

---

## 0. Structure of the model

For a primary of mass $M_1$, the sampler draws one of two outcomes:

- **a single star**, with probability $1-f_\mathrm{bin}(M_1)$;
- **one inner binary**, with joint density

$$
p(\log P, q, e \mid M_1) \;=\;
\underbrace{f_{\log P;q>0.3}(M_1,\log P)\; p_\mathrm{inner}(M_1,\log P)}_{\text{period, Sect. 4}}
\;\times\;
\underbrace{p_q(q \mid M_1,\log P)}_{\text{mass ratio, Sect. 2}}
\;\times\;
\underbrace{p_e(e \mid M_1,\log P)}_{\text{eccentricity, Sect. 3}} .
$$

Here $p_q$ is normalised to 1 over $0.3<q<1$, as in the paper, so its integral over $0.1<q<1$
is the factor $\mathcal{Q}$ of Sect. 2.3. The density integrates to $f_\mathrm{bin}(M_1)$.

### What is introduced where

The paper fits five parameters as functions of $(M_1, P)$. Three describe the shape of
the mass-ratio distribution, one the eccentricity distribution, and one how many
companions there are per decade of period. The document first gives each parameter's
fit, then shows how the parameters are combined into a distribution, and finally how
the code samples from it.

| Symbol | Meaning | Paper | Code | Defined in | Used in |
|---|---|---|---|---|---|
| $\mathcal{F}_\mathrm{twin}$ | excess fraction of near-equal-mass twins, $0.95<q<1$ | Eqns. 5–7 | `f_twin` | Sect. 1.1 | $p_q$, Sect. 2 |
| $\gamma_\mathrm{largeq}$ | power-law slope of $p_q$ over $0.3<q<1$ | Eqns. 9–11 | `gamma_largeq` | Sect. 1.2 | $p_q$, Sect. 2 |
| $\gamma_\mathrm{smallq}$ | power-law slope of $p_q$ over $0.1<q<0.3$ | Eqns. 13–15 | `gamma_smallq` | Sect. 1.3 | $p_q$, Sect. 2 |
| $p_q$ | mass-ratio distribution built from the three above | Eqn. 2, Fig. 2 | `q_components`, `q_pdf` | Sect. 2 | sampling q (Sect. 5), PDF (Sect. 6) |
| $\mathcal{Q}$ | $\int_{0.1}^{1}p_q\,dq$: converts $q>0.3$ frequencies to $q>0.1$ | Sect. 9.3 | `q_factor` | Sect. 2.3 | period density, Sect. 4.3 |
| $e_\mathrm{max}$ | Roche-lobe limit on $e$ | Eqn. 3 | `e_max`, `is_circular` | Sect. 3.1 | $p_e$, Sect. 3.3 |
| $\eta$ | power-law slope of $p_e$ | Eqns. 17–18 | `eta` | Sect. 3.2 | $p_e$, Sect. 3.3 |
| $p_e$ | eccentricity distribution | Sect. 2 | `e_pdf` | Sect. 3.3 | sampling e (Sect. 5), PDF (Sect. 6) |
| $f_{\log P;q>0.3}$ | companions with $q>0.3$ per decade of $P$ | Eqns. 20–23 | `f_logP_q03` | Sect. 4.1 | period density, Sect. 4.3 |
| $p_\mathrm{inner}$ | probability that a companion is the inner binary | Moe's code (not in paper) | `prob_inner_binary` | Sect. 4.2 | period density, Sect. 4.3 |
| $g$, $f_\mathrm{bin}$ | inner-binary period density, and its integral (the binary fraction) | — | `inner_binary_density`, `binary_fraction` | Sect. 4.3 | sampling P (Sect. 5), PDF (Sect. 6) |

The rest of the document is ordered as follows:
- **Sects. 1–4:** the fits, each shown next to its code.
- **Sect. 5:** how the sampler draws from them.
- **Sect. 6:** the PDF.
- **Sect. 7:** a summary of where the code differs from the paper.
- **Sect. 8:** the changes from the previous implementation.

### Mass interpolation (used by several equations)

The paper gives fits at a few anchor masses and says to "interpolate according to
$M_1$" in between. The code interpolates **linearly in $\log M_1$**, as Moe's own code
does, and holds the value constant outside the anchors (`_wlog`, `_mass_interp`).
- The $\gamma$'s are interpolated between 1.2, 3.5 and 6 $M_\odot$.
- $\eta$ is interpolated between 3 and 7 $M_\odot$.

Primaries outside 0.8–40 $M_\odot$ use the statistics of 0.8 or 40 $M_\odot$
(`Moe_17_PsandQs._m_stats`).

---

## 1. Mass-ratio parameters

### 1.0 What the three parameters describe

For a given $M_1$ and $P$, the paper models the mass-ratio distribution over
$0.1<q<1$ with three numbers (paper Sect. 2 and Fig. 2):

- **$\gamma_\mathrm{largeq}$**, the power-law slope $p_q\propto q^{\gamma}$ over $0.3<q<1$.
  $\gamma=0$ is a flat distribution; $\gamma=-2.35$ is what random pairing from a
  Salpeter IMF would give. Most surveys detect companions down to $q\approx0.3$, so
  this slope is the best constrained of the three.
- **$\gamma_\mathrm{smallq}$**, the power-law slope over $0.1<q<0.3$. The distribution
  is required to be **continuous at $q=0.3$**, so this slope only sets how quickly
  the density changes below 0.3, not its level there. It mainly controls how many
  small-q companions there are.
- **$\mathcal{F}_\mathrm{twin}$**, an *excess* fraction of twins with $q>0.95$, measured
  relative to the $q>0.3$ power law. It is a narrow, flat bump on top of the power
  law, not the total fraction of systems above 0.95.

A sketch of a typical shape (negative slopes; the paper's Fig. 2 shows an example
with $\gamma_\mathrm{smallq}=0.5$, $\gamma_\mathrm{largeq}=-1$, $\mathcal{F}_\mathrm{twin}=0.1$):

```
 p_q
  │ *
  │   *                                         ┌───┐  ← twin excess F_twin,
  │     *   q^γ_smallq                          │   │    flat over 0.95–1,
  │        *                                    │   │    added on top
  │           *  *   q^γ_largeq                 │   │
  │                *  *   *   *   *   *   *   * │   │
  │                                             │   │
  └─────────────────┬───────────────────────────┬───┬──── q
 0.1               0.3                        0.95  1
   ←  small q  →   ↑ continuous here   ←  large q  →
```

Each parameter depends on $(M_1,P)$. The paper fits them at a few anchor masses as
piecewise-linear functions of $\log P$, and interpolates in $\log M_1$ in between
(Sect. 0). Overall, binaries become weighted towards smaller $q$ at longer periods:
both slopes fall with $\log P$, and twins appear only at short periods for massive
primaries.

This section gives the three fits (Sects. 1.1–1.3). Sect. 2 assembles them into
$p_q$, and explains the normalisation convention: the $q>0.3$ part integrates to 1,
and the $0.1<q<0.3$ part adds companions on top of that.

### 1.1 Excess twin fraction, Eqns. 5–7

$$
\mathcal{F}_\mathrm{twin}(M_1,\log P)=
\begin{cases}
\mathcal{F}_{\mathrm{twin};\log P<1} & \log P<1\\[2pt]
\mathcal{F}_{\mathrm{twin};\log P<1}\left(1-\dfrac{\log P-1}{\log P_\mathrm{twin}-1}\right) & 1\le\log P<\log P_\mathrm{twin}\\[6pt]
0 & \log P\ge\log P_\mathrm{twin}
\end{cases}
\tag{5}
$$

$$
\mathcal{F}_{\mathrm{twin};\log P<1}(M_1)=0.30-0.15\log M_1
\tag{6}
$$

$$
\log P_\mathrm{twin}(M_1)=
\begin{cases}8.0-M_1/M_\odot & M_1\le 6.5\,M_\odot\\ 1.5 & M_1>6.5\,M_\odot\end{cases}
\tag{7}
$$

`f_twin`:
```python
ft1 = 0.30 - 0.15 * np.log10(m)                       # Eqn. 6
lptwin = np.where(m <= 6.5, 8.0 - m, 1.5)             # Eqn. 7
return np.where(lp < 1.0, ft1,                         # Eqn. 5
                np.where(lp < lptwin,
                         ft1 * (1.0 - (lp - 1.0) / (lptwin - 1.0)), 0.0))
```

### 1.2 Large-q slope $\gamma_\mathrm{largeq}$, Eqns. 9–11

$$
\gamma_\mathrm{largeq}(M_1<1.2)=\begin{cases}-0.5 & 0.2\le\log P<5\\ -0.5-0.3(\log P-5) & 5\le\log P<8\end{cases}
\tag{9}
$$

$$
\gamma_\mathrm{largeq}(M_1=3.5)=\begin{cases}-0.5 & \log P<1\\ -0.5-0.2(\log P-1) & 1\le\log P<4.5\\ -1.2-0.4(\log P-4.5) & 4.5\le\log P<6.5\\ -2.0 & 6.5\le\log P<8\end{cases}
\tag{10}
$$

$$
\gamma_\mathrm{largeq}(M_1>6)=\begin{cases}-0.5 & \log P<1\\ -0.5-0.9(\log P-1) & 1\le\log P<2\\ -1.4-0.3(\log P-2) & 2\le\log P<4\\ -2.0 & 4\le\log P<8\end{cases}
\tag{11}
$$

`gamma_largeq`: the three equations are `g12`, `g35` and `g6`, combined with
`_mass_interp(m, g12, g35, g6, 1.2, 3.5, 6.0)`.

### 1.3 Small-q slope $\gamma_\mathrm{smallq}$, Eqns. 13–15

$$
\gamma_\mathrm{smallq}(M_1<1.2)=0.3
\tag{13}
$$

$$
\gamma_\mathrm{smallq}(M_1=3.5)=\begin{cases}0.2 & \log P<2.5\\ 0.2-0.3(\log P-2.5) & 2.5\le\log P<5.5\\ -0.7-0.2(\log P-5.5) & 5.5\le\log P<8\end{cases}
\tag{14}
$$

$$
\gamma_\mathrm{smallq}(M_1>6)=\begin{cases}0.1 & \log P<1\\ 0.1-0.15(\log P-1) & 1\le\log P<3\\ -0.2-0.5(\log P-3) & 3\le\log P<5.6\\ -1.5 & 5.6\le\log P<8\end{cases}
\tag{15}
$$

`gamma_smallq`: again `g12`, `g35` and `g6`, interpolated in the same way.

---

## 2. The mass-ratio distribution $p_q$

### 2.1 Paper definition (Sect. 2, Eqn. 2 and Fig. 2)

$p_q\propto q^\gamma$ is a broken power law with

- $\gamma=\gamma_\mathrm{smallq}$ for $0.1<q<0.3$;
- $\gamma=\gamma_\mathrm{largeq}$ for $0.3<q<1$;

and it is **continuous at $q=0.3$**. On top of the power law there is an *excess* fraction
$\mathcal{F}_\mathrm{twin}$ of systems spread uniformly over $0.95<q<1$.
$p_q$ is normalised so that $\int_{0.3}^{1}p_q\,dq=1$.

Written out, with
$c=(1-\mathcal{F}_\mathrm{twin})\big/\int_{0.3}^{1}q^{\gamma_l}\,dq$:

$$
p_q(q)=
\begin{cases}
c\;0.3^{\,\gamma_l-\gamma_s}\;q^{\gamma_s} & 0.1\le q<0.3\\[4pt]
c\;q^{\gamma_l}\;+\;\dfrac{\mathcal{F}_\mathrm{twin}}{0.05}\,\mathbb{1}[q\ge0.95] & 0.3\le q\le1
\end{cases}
$$

- The factor $0.3^{\gamma_l-\gamma_s}$ makes the two power laws meet at $q=0.3$.
- The large-q part integrates to $(1-\mathcal{F}_\mathrm{twin})+\mathcal{F}_\mathrm{twin}=1$.

### 2.2 The three pieces in closed form (`q_components`)

```python
c = (1.0 - ft) / _powint(Q_SPLIT, 1.0, gl)
w_small = c * Q_SPLIT**(gl - gs) * _powint(Q_MIN, Q_SPLIT, gs)
return w_small, 1.0 - ft, ft, gs, gl, c        # w_s, w_l, w_t, ...
```

with the power-law integral (`_powint`)

$$
I(a,b,\gamma)=\int_a^b x^\gamma\,dx=
\begin{cases}\dfrac{b^{\gamma+1}-a^{\gamma+1}}{\gamma+1} & \gamma\ne-1\\[6pt] \ln(b/a) & \gamma=-1\end{cases}
$$

The three pieces and their masses are:

| Piece | Range | Density | Mass |
|---|---|---|---|
| small-q | 0.1–0.3 | $c\,0.3^{\gamma_l-\gamma_s}q^{\gamma_s}$ | $w_s=c\,0.3^{\gamma_l-\gamma_s}I(0.1,0.3,\gamma_s)$ |
| large-q | 0.3–1 | $c\,q^{\gamma_l}$ | $w_l=1-\mathcal{F}_\mathrm{twin}$ |
| twin | 0.95–1 | $\mathcal{F}_\mathrm{twin}/0.05$ | $w_t=\mathcal{F}_\mathrm{twin}$ |

### 2.3 The factor $\mathcal{Q}$: from $q>0.3$ to $q>0.1$

The paper defines companion frequencies for $q>0.3$, and converts them to $q>0.1$ by
integrating $p_q$ down to 0.1 (Sect. 9.3):

$$
f_{\log P;q>0.1}=f_{\log P;q>0.3}\;\mathcal{Q},\qquad
\mathcal{Q}=\int_{0.1}^{1}p_q\,dq = w_s+w_l+w_t = 1 + w_s .
$$

In the code: `q_factor` returns $\mathcal{Q}$, and `f_logP_q01` multiplies
`f_logP_q03` by it.

#### Why normalise over $0.3<q<1$, and where the $0.1<q<0.3$ part goes

The normalisation $\int_{0.3}^{1}p_q\,dq=1$ looks odd, because it leaves the small-q
piece outside the normalisation. It follows from what the paper measures. The robust
observable is $f_{\log P;q>0.3}$, the number of companions per decade of period with
$q>0.3$, since most surveys are complete down to $q\approx0.3$. Normalising $p_q$ over
0.3–1 keeps that number fixed whatever the slopes are.

The small-q piece is then **added on top**, so $p_q$ over $0.1<q<1$ integrates to
$\mathcal{Q}=1+w_s>1$. In other words, $p_q$ is a companion *count* relative to the
$q>0.3$ companions, not a probability density over 0.1–1. Its extra mass shows up in
the companion frequency: $f_{\log P;q>0.3}\,w_s$ companions per decade with
$0.1<q<0.3$ come in addition to the $q>0.3$ ones.

The paper's own example (Sect. 9.3) shows this. For $f_{\log P;q>0.3}=0.14$ and a
flat distribution ($\gamma_s=\gamma_l=\mathcal{F}_\mathrm{twin}=0$),
$w_s=0.2/0.7\approx0.29$, so $\mathcal{Q}\approx1.29$ and $f_{\log P;q>0.1}=0.18$.

If $p_q$ were instead normalised over 0.1–1, a steeper $\gamma_\mathrm{smallq}$ would move
companions out of $q>0.3$. That would change $f_{\log P;q>0.3}$, the quantity that
was measured and fitted.

The sampler needs a proper probability density, so $p_q$ is used twice:

1. **In the period density and binary fraction** (Sect. 4.3), as
   $f_{\log P;q>0.3}\,\mathcal{Q}$. More small-q companions mean a higher binary
   fraction, and more weight at periods where $\mathcal{Q}$ is large: long periods
   for massive primaries, where $\gamma_\mathrm{smallq}$ is strongly negative.
2. **When drawing q** for a binary at a given $P$, as $p_q/\mathcal{Q}$, which
   integrates to 1 over $0.1<q<1$. `q_pdf` divides by the summed masses, and
   `_sample_q` draws $u\in[0,w_s+w_l+w_t)$.

In the joint density the two factors of $\mathcal{Q}$ cancel:
$(f_{\log P;q>0.3}\,\mathcal{Q}\,p_\mathrm{inner})\times(p_q/\mathcal{Q})=f_{\log P;q>0.3}\,p_\mathrm{inner}\,p_q$.
This is the form used in Sect. 0.

### 2.4 Low-mass extension (not in the paper)

Below $0.8\,M_\odot$ the q distribution is truncated at $q_\mathrm{min}=M_\mathrm{min}/M_1$
and renormalised; $M_\mathrm{min}=0.08\,M_\odot$ is the minimum stellar mass. The code
cuts each of the three pieces at $q_\mathrm{min}$ (`Moe_17_PsandQs._q_masses`).

---

## 3. Eccentricity

### 3.1 Maximum eccentricity, Eqn. 3

$$
e_\mathrm{max}(P)=1-\left(\frac{P}{2\,\mathrm{d}}\right)^{-2/3}\quad(P>2\,\mathrm{d}),
\qquad\text{all binaries with }P\le2\,\mathrm{d}\text{ are circular.}
\tag{3}
$$

- `e_max` returns 0 for $P\le2$ d.
- `is_circular` is true for $P\le2$ d, and also just above 2 d, where
  $0.8\,e_\mathrm{max}$ would fall below the lower cut-off $e_\mathrm{min}$ of Sect. 3.3.

### 3.2 Power-law slope $\eta$, Eqns. 17–18

$$
\eta(M_1<3)=0.6-\frac{0.7}{\log P-0.5}\quad(0.5<\log P<6)
\tag{17}
$$

$$
\eta(M_1>7)=0.9-\frac{0.2}{\log P-0.5}\quad(0.5<\log P<5)
\tag{18}
$$

These are interpolated in $\log M_1$ across 3–7 $M_\odot$. Both fits diverge at $\log P=0.5$.
**The code holds $\eta$ at its $\log P=0.7$ value below 0.7** ($-2.9$ late-type,
$-0.1$ early-type). This choice comes from Moe's code; the paper does not specify it.

`eta`:
```python
x = np.maximum(lp, 0.7) - 0.5
e3 = 0.6 - 0.7 / x                    # Eqn. 17
e7 = 0.9 - 0.2 / x                    # Eqn. 18
return e3 + (e7 - e3) * _wlog(m, 3.0, 7.0)
```

### 3.3 The eccentricity distribution $p_e$

**Paper (Sect. 2):** $p_e\propto e^\eta$ on $0<e<e_\mathrm{max}$.

**The code** adds two things from Moe's own code:

- **a linear turnover:** $p_e\propto e^\eta$ on $[e_\mathrm{min},\,0.8\,e_\mathrm{max}]$, then a
  straight line down to 0 at $e_\mathrm{max}$ (continuous at $0.8\,e_\mathrm{max}$);
- **a lower cut-off $e_\mathrm{min}=10^{-4}$** (`E_MIN`): for $\eta\le-1$ (short periods)
  $e^\eta$ cannot be normalised from $e=0$.

$$
p_e(e)\propto
\begin{cases}
e^\eta & e_\mathrm{min}\le e<e_8\\[2pt]
e_8^{\,\eta}\,\dfrac{e_\mathrm{max}-e}{e_\mathrm{max}-e_8} & e_8\le e\le e_\mathrm{max}
\end{cases},
\qquad e_8\equiv0.8\,e_\mathrm{max}.
$$

The normalisation is $A_\mathrm{pow}+A_\mathrm{tri}$, with
$A_\mathrm{pow}=I(e_\mathrm{min},e_8,\eta)$ and
$A_\mathrm{tri}=\tfrac12(e_\mathrm{max}-e_8)\,e_8^{\,\eta}$ (`e_pdf`).
Circular binaries have $e=0$ exactly, as the paper states.

---

## 4. Period distribution and binary fraction

### 4.1 Companion frequency for $q>0.3$, Eqns. 20–23

The frequency at three anchor periods depends only on $M_1$:

$$
f_{\log P<1;q>0.3}=0.020+0.04\log M_1+0.07(\log M_1)^2
\tag{20}
$$

$$
f_{\log P=2.7;q>0.3}=0.039+0.07\log M_1+0.01(\log M_1)^2
\tag{21}
$$

$$
f_{\log P=5.5;q>0.3}=0.078-0.05\log M_1+0.04(\log M_1)^2
\tag{22}
$$

Writing $f_1, f_{2.7}, f_{5.5}$ for these, with $\alpha=0.018$ and $\Delta=0.7$
(`ALPHA`, `DLOGP`):

$$
f_{\log P;q>0.3}(M_1,\log P)=
\begin{cases}
f_1 & 0.2\le\log P<1\\[2pt]
f_1+\dfrac{\log P-1}{1.7-\Delta}\left(f_{2.7}-f_1-\alpha\Delta\right) & 1\le\log P<2.7-\Delta\\[8pt]
f_{2.7}+\alpha(\log P-2.7) & 2.7-\Delta\le\log P<2.7+\Delta\\[4pt]
f_{2.7}+\alpha\Delta+\dfrac{\log P-2.7-\Delta}{2.8-\Delta}\left(f_{5.5}-f_{2.7}-\alpha\Delta\right) & 2.7+\Delta\le\log P<5.5\\[8pt]
f_{5.5}\,e^{-0.3(\log P-5.5)} & 5.5\le\log P<8
\end{cases}
\tag{23}
$$

`f_logP_q03` implements Eqns. 20–22 as `f1`, `f27`, `f55`, and Eqn. 23 with `np.select`.
Below $\log P=0.2$ the code extends $f_1$ flat; the paper's fit starts at 0.2. The
default `orbital_period_min` is $10^{0.15}$ d, as in Moe's code. A shorter minimum
period (POSYDON's ini default is 0.75 d) extends the flat part further; see Sect. 4.3.

### 4.2 Inner-binary probability (not in the paper)

The paper counts *all* companions. Its multiplicity frequency
$f_\mathrm{mult}=\int f_{\log P}\,d\log P$ (Eqn. 1) reaches about 2 for O stars. It
estimates single, binary, triple and quadruple fractions with a Poisson model
(Eqns. 28–29, Fig. 39).

The sampler draws **one** companion, the inner binary. To do so it uses an
approximation from Moe's code that the paper does not contain (`prob_inner_binary`):

$$
p_\mathrm{inner}(M_1,\log P)=
\begin{cases}1 & \log P\le1.5\\ \max\!\left[0,\;1-0.11(\log P-1.5)^{1.43}(M_1/10)^{0.56}\right] & \log P>1.5\end{cases}
$$

According to Moe, it was tuned so that the binary fraction matches the multiplicity
statistics of Sect. 9.4 (41 % at 1 $M_\odot$, 96 % at 28 $M_\odot$). It is within about
0.03 of the paper's Poisson-model fraction of primaries with at least one companion,
$1-\mathcal{F}_{n=0}$.

### 4.3 Inner-binary period density and binary fraction

$$
g(M_1,\log P)=f_{\log P;q>0.3}\,\mathcal{Q}\,p_\mathrm{inner},\qquad
f_\mathrm{bin}(M_1)=\int_{\log P_\mathrm{min}}^{\log P_\mathrm{max}} g\,d\log P .
$$

- $g$ is `inner_binary_density`.
- $f_\mathrm{bin}$ is integrated once in `Moe_17_PsandQs.__init__`, on 600 nodes in
  $\log M_1$ (plus the kinks of the mass interpolation) and 4001 nodes in $\log P$.
  `binary_fraction` interpolates it in $\log M_1$; the table and the sampled
  fraction agree to about $10^{-4}$.
- A primary whose inner companion would lie outside
  [`orbital_period_min`, `orbital_period_max`] is a single star.

**Binary fraction above 1.** For `orbital_period_min` below about 1.4 d, the flat
extension of Eqn. 23 makes the integral exceed 1 for massive primaries (1.075 at
40 $M_\odot$ for 0.75 d). The code then caps $f_\mathrm{bin}$ at 1 and divides the
period density by the uncapped integral, which keeps its shape. It also issues an
`InappropriateValueWarning`.

The low-mass extension multiplies $f_\mathrm{bin}$ by
$\mathrm{interp}(\log M_1;\,[\log M_\mathrm{min},\log0.8],\,[0,1])$ below $0.8\,M_\odot$
(`_lowmass_scale`).

---

## 5. Sampling

The code uses **inverse-transform sampling**: if $F$ is the CDF of $x$ and $u\sim U(0,1)$, then
$x=F^{-1}(u)$ has density $F'$. Every CDF is written down analytically, at the drawn
$(M_1,P)$ of each system; there are no grids. `__call__` first decides which systems
are binaries ($u<f_\mathrm{bin}(M_1)$, or all of them with `all_binaries=True`). It
then draws $\log P$, then $q$ and $e$ given that period.

**The basic piece: a truncated power law** (`_powinv`). For $p(x)\propto x^\gamma$ on
$[a,b]$, with $k=\gamma+1$:

$$
F(x)=\frac{x^k-a^k}{b^k-a^k}
\quad\Longrightarrow\quad
x=F^{-1}(v)=\left(a^k+v\,(b^k-a^k)\right)^{1/k},
\qquad (\gamma=-1:\; x=a\,(b/a)^v).
$$

**Mass ratio** (`_sample_q`). $p_q$ is a mixture of the three pieces of Sect. 2.2.
Sampling a mixture takes two uniform numbers:

1. **Choose a piece.** Draw $u\sim U(0,\,w_s+w_l+w_t)$: small-q if $u<w_s$, large-q if
   $u<w_s+w_l$, twin otherwise. Each piece is chosen with probability equal to its mass.
2. **Place q within the piece.** Draw $v\sim U(0,1)$ and use
   - $q=F^{-1}_{[0.1,0.3],\gamma_s}(v)$ for the small-q piece,
   - $q=F^{-1}_{[0.3,1],\gamma_l}(v)$ for the large-q piece,
   - $q=0.95+0.05\,v$ for the twins.

This gives exactly $p_q$, because the density of a mixture is the mass-weighted sum of
its pieces, and the twin piece overlaps the large-q piece additively, as in Sect. 2.1.
For $M_1<0.8\,M_\odot$ the lower bounds become $\max(a,q_\mathrm{min})$ and the masses
shrink accordingly (`_q_masses`).

```python
u = rng.random(M1.size) * (w_s + w_l + w_t)
v = rng.random(M1.size)
q = np.where(u < w_s, _powinv(v, q_lo_small, Q_SPLIT, gs),
     np.where(u < w_s + w_l, _powinv(v, q_lo_large, 1.0, gl),
              lo_twin + v * (1.0 - lo_twin)))
```

**Eccentricity** (`_sample_e`). A two-piece mixture: the power law
$[e_\mathrm{min},e_8]$ with mass $A_\mathrm{pow}$, and the triangle $[e_8,e_\mathrm{max}]$ with
mass $A_\mathrm{tri}$. The triangle's density is $\propto(e_\mathrm{max}-e)$, so

$$
F(e)=1-\left(\frac{e_\mathrm{max}-e}{e_\mathrm{max}-e_8}\right)^2
\quad\Longrightarrow\quad
e=e_\mathrm{max}-(e_\mathrm{max}-e_8)\sqrt{1-v}
$$

The code uses $\sqrt v$ instead of $\sqrt{1-v}$; both are uniform. Circular binaries
get $e=0$.

**Period** (`_sample_logP`). $g(\log P)$ contains the piecewise Eqn. 23, $\mathcal{Q}(\log P)$
(itself a function of $\gamma_s,\gamma_l,\mathcal{F}_\mathrm{twin}$) and the
$(\log P-1.5)^{1.43}$ term. Its integral has no clean inverse, so the code uses
**rejection sampling**:

1. Draw $\log P\sim U(\log P_\mathrm{min},\log P_\mathrm{max})$ and $y\sim U(0,H)$, where
   $H\ge\max_{\log P} g$ is the per-$M_1$ envelope (tabulated maximum + 2 %).
2. Accept if $y<g(\log P)$; redraw the rejected systems.

The acceptance probability is $g(\log P)/H$, so accepted values have density
$\propto g(\log P)$ exactly. If $g$ ever exceeds $H$, the code raises an error rather
than returning biased samples.

All random numbers come from the `RNG` passed to `__call__`, or from the generator's
own `RNG` if none is passed.

---

## 6. The PDF

The sampler draws $\log P$ from $g$, then $q$ from $p_q(\cdot\mid\log P)$, then $e$ from
$p_e(\cdot\mid\log P)$. The joint density is therefore the product of the same
formulas, **without** the inversion step:

$$
p(\log P,q,e\mid M_1)=
\underbrace{\frac{g(M_1,\log P)}{\max(1,f_\mathrm{bin,raw})}}_{\texttt{logP\_pdf}}
\;\times\;
\underbrace{\frac{p_q(q)}{w_s+w_l+w_t}}_{\texttt{q\_pdf}}
\;\times\;
\underbrace{\frac{p_e(e)}{A_\mathrm{pow}+A_\mathrm{tri}}}_{\texttt{e\_pdf}}
$$

where $f_\mathrm{bin,raw}$ is the uncapped integral of $g$ (Sect. 4.3).

| Factor | Method | Normalisation |
|---|---|---|
| period | `logP_pdf` | integrates over $\log P$ to $f_\mathrm{bin}(M_1)$ (0 outside the period range) |
| mass ratio | `q_pdf` | integrates to 1 |
| eccentricity | `e_pdf` | integrates to 1; circular binaries have $e=0$ exactly and `pdf` leaves this factor out |
| joint | `pdf` | integrates to $f_\mathrm{bin}(M_1)$; singles carry $1-f_\mathrm{bin}$ |

The densities are per decade of $P$, the same convention as the other POSYDON period
distributions (`Sana12Period.pdf`, `PowerLawPeriod.pdf`).

Note that $g/\mathcal{Q}\cdot p_q = f_{\log P;q>0.3}\,p_\mathrm{inner}\,p_q$ with the paper's
normalisation of $p_q$. Section 0 writes it that way; the code keeps $\mathcal{Q}$
inside $g$ and divides $p_q$ by it.

The unit tests check that `pdf` integrates to `binary_fraction`, including the
low-mass and capped cases. They also check that the sampled log P histogram agrees
with `logP_pdf`.

---

## 7. Summary: paper vs. code

| Item | Paper | `Moe_17_PsandQs` |
|---|---|---|
| Eqns. 3, 5–7, 9–11, 13–15, 17–18, 20–23 | defined | implemented as module-level functions |
| Mass interpolation | "according to $M_1$" | linear in $\log M_1$ |
| $\eta$ for $\log P<0.7$ | undefined (diverges at 0.5) | held at the $\log P=0.7$ value |
| $e$ for $P\le2$ d | circular | $e=0$ |
| $e$ shape | $e^\eta$ on $0$–$e_\mathrm{max}$ | + linear turnover above $0.8e_\mathrm{max}$, $e\ge10^{-4}$ (Moe's code) |
| Period range | $0.2<\log P<8$ | `orbital_period_min`–`orbital_period_max` (default $10^{0.15}$–$10^8$ d), flat extension below 0.2, $f_\mathrm{bin}$ capped at 1 |
| Companions | all (triples count twice) | one inner binary via $p_\mathrm{inner}$ (Moe's code) |
| $M_1$ outside 0.8–40 | not covered | 0.8/40 $M_\odot$ statistics; below 0.8 $M_\odot$ the binary fraction falls to 0 at $M_\mathrm{min}=0.08$ and $q\ge M_\mathrm{min}/M_1$ |

---

## 8. Changes from the previous implementation

The previous `Moe_17_PsandQs` translated Moe's IDL code. It tabulated cumulative
distributions on grids of 101 (M1) × 158 (log P) × 91 (q) × 200 (e) nodes, in a
Python loop (~12 s to build), and drew each system in a Python loop with `np.interp`
(~18 µs per system). The equations were correct, but sampling from the grids
introduced several problems:

| Previous behaviour | Effect | Now |
|---|---|---|
| CDFs built with `cumsum(f) - f[0]` (right-endpoint sums) | each cell got the density at its right edge; the twin excess started at q = 0.94 instead of 0.95 (8× too much mass in 0.94–0.95), lowering the fitted $\mathcal{F}_\mathrm{twin}$ by ~15 % | closed-form CDFs |
| q and e taken from the nearest grid period and nearest grid M1 | 0.3 % of binaries had $e>e_\mathrm{max}(P)$; 6 % of $P\le2$ d binaries had $e>0.01$ | evaluated at each system's $(M_1,P)$ |
| One random number chose both binarity and period; the period CDF ended at $f_\mathrm{bin}$ | for $f_\mathrm{bin}>1$ (e.g. `orbital_period_min = 0.75`), periods above the CDF = 1 point were never drawn (log P > 3.2 at 40 M☉); with `all_binaries=True` (the old default), $1-f_\mathrm{bin}$ of systems got exactly `orbital_period_max` | binarity and period drawn separately; $f_\mathrm{bin}$ capped at 1 with the period shape kept; `all_binaries` defaults to False |
| Eccentricity turnover found by index on the e grid | `orbital_period_min = 1.0` raised a `ValueError` | analytic turnover |
| `independent_sample` cached one generator in a global, created with the first `kwargs` | later calls ignored a new `RNG` (same seed → different q, P, e) and a new period range | generators cached per period range; `RNG` passed on every call |
| `M_min` set to `primary_mass_min` | for primaries below 0.8 M☉, wrong binary fraction and q floor | stellar minimum 0.08 M☉ |
| Moe with `orbital_scheme='separation'` | periods in days were returned as separations in R☉ | `ValueError` |
| Moe periods with a non-Moe `eccentricity_scheme` | single stars got an eccentricity | NaN for single stars |
| A metallicity was drawn per system (and discarded) | extra random numbers | removed |

The method signatures changed:
- `__call__(M1, all_binaries=False, RNG=None)` returns `(M2, P, e)`.
- `__init__` takes `orbital_period_min`, `orbital_period_max`, `M_min` and `RNG`; other
  keyword arguments, such as the old grid sizes, are ignored.
- New methods: `binary_fraction`, `pdf`, `logP_pdf`, `q_pdf` and `e_pdf`.
