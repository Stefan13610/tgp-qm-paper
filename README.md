# Theory of Generated Space — Emergent quantum mechanics

[![DOI](https://img.shields.io/badge/DOI-pending-lightgrey.svg)](#)

This repository will accompany the preprint:

> **Emergent quantum mechanics from the TGP substrate: measurement, Born rule, superposition, entanglement, spin–statistics, and decoherence from a single soliton ODE**
> M. Serafin, 2026. (Zenodo DOI to be minted on release)

It is the quantum-mechanics companion to the TGP core paper:

> **Theory of Generated Space: A minimal core of axioms, substrate, and effective field** —
> Zenodo DOI [10.5281/zenodo.19670324](https://doi.org/10.5281/zenodo.19670324),
> repository [Stefan13610/tgp-core-paper](https://github.com/Stefan13610/tgp-core-paper).

Starting from the same substrate soliton ordinary differential equation
that fixes the TGP core, we show that quantum mechanics is not
fundamental but emergent: measurement, the Born rule, superposition,
Bell-type entanglement, spin $\tfrac{1}{2}$, spin–statistics, and
decoherence all follow from the self-referential dynamics of the
substrate field $\Phi$.

**Headline numerical results** (verified, see `research/`):

- **$\hbar = \pi\,\chi\,A_{\mathrm{tail}}$** derived from the scaling
  symmetry of the soliton ODE: $A(\alpha)=A_{0}/\sqrt{\alpha}$ with
  fitted exponent $-0.49999$ to $R^{2}=0.9999999955$.
- **$\hbar(\Phi)=\hbar_{0}\sqrt{\Phi_{0}/\Phi}$** — a position-dependent
  Planck constant; testable shift $\Delta\hbar/\hbar\approx-3.5\times
  10^{-10}$ near Earth (atom interferometry at different altitudes).
- **$\Delta x\,\Delta p \geq \hbar$** from three independent derivations
  (oscillation period + Nyquist; Fisher information + Cramér–Rao;
  energy minimisation).
- **Born rule $P\propto|\psi|^{2}$** from the measurement asymmetry:
  detector back-reaction $\langle\Delta A_{\mathrm{det}}^{2}\rangle\sim
  A_{\mathrm{part}}^{2.028}$ with $R^{2}=0.99995$, CV $=2.3\%$.
- **Superposition** as the linear limit of the substrate ODE
  (error $\sim\varepsilon^{0.98}$, $R^{2}=0.999$); non-linear correction
  $\mathrm{d}A/A = 0.392\,\varepsilon^{1.01}$ sets the classical
  boundary $\varepsilon_{\mathrm{crit}}(10\%)=0.26$.
- **CHSH $=2\sqrt{2}$ (Tsirelson bound)** from a substrate phasor on
  the topological singlet $\pi_{3}(S^{3})=\mathbb{Z}$:
  $E(a,b)=-\cos(a-b)$ with $\max\text{diff}=3.3\times10^{-16}$;
  no-signaling verified.
- **Spin $\tfrac{1}{2}$** from hedgehog-vielbein winding $B=1$, universal
  across $g_{0}\in[0.3,0.95]$; analytic proof
  $B=2\int_{0}^{1}\sin^{2}(\pi u)\,du = 1$.
- **Spin–statistics automatic:** exchange phase $(-1)^{B}$;
  Fermi–Dirac and Bose–Einstein distributions derived; anyons in 2D
  from $\pi_{1}(\text{config})=\mathbb{Z}$.
- **Decoherence** from three unified channels ($\hbar(\Phi)$ damping,
  non-linear mode mixing, environmental back-reaction) plus quantum
  Darwinism as a geometric tail imprint on the substrate.

**Open problems explicitly flagged**: Lean 4 formalisation of CHSH
derivation (Q4); connection to the core-paper gauge-emergence problem
OP-8; second quantisation of $\Phi$ beyond mean-shift linearisation
(core-paper OP-10); direct detection of $\hbar(\Phi)$ altitude
dependence.

## Repository contents

```
paper/
  tgp_qm.tex                       — LaTeX source
  tgp_qm.pdf                       — compiled preprint

research/                          — numerical support cited in the paper
  qm_foundations/                  — q0: five analytic theorems (Born, hbar, scaling)
  qm_measurement/                  — q1: four scripts, three derivations of Dx*Dp >= hbar
  qm_superposition/                — q3: linearisation + NL correction
  qm_entanglement/                 — q4: CHSH = 2*sqrt(2), Tsirelson bound
  qm_spin/                         — q5: B=1 from pi_3(S^3)
  qm_statistics/                   — q6: (-1)^B exchange, FD/BE, anyons
  qm_decoherence/                  — q7: three paths + quantum Darwinism
```

Each `research/<folder>/` contains its own `README.md` describing the
sub-problem, scripts and outcomes.

## Build

```
cd paper
pdflatex tgp_qm.tex
pdflatex tgp_qm.tex        # second pass for cross-references
pdflatex tgp_qm.tex        # third pass if needed
```

Requirements: any modern LaTeX distribution with `amsmath`, `amssymb`,
`amsthm`, `mathtools`, `longtable`, `enumitem`, `hyperref`, `caption`.

## Reproducibility of the numerical support

All Python scripts under `research/` run standalone:

```
python research/qm_foundations/q0_analytical.py
python research/qm_foundations/q0_hbar_scaling.py
python research/qm_measurement/q1_self_referential.py
python research/qm_measurement/q1_back_reaction.py
python research/qm_measurement/q1_born_detector.py
python research/qm_measurement/q1_uncertainty_bound.py
python research/qm_superposition/q3_superposition.py
python research/qm_entanglement/q4_entanglement.py
python research/qm_entanglement/q4_substrate_bell_chsh.py
python research/qm_spin/q5_spin.py
python research/qm_statistics/q6_statistics.py
python research/qm_decoherence/q7_decoherence.py
```

Only `numpy`, `scipy` and `matplotlib` are required.

## Relationship to the core paper, sector closures, and the full workshop

This repository is a deliberately narrow slice:

- The **core paper** (axioms, substrate, effective metric, PPN
  parameters, gravitational-wave propagation, 15 proofs + 10 open
  problems) lives at
  [Stefan13610/tgp-core-paper](https://github.com/Stefan13610/tgp-core-paper)
  and is cited here via its Zenodo DOI.
- The **particle-sector closure paper** (charged leptons, Koide,
  Cabibbo) lives at
  [Stefan13610/tgp-leptons-paper](https://github.com/Stefan13610/tgp-leptons-paper)
  (Zenodo DOI [10.5281/zenodo.19706861](https://doi.org/10.5281/zenodo.19706861))
  and is independent of the present work.
- The **superconductivity closure paper** ($T_{c}$ across five families)
  lives at
  [Stefan13610/tgp-sc-paper](https://github.com/Stefan13610/tgp-sc-paper)
  (Zenodo DOI [10.5281/zenodo.19670557](https://doi.org/10.5281/zenodo.19670557)).
- The **full TGP research workshop** — long-form companion manuscript,
  cosmology channels (Hubble, DESI, S8), neutrino MSW, galaxy / mass
  scaling, UV completion — is kept separately at
  [Stefan13610/TGP](https://github.com/Stefan13610/TGP).

That is where development happens. This repository is the stable,
paper-aligned snapshot for the emergent-QM sector only.

## Citation

Citation keys will be finalised on release. Provisional BibTeX:

```bibtex
@misc{Serafin2026TGPQM,
  author       = {Serafin, Mateusz},
  title        = {{Emergent quantum mechanics from the TGP substrate:
                   measurement, Born rule, superposition, entanglement,
                   spin--statistics, and decoherence from a single
                   soliton ODE}},
  year         = {2026},
  publisher    = {Zenodo},
  note         = {DOI to be minted on release}
}

@misc{Serafin2026TGPCore,
  author       = {Serafin, Mateusz},
  title        = {{Theory of Generated Space: A minimal core of axioms,
                   substrate, and effective field}},
  year         = {2026},
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.19670324},
  url          = {https://doi.org/10.5281/zenodo.19670324}
}
```

## License

The paper text and the accompanying numerical code are released under
[CC BY 4.0](LICENSE).
