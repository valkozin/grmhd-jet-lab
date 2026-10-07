# JET — GRMHD data viewer

The primary page now replays real computed Kerr GRMHD snapshots. Axisymmetric iharm2d_v4 96×96 ran to t=1500; the separate HARMPI 64×32×16 3D pilot is saved to t=150.012. Both are coarse and not convergence-validated. Source code, selected unmodified raw snapshots, checkpoints, diagnostics, exact commits, local HARMPI fixes and reproduction instructions are in dist/grmhd-repro.zip. The 3D pilot was stopped before its configured tf=300 and is not an established jet.

The renderer is a coordinate-space volume view with illustrative emission, not a relativistic radiative-transfer image. GLSL and interactions have not been visually tested in a browser in this environment. JS syntax and exported data were checked.

Previous browser MHD model preserved under dist/legacy.

# JET — axisymmetric magnetic-tower laboratory

A browser-local Newtonian ideal-MHD experiment, NOT a GRMHD black-hole simulation. No build or external JS dependencies. Serve `dist/` as static files. Worker computes finite-volume evolution; WebGL2 displays a revolved slice with illustrative emissivity. The interface and full physical limitations are in English.

## Model
Cylindrical (r,z), all velocity and magnetic components, gamma=5/3, mu0=1. Rusanov flux, minmod reconstruction, first-order Euler time integration, CFL=.27 with a two-direction bound, GLM divergence cleaning. Geometric cylindrical sources included. Uniform initial atmosphere and vertical B. Prescribed odd toroidal field source represents an unresolved engine; consistent incremental magnetic energy and a small mass/thermal source. No injected axial momentum. No gravity, Kerr metric, accretion disk evolution, relativity, radiation transfer, or 3D instabilities.

## Validation
`node verify.cjs` checks static uniform equilibrium; finite driven runs; reflection symmetry; outward mass flux; absence of pressure-floor energy in selected runs; engine shutoff. These are regression sanity checks, not a convergence study or astrophysical validation. Static JS syntax was checked. A numerical snapshot was inspected as scalar-field plots. Browser/WebGL rendering and interaction QA was unavailable in the authoring environment.

Open boundaries can reflect waves. Pressure floors and GLM cleaning have conservation limitations. Grid results are qualitative, and the central schematic dark region is not a calculated black-hole shadow. Rendered brightness is not a measured spectrum. See the in-app model panel for equations, normalization, sources, diagnostics and references.

## Run locally

From the repository root:

```sh
python3 -m http.server 8000 --directory dist
```

Open http://localhost:8000 in a WebGL2-capable browser. No build step is required. The primary viewer replays precomputed data; it does not run GRMHD in the browser. The legacy MHD solver runs locally in a Web Worker. Reproduction scripts and third-party source licenses are inside `dist/grmhd-repro.zip`.
