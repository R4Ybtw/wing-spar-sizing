# Wing spar sizing for a small UAV

What is the lightest aluminium tube spar that lets a 4 kg UAV pull 4 g with a 1.5 factor of safety, keep tip deflection under 5% of the half-span, and still fit inside the wing? And how much does the answer depend on the aerodynamic load model?

The design chain is: spanwise lift from lifting-line theory, then bending moment from statics, then deflection from a beam finite-element model, then a sweep over tube sizes. Each link is cross-checked by an independent method:

| Link | Cross-check |
|---|---|
| Spanwise lift (lifting line) | SolidWorks Flow Simulation, with the half-wing split into ten strips |
| Deflection and stress (beam FEM) | SolidWorks Simulation FEA of the tube |
| Tube stiffness (beam theory) | A physical bending test |

*Work in progress: the load model, sizing and analysis tools are complete and verified; the beam FEM, CFD, FEA and bench-test results are being added.*

## The aircraft

| Item | Value |
|---|---|
| Mass | 4.0 kg (weight 39.2 N) |
| Wing | Rectangular, 2.4 m span, 0.25 m chord, NACA 4412 (aspect ratio 9.6) |
| Load case | Limit load factor 4 at 25 m/s; ultimate load factor 6 (limit × 1.5) |
| Spar | 6061-T6 round tube at the quarter chord. The airfoil is 29.7 mm deep there, so the tube can be at most 25.7 mm |

All inputs live in `data/aircraft.json`.

## Results with the current inputs

| Quantity | Value |
|---|---|
| Wing angle of attack at the design point | 4.1° |
| Span efficiency | 0.92 |
| Root bending moment, ultimate | 65.3 N·m (lifting line), 59.9 N·m (elliptic assumption). The rectangular wing carries more lift outboard |
| Lightest feasible spar | 25 × 0.7 mm tube, 346 g for both half-spars |
| Tip deflection at limit load | 55 mm (limit 60 mm) |
| Root bending stress at ultimate load | 207 MPa (allowable 310 MPa) |
| Binding requirement | Fit inside the airfoil (3% margin), then stiffness (9%) |

## Approach

- **Loads:** Prandtl lifting-line theory, solved by Glauert's Fourier-series method. The load factor fixes the total lift; the lifting line sets only how it is distributed along the span.
- **Statics:** shear force and bending moment are exact for a cantilever, which is statically determinate. Stresses therefore come from statics.
- **Deflection:** an Euler–Bernoulli beam FEM with Hermite cubic elements and consistent load vectors.
- **Sizing:** a sweep over tube diameters and wall thicknesses against five requirements:
  - ultimate strength
  - yield at limit load
  - tip deflection
  - local wall buckling (classical cylinder buckling stress with a 0.3 knockdown factor)
  - fit inside the airfoil
- **Bench test:** the deflection of a real tube is measured at several lengths. Fitting the results separates the tube's bending stiffness EI from the rotational flexibility of the clamp.

## Verification

`tests/test_wing.m` checks:

| Component | Checked against |
|---|---|
| Lifting line, elliptic wing | Exact theory (CL, and span efficiency e = 1) |
| Lifting line, rectangular wing | Convergence in the number of terms, and Helmbold's formula (within 2%) |
| Statics | Exact root moments for uniform and elliptic loads |
| Section properties | Closed-form values |
| Load case | Closed-form values |
| Bench-test fit | Synthetic data with a known answer |
| Beam FEM, uniform load | qL⁴/8EI, exact with a single element |
| Beam FEM, tip point load | PL³/3EI |
| Beam FEM, elliptic load | An independent virtual-work integral |
| Beam FEM, tapered beam | Second-order convergence |

## Simulation setup

| Item | Setting |
|---|---|
| CFD | Flow Simulation, 3D half-wing with a symmetry plane at the root; 25 m/s at 4.1°; separate X and Y force goals on each of the ten strips; three-level mesh study on total lift |
| FEA | SolidWorks Simulation static study of the tube: fixed root face; each strip's lift applied to its split face; second-order solid mesh at 3 mm and 2 mm. Stress is read 100 mm from the root, clear of the fixture singularity |
| Expected FEA answer (25 × 0.7 mm) | Tip deflection 55.3 mm at limit load; bending stress 171 MPa at ultimate load, 100 mm from the root |

The settings actually used for each run are recorded in `docs/run_log.md`.

## Repository layout

| Path | Contents |
|---|---|
| `data/aircraft.json` | Aircraft, load case and material |
| `src/load_aircraft.m`, `src/wing_loads.m` | Load case and spanwise lift at limit and ultimate load |
| `src/lifting_line.m` | Prandtl lifting line |
| `src/beam_fem.m` | Euler–Bernoulli beam FEM |
| `src/cantilever_statics.m`, `src/exact_tip_deflection.m` | Exact statics, and the virtual-work deflection used in testing |
| `src/size_spar.m` | Design sweep |
| `src/fea_loads_table.m` | Strip loads for the FEA, with the expected answer |
| `src/compare_cfd_loads.m` | Lifting line against the CFD strips |
| `src/bench_test.m` | Fit of measured deflections |
| `tests/test_wing.m` | Verification suite |
| `cad/` | Half-wing and spar tube split into strips (STEP), plus the script that builds them |
| `docs/run_log.md` | Run-by-run record of settings and results |

## Reproducing the analysis

```matlab
addpath tests; test_wing          % verification
run('src/size_spar.m')            % design sweep
run('src/fea_loads_table.m')      % FEA loads and expected answer
```

The code runs in MATLAB and in GNU Octave.

## Limitations

- **Beam model of the spar only:** the skin and ribs carry no load, and torsion is not modelled.
- **One load case:** a single static maneuver load, with no gust loads or fatigue.
- **No weight relief:** the wing's own weight is not subtracted, which is conservative.
- **Crude buckling check:** local buckling is estimated with a knockdown factor.

## References

- J. D. Anderson, *Fundamentals of Aerodynamics*, McGraw-Hill: lifting-line theory and Glauert's method.
- D. L. Logan, *A First Course in the Finite Element Method*, Cengage: beam elements.
- NASA SP-8007, *Buckling of Thin-Walled Circular Cylinders*: knockdown factors.
- Metallic Materials Properties Development and Standardization (MMPDS): design allowables for 6061-T6.

## Tools and AI use

MATLAB, SolidWorks Flow Simulation and SolidWorks Simulation. Apart from `src/beam_fem.m`, which is my own, the code was first drafted with an AI assistant, then read, modified and run by me. The simulations, the bench test and the conclusions are my own.

## Licence

Code: MIT (see `LICENSE`). Text and figures: CC BY 4.0.
