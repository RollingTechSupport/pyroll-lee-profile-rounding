# PyRolL Profile Bulging / Rounding Plugin

A [PyRolL](https://pyroll.readthedocs.io) plugin predicting the free-surface ("bulge" /
rounding) shape of a profile's cross-section after a rolling pass, i.e. the part of the
boundary not in contact with the rolls. It focuses on **three-roll passes** (both flat-roll
mills and Kocks mills with a curved/grooved roll surface), implementing the analytical models
of:

- J.H. Min, H.C. Kwon, Y. Lee, J.S. Woo, Y.T. Im (2003), *"Analytical model for prediction of
  deformed shape in three-roll rolling process"*, J. Mater. Process. Technol. 140, 471-477 --
  flat-roll three-roll mills.
- S.-M. Byon, S.-R. Kim, T.-Y. Kim, Y. Lee (2017), *"An approximate model to predict the
  surface profile of material sections in a 3-roll rolling process"*, J. Mech. Sci. Technol.
  31(7), 3489-3497 -- Kocks mills (curved/grooved three-roll mills).

For completeness, it also keeps (unchanged) the round-oval-round two-roll model of Lee and
Choi (2000) / Lee and Goldhahn (2001) / Lee (2002), and a circular-arc construction for
square-diamond-square / oval-square two-roll passes.

See [`docs/docs.pdf`](docs/docs.pdf) for the full model description (equations, dispatch
table, validation figures) and [`docs/docs.tex`](docs/docs.tex) for its LaTeX source.

## Model scope

All models here predict the free-surface shape **given** the maximum width of the profile
after rolling; predicting that width (spread prediction) is a separate concern handled by
other PyRolL plugins such as
[`pyroll-wusatowski-spreading`](https://github.com/pyroll-project/pyroll-wusatowski-spreading).

## Usage

```python
from pyroll.core import Profile, Roll, ThreeRollPass, FlatGroove
import pyroll.profile_bulging  # registers the post-processors, no further setup needed

in_profile = Profile.round(diameter=35.8e-3, temperature=1323.15, strain=0,
                            material=["C45", "steel"], flow_stress=100e6, length=1,
                            density=7.5e3, specific_heat_capacity=690, thermal_conductivity=23)

roll_pass = ThreeRollPass(
    roll=Roll(groove=FlatGroove(r1=5e-3, usable_width=45e-3, pad_angle=30),
              nominal_radius=185e-3, rotational_frequency=248 / 60),
    inscribed_circle_diameter=29e-3,
)

out_profile = roll_pass.solve(in_profile)
print(out_profile.bulge_radius, out_profile.cross_section.area)
```

The predicted, post-processed profile is the **return value** of `solve()` (as with any
PyRolL post-processor); `roll_pass.out_profile` itself reflects the pre-post-processing
(un-bulged) geometry.

## Validation

The three-roll models are validated against the source papers' own figures by reproducing
their process parameters and comparing the predicted cross-section against boundary
coordinates digitized from the printed figures (`tests/data/`). Running the test suite
regenerates side-by-side comparison plots under `tests/output/` for human visual review:

```
hatch run test:all
```

See `tests/test_three_roll_min_flat.py`, `tests/test_three_roll_byon_koval.py` and
`tests/test_three_roll_bulge_formulas.py`.

## Usage of the Preconfigured Hatch Scripts

To build the docs use (needs a working LaTeX installation, incl. `siunitx`, `subcaption`,
`biblatex`/`biber` and `minted`/`pygments`)

    hatch run docs:build

To run the tests use

    hatch run test:all

The test environment is preconfigured to test with Python 3.9, 3.10 and 3.11.
The tests are skipped for the respective version, if it can not be found.
