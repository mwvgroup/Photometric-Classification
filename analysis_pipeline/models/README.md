# SN 1991bg Source Classes

This directory provides `sncosmo.Source` objects for a custom SN 1991bg-like 
model. Multiple versions of source classes are provided. Each version is
contained in it's own file entitled with the version name of the contained
source.

## File List:

| File Name | Description  |
|-----------|--------------|
| template.npy | The wavelength and phase dependent flux values of the unproved 91bg model. |
| _utils.py | Various utilities shared by different source versions. |
| _color_interpolation.py | The original ported version of the 91bg model into SNCosmo. Flux is determined by using a rectangular bivariate spline to for a given stretch value and then linearly interpolating for color. Stretch values are internally handeled as normalized to a value of 0.65. |
| _salt2_phase.py | The same as the `color_interpolation` version, but with Stretch values are internally normalized to a value of 1 and the phase range limited to only cover -18 to 50 days. This more closely matches the salt2 model which covers -20 to 50 days. |

