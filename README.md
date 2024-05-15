# Binaural Rendering of Baffled Microphone Arrays
This repository contains MATLAB code for the binaural rendering of spherical microphone arrays (**SMA**s), equatorial microphone arrays (**EMA**s), and non-spherical microphone arrays (**XMA**s).

This version of the code (release [**v2024.JAES**](https://github.com/HaHeho/baffled-arrays-to-binaural/releases/tag/v2024.JAES)) was utilized in the manuscript
> TODO_CITATION
> <br>
> [[TODO_pdf]](TODO_PDF)
> [[audio examples]](http://www.ta.chalmers.se/research/audio-technology-group/audio-examples/jaes-2024a/)
> [[data]](https://zenodo.org/records/8206571)
> [[experiment resources]](https://doi.org/10.5281/zenodo.10901444)
* [README](https://github.com/HaHeho/baffled-arrays-to-binaural/tree/v2024.JAES)
* [Download](https://github.com/HaHeho/baffled-arrays-to-binaural/releases/tag/v2024.JAES)
* <details>
  <summary>Clone</summary>
  
  ```
  git clone --recurse-submodules https://github.com/HaHeho/baffled-arrays-to-binaural.git --branch v2024.JAES
  ```
</details>

An older version of the code (release [**v2023.DAGA**](https://github.com/HaHeho/baffled-arrays-to-binaural/releases/tag/v2023.DAGA)) was utilized in the manuscript
> H. Helmholz, T. Deppisch, and J. Ahrens, “End-to-End Magnitude Least Squares Binaural Rendering for Equatorial Microphone Arrays,” in Fortschritte der Akustik -- DAGA 2023, 2023, pp. 1679–1682.
> <br>
> [[pdf]](https://research.chalmers.se/publication/535525/file/535525_Fulltext.pdf)
> [[audio examples]](http://www.ta.chalmers.se/research/audio-technology-group/audio-examples/daga-2023a/)
* [README](https://github.com/HaHeho/baffled-arrays-to-binaural/tree/v2023.DAGA)
* [Download](https://github.com/HaHeho/baffled-arrays-to-binaural/releases/tag/v2023.DAGA)
* <details>
  <summary>Clone</summary>
  
  ```
  git clone --recurse-submodules https://github.com/HaHeho/baffled-arrays-to-binaural.git --branch v2023.DAGA
  ```
</details>

## Setup
Clone the repository locally. Make sure to include the submodules with
```
git clone --recurse-submodules https://github.com/HaHeho/baffled-arrays-to-binaural.git
```

Alternatively, download the code and manually add the following additional tools to the [`dependencies/`](dependencies) folder:
* [AKtools](https://github.com/f-brinkmann/AKtools.git)
* [eMagLS rendering scripts](https://github.com/thomasdeppisch/eMagLS.git)
* [Ambisonic Encoding toolbox](https://github.com/AppliedAcousticsChalmers/ambisonic-encoding.git) (included in eMagLS)
* [Array Response Simulator toolbox](https://github.com/polarch/Array-Response-Simulator.git) (included in eMagLS)
* [Spherical Harmonic Transform Library](https://github.com/polarch/Spherical-Harmonic-Transform.git) (included in eMagLS)
* [SOFA Toolbox](https://github.com/sofacoustics/SOFAtoolbox.git)

---
Retrieve the following additional dependencies via the MATLAB Add-Ons menu or from Mathworks File Exchange:
* `natsortfiles()` [[file exchange]](https://se.mathworks.com/matlabcentral/fileexchange/47434-natural-order-filename-sort)

Furthermore, manually add the necessary acoustic impulse responses (external source files are subject to their corresponding licenses) to the project structure:
* Microphone array impulse responses (see [`resources/ARIR_processed/README`](resources/ARIR_processed/README.md))
* **XMA** calibration and equalization filters (see [`resources/ARIR_processed/README`](resources/ARIR_processed/README.md))
* Head-related impulse responses (see [`resources/HRIR_KEMAR/README`](resources/HRIR_KEMAR/README.md))

## Usage
Script [`x4_Preprocess_HRTFs.m`](x4_Preprocess_HRTFs.m):
* Pre-process an HRTF data set in various ways before employing it for binaural rendering. At this point, we utilize a spatially subsampled set of the published HRTF to lower computation resources for all binaural renderings without losing spatial resolution.
* The HRTFs are circ-shifted so that the signals start with the lowest energy sample out of all incidence directions and both ears. The data is then time-windowed to fade the start and the end of the HRTFs to further limit artifacts from time discontinuities. The resulting data set is exported with the suffix `_adjusted`, then utilized in the remainder of the binaural rendering.
* Additional versions of the data set are exported employing direction-independent equalizations by the HRTFs diffuse-field response (suffix `_diffuse`), by the former superimposed with the Harman curve (suffix `_harman`), and by the diffuse-field response of another measurement of the same HRTF (suffix `_reference`). Thereby, respective minimum-phase and linear-phase variants of the generated equalization filters are exported as impulse responses and plotted.

Script [`x5_Render_Arrays.m`](x5_Render_Arrays.m):
* Perform binaural rendering of SMA, EMA and XMA impulse response data sets at a specified spherical harmonics order and with a desired HRTF set.
* Thereby, a reference binaural room impulse response data set (either from a dummy head measurement or from a former high-resolution array rendering) can be specified. The reference data is compared against the currently rendered binaural ear signals by generating extensive plots to evaluate time domain and frequency domain differences individually at all rendered head orientations.

Script [`x5a_Render_Arrays_Batch.m`](x5a_Render_Arrays_Batch.m):
* Consecutively perform an arbitrary combination of binaural renderings from microphone array impulse responses and varying parameter sets.

Script [`x5b_Compare_Rendering_Differences.m`](x5b_Compare_Rendering_Differences.m):
* Compare arbitrary combinations of (rendered) binaural room impulse responses by plotting various time domain and averaged frequency domain representations. This is helpful for directly comparing similar rendering configurations and validating the rendering method or detecting potentially flawed data sets.

Script [`x5c_Compare_Rendering_Levels.m`](x5c_Compare_Rendering_Levels.m):
* Compare arbitrary combinations of (rendered) binaural room impulse responses by plotting various resulting signal levels. This is helpful for the preparation of the rendered BRIRs for the perceptual comparison in a user study, where all stimuli should ideally be loudness normalized.

Script [`x5d_Compare_Rendering_EDC.m`](x5d_Compare_Rendering_EDC.m):
* Generate an interactive plot to compare arbitrary combinations of (rendered) binaural room impulse responses by their Energy Decay Curve. This is helpful for directly comparing similar rendering configurations and validating the rendering method or detecting potentially flawed data sets.

Script [`x5e_Compare_Rendering_Azimuths.m`](x5e_Compare_Rendering_Azimuths.m):
* Generate an interactive plot to compare arbitrary combinations of rendered binaural room impulse responses by their azimuth alignment angle applied during binaural rendering. This is helpful for directly comparing similar rendering configurations and validating the rendering method or detecting potentially flawed data sets.

## Contributing
See [CONTRIBUTING](CONTRIBUTING.md) and [CODE_OF_CONDUCT](CODE_OF_CONDUCT.md) for full details.

## Acknowledgment
We thank [Reality Labs Research by Meta](https://research.facebook.com/publications/research-area/augmented-reality-virtual-reality/) for funding this project.

## License
This software is licensed under a Non-Commercial Software License (see [LICENSE](LICENSE) for full details).
