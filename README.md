# ship-waves
Fourier Galerkin method for ship waves predicts ship wakes generated by a moving pressure disturbance;
Toolbox contains:
shipWaves - a function which solves equations of motion for given Vessel, NumModel and Constants parameters - a whole time history is available;
shipWaves_lm - low memory version of shipWaves - only last time instant results are available;
vessel - a function creating Vessel structure with ship parameters;
constants - a function creating Constants structure with physical constants parameters;
numModel - a function createing NumModel structure with numerical model parameters;
movingPressureAmps - a function calculating moving pressure amplitudes, which are defined for particular vessel and domain parameters;
movingPressure - a function calculating moving pressure values of a vessel in physical domain;
freeSurfaceElevation - a function calculating free-surface elevation in a fluid domain for given values calculated by shipWaves or shipWaves_lm;
example_case - a numerical example of circular motion of a catamaran ship generating waves in deep water;
example_case_lm - % a numerical example of circular motion of a catamaran ship generating waves in finite-depth water - low memory version.
Reference: M. Paprota. 2023. A Fourier Galerkin method for ship waves. Ocean Engineering, 271, 113796
