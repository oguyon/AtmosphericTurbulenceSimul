# AtmosphereModel

**Architectural Layer**: Level 1 (Physical Atmospheric Environment & Dispersion)

## Purpose

The `AtmosphereModel` module computes thermodynamic, chemical composition, and optical refractive index profiles for Earth's atmosphere from sea level to 100 km altitude. It integrates the NRLMSISE-00 empirical atmospheric model with multi-species Lorentz-Lorenz dispersion relations to evaluate polychromatic refractivity and optical ray bending.

## Components

- `AtmosphereModel.c` / `AtmosphereModel.h`: Public API orchestrator, global atmospheric state allocation, and milk CLI command registration.
- `atmmod_types.h`: Internal data structures, shared atmospheric parameters, and helper prototypes.
- `atmmod_ria_table.c`: File parser and memory loader for Refractive Index and Absorption (RIA) species tables.
- `atmmod_gas_mixture.c`: Multi-species Lorentz-Lorenz refractive index evaluator for 14 atmospheric gases ($\text{N}_2$, $\text{O}_2$, $\text{Ar}$, $\text{H}_2\text{O}$, $\text{CO}_2$, $\text{Ne}$, $\text{He}$, $\text{CH}_4$, $\text{Kr}$, $\text{H}_2$, $\text{O}_3$, $\text{N}$, $\text{O}$, $\text{H}$).
- `atmmod_standard_model.c`: Standard atmospheric profile generator calling NRLMSISE-00 (`gtd7`), water vapor saturation calculations (IAPWS-95 formulation), profile saving (`AtmosphereModel_save_stdAtmModel`), and loading (`AtmosphereModel_load_stdAtmModel`).
- `atmmod_refraction.c`: Numerical ray tracing through stratified spherical atmospheric shells to determine angular deviation and total flux transmission.
- `atmmod_config_parser.c`: Configuration file parser (`AtmosphereModel_Create_from_CONF`) to instantiate a site model from disk configuration.

## Public Headers

- `AtmosphereModel.h`:
  - `int init_AtmosphereModel(void)`
  - `double AirMixture_N(...)`
  - `float AtmosphereModel_stdAtmModel_N(float alt, float lambdaum, int mode)`
  - `double AtmosphereModel_H2O_Saturation(double T)`
  - `int AtmosphereModel_save_stdAtmModel(char *fname)`
  - `int AtmosphereModel_build_stdAtmModel(char *fname)`
  - `int AtmosphereModel_load_stdAtmModel(char *fname)`
  - `int AtmosphereModel_Create_from_CONF(char *CONFFILE, float slambda)`
  - `double AtmosphereModel_RefractionPath(double lambda, double Zangle, int WritePath)`
  - `int ATMOSPHEREMODEL_loadRIA_readsize(char *fname)`
  - `int ATMOSPHEREMODEL_loadRIA(char *fname, double *lptr, double *RIptr, double *absptr)`

## Dependencies

- `CLIcore` (command parsing and configuration parameter extraction)
- `OpticsMaterials` (material dispersion and refractive index lookups)
- Vendored `nrlmsise-00.20131225` (empirical density and temperature model)
- C Standard Library (`math.h`, `stdio.h`, `stdlib.h`, `string.h`)
