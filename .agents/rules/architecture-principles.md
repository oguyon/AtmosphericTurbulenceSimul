---
trigger: always_on
---

# Architecture Principles

- **Minimize Cross-Module Dependencies**: Avoid introducing circular or unstructured dependencies
  between modules (e.g. `AtmosphereModel`, `AtmosphericTurbulence`, `OpticsMaterials`,
  `WFpropagate`). Refactor towards a layered architecture where modules interact only through
  well-defined public APIs declared in their top-level headers.
- **Pure Compute vs. CLI Wrappers**: Keep simulation engines (Fourier phase synthesis, Fresnel
  propagation, linear predictors) pure and decoupled from interactive terminal interfaces or
  direct user prompt logic.
- **Adhere to Module Boundaries**: Each submodule directory contains its own `README.md` documenting
  its role, public headers, and dependencies. Respect these boundaries when adding new features.
