---
description: Audit type consistency in math operations and loops.
---

# Check Type Consistency

Verify numerical and type consistency in `milkatmturb`:
- [ ] Loop index types match bound types (e.g. `long ii = 0; ii < npupil; ii++`).
- [ ] Math functions match datatypes (`sqrtf`, `sinf` for `float`; `sqrt`, `sin` for `double`).
- [ ] Complex FFTW types match single vs double: `fftwf_complex` and `fftwf_*` for single precision;
  `fftw_complex` and `fftw_*` for double precision.
- [ ] Multiplications that can overflow are cast to wider types (e.g. `(size_t)xsize * ysize`).
- [ ] No implicit double promotions in single-precision wavefront calculations.
