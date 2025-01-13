# Test Library Requirements

This is the specification of the `fp_test` library.
The requirements of this software are treated differently from a software system.

1. `fp_test` shall perform range test giving the following results:

    - Report date
    - Histogram by range
    - Y ranges




Report date: 2024-12-30 17:24:01


[math_sqrt] MAXFINITE_ERR 2.22044604925031308e-16, MAXERR 2.22044604925031308e-16 at 5.32493416443430487e-44
MEAN ERR 3.96173944761102348e-08 in [-3.40282346638528860e+38, 3.40282346638528860e+38]

Relative Error Histogram of 4278190082 points
-------------------------------------------
|err| == ZERO: 95.582%
|err| <= E-15: 4.418e+00%
[INFO] Running test log ranges for math_sqrt

 math_sqrt [-inf , inf], spec-v2 format: [%+.17e (0x%08X), %+.17e (0x%08X)] -> [%+.17e, %+.17e] %s
=============================================================================================================================
[-inf (0xFF800000), -1.40129846432481707e-45 (0x80000001)] -> [+nan, +nan] NIL
[-0.00000000000000000e+00 (0x80000000), +0.00000000000000000e+00 (0x00000000)] -> [-0.00000000000000000e+00, +0.00000000000000000e+00] ZERO
[+1.40129846432481707e-45 (0x00000001), +3.40282346638528860e+38 (0x7F7FFFFF)] -> [+3.74339206650921624e-23, +1.84467429741979238e+19] FINITE
[+inf (0x7F800000), +inf (0x7F800000)] -> [+inf, +inf] INF
