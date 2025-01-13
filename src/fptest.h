#ifndef FPTEST_H
#define FPTEST_H

/* Boolean */
#if ((__STDC_VERSION__ + 0) >= 199901L)
    typedef _Bool bool_t;
#else
    typedef char bool_t;
#endif

/* Character */
typedef char         char_t;
typedef const char * cstr_t;

/* Integers */
typedef unsigned char       u8_t;
typedef unsigned short      u16_t;
typedef unsigned int        u32_t;
typedef unsigned long long  u64_t;

typedef signed   char       i8_t;
typedef signed   short      i16_t;
typedef signed   int        i32_t;
typedef signed   long long  i64_t;

/* Floating Points */
typedef float               f32_t;
typedef double              f64_t;
typedef long double         f80_t;

/* Function Pointers */
typedef f32_t  ( *f3232_t   ) ( f32_t );
typedef f64_t  ( *f6464_t   ) ( f64_t );
typedef f80_t  ( *f8080_t   ) ( f80_t );
typedef void   ( *fvoid64_t ) ( f64_t );

/* Composites */
typedef struct { f32_t x0; f32_t x1;            } fp32_vec2_t;
typedef struct { f32_t x0; f32_t x1; f32_t x2;  } fp32_vec3_t;

typedef struct { f64_t x0; f64_t x1;            } fp64_vec2_t;
typedef struct { f64_t x0; f64_t x1; f64_t x2;  } fp64_vec3_t;

/* Enums */
typedef enum
{
    NAMED_FP_INF,       /* -inf */
    NAMED_FP_MAXNORM,   /* 1e+38, 1.e+308 */
    NAMED_FP_ZERO,      /* +0.0 */
    NAMED_FP_MINNORM,   /* 1e-38, 1.e-308 */
    NAMED_FP_MINSUBN,   /* 2**(emin-p) = 2**(-126-23) ~ 1.4e-45 (f32_t) , ~2.47e-324 (f64_t),  */
} named_fp_t;

/*
 * identifes the type of floating-point value according to the union of sets:
 * { NAN } U { ZERO } U { NORMAL } U { SUBNORMAL } U { INFINITE }
 */

typedef enum
{
    FP_TYPE_NAN,
    FP_TYPE_ZERO,
    FP_TYPE_NORMAL,
    FP_TYPE_SUBNORMAL,
    FP_TYPE_INFINITE,
} fpset_t;

/*
 * Triplet structure
 */

typedef struct
{
    u64_t M; /* integral significand: M = m·B^{-p+1}, with m: 1 <= m < B (if normal), or 0 < m < 1 (if subnormal) */
    i16_t e; /* exponent: emin <= e <= emax, emin = 1 - emax, with 255 reserved for f32_t qNaN, and 2047 for f64_t qNaN */
    i16_t s; /* { -1, +1 } */
} fptriplet_t;

/* You should define ADD_EXPORTS *only* when building the DLL. */
#ifdef ADD_EXPORTS
  #define ADDAPI __declspec(dllexport)
#else
  #define ADDAPI __declspec(dllimport)
#endif

/* Define calling convention in one place, for convenience. */
#define ADDCALL __cdecl

/* Make sure functions are exported with C linkage under C++ compilers. */

#ifdef __cplusplus
extern "C"
{
#endif

typedef struct
{
    bool_t debug_next_x_inside_boundaries;
    bool_t log_completation;
    bool_t log_completation_remove_history;
} fpenv_t;

extern ADDAPI fpenv_t g_fpenv;

/* Functions */
ADDAPI extern f64_t         ADDCALL fp32_ulps( f32_t approx, f32_t exact );
ADDAPI extern f64_t         ADDCALL fp64_ulps( f64_t approx, f64_t exact );

ADDAPI extern bool_t        ADDCALL fp32_equals( f32_t a, f32_t b );
ADDAPI extern bool_t        ADDCALL fp64_equals( f64_t a, f64_t b );
ADDAPI extern bool_t        ADDCALL fp32_equals_sign( f32_t a, f32_t b );
ADDAPI extern bool_t        ADDCALL fp64_equals_sign( f64_t a, f64_t b );

ADDAPI extern f32_t         ADDCALL f32_geom_step_real_line( f32_t x, f32_t frac );
ADDAPI extern f64_t         ADDCALL f64_geom_step_real_line( f64_t x, f64_t frac );

ADDAPI extern fpset_t       ADDCALL f32_get_subset( f32_t x );
ADDAPI extern fpset_t       ADDCALL f64_get_subset( f64_t x );

ADDAPI extern cstr_t        ADDCALL fp64_get_subset_name( f64_t x );
ADDAPI extern cstr_t        ADDCALL fp32_get_subset_name( f32_t x );

ADDAPI extern f32_t         ADDCALL f32_get_named_fp_in_real_line( named_fp_t point );
/* TODO: ADDAPI extern f32_t         ADDCALL f32_get_named_fp_in_real_line( named_fp_t point ); */

ADDAPI extern cstr_t        ADDCALL f32_sprint_digits_radix2( char_t buff [ 32 + 2 + 1 ], char_t sepparator, f32_t x );
ADDAPI extern cstr_t        ADDCALL f64_sprint_digits_radix2( char_t buff [ 64 + 2 + 1 ], char_t sepparator, f64_t x );

ADDAPI extern f32_t         ADDCALL f32_next_float( f32_t x );
ADDAPI extern f64_t         ADDCALL f64_next_float( f64_t x );

ADDAPI extern f32_t         ADDCALL f32_mount_bitfields( u32_t s, u32_t e, u32_t m );
ADDAPI extern f64_t         ADDCALL f64_mount_bitfields( u64_t s, u64_t e, u64_t m );

ADDAPI extern fptriplet_t   ADDCALL f32_get_triplet( f32_t x );
ADDAPI extern fptriplet_t   ADDCALL f64_get_triplet( f64_t x );

ADDAPI extern i16_t         ADDCALL fp32_get_exp( f32_t x );
ADDAPI extern i16_t         ADDCALL fp64_get_exp( f64_t x );

ADDAPI extern f32_t         ADDCALL f32_eval_triplet( fptriplet_t x );
ADDAPI extern f64_t         ADDCALL f64_eval_triplet( fptriplet_t x );

ADDAPI extern f32_t         ADDCALL f32_set_exp( f32_t x, i16_t n );
ADDAPI extern f64_t         ADDCALL f64_set_exp( f64_t x, i16_t n );

ADDAPI extern fp32_vec2_t   ADDCALL fp32_find_control_boundaries( f32_t at_x, const f32_t * control_points,  f64_t boundary_semi_length_ulps, i32_t n );
ADDAPI extern fp64_vec2_t   ADDCALL fp64_find_control_boundaries( f64_t at_x, const f64_t * control_points,  f64_t boundary_semi_length_ulps, i32_t n );

ADDAPI extern f32_t         ADDCALL fp32_next_x( f32_t x, f32_t frac, const f32_t * control_points, f64_t round_ulps, i32_t n );
ADDAPI extern f64_t         ADDCALL fp64_next_x( f64_t x, f64_t frac, const f64_t * control_points, f64_t round_ulps, i32_t n );

ADDAPI extern void          ADDCALL fp32_range_analyzer( cstr_t desc, f3232_t lhs, f3232_t rhs, f3232_t next_x, f64_t accept, f64_t reject, u8_t state[ /* TBD: Size in bytes of the state */ ] );
ADDAPI extern void          ADDCALL fp64_range_analyzer( cstr_t desc, f6464_t lhs, f6464_t rhs, f6464_t next_x, f64_t accept, f64_t reject, u8_t state[ /* TBD: Size in bytes of the state */ ] );

/* PUBLIC UTILITIES */
#define NULLPTR ((void*) 0)
#define TRUE  ((bool_t) 1)
#define FALSE ((bool_t) 0)
#define SIZEOF(array) ( sizeof( array ) / sizeof( (array)[0] ) )

#ifdef __cplusplus
} // __cplusplus defined.
#endif

#endif/*FPTEST_H*/
