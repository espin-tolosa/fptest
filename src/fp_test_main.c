#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "fp_test_internals.h"

fpenv_t g_fpenv = { 0 }; /* Initialize floating-point environment */

static f64_t    f32_ulps_impl ( f32_t approx, f80_t diff );
static f64_t    f64_ulps_impl ( f64_t approx, f80_t diff );

typedef struct
{
    i32_t xexp;
    u64_t count_points_in_the_range;
} fp_progress_state_t;

static u32_t    fp32_progress_bar_status( f32_t x );
static void     fp32_progress_bar( fpenv_t fpenv, f32_t x_curr, f32_t x_min, f32_t x_max, fp_progress_state_t * pstate );

static u32_t    fp64_progress_bar_status( f64_t x );
static void     fp64_progress_bar( fpenv_t fpenv, f64_t x_curr, f64_t x_min, f64_t x_max, fp_progress_state_t * pstate );

static void     console_erase_last_n_chars( u8_t n, fpenv_t fpenv );

#define LINE_SEPARATOR "========================================================================================="

/* >> ULPs HISTOGRAM */
#define ULPS_NHIST 20
#define ULPS_HIST_INIT \
static u64_t m_hist[ ULPS_NHIST ]; \
static f64_t m_hist_xl_negx[ ULPS_NHIST ]; \
static f64_t m_hist_xr_negx[ ULPS_NHIST ]; \
static f64_t m_hist_xl_posx[ ULPS_NHIST ]; \
static f64_t m_hist_xr_posx[ ULPS_NHIST ]; \
for(u32_t i = 0; i < SIZEOF(m_hist); i++){ m_hist[i]=0; m_hist_xl_negx[i]=0.0/0.0; m_hist_xr_negx[i]=0.0/0.0; m_hist_xl_posx[i]=0.0/0.0; m_hist_xr_posx[i]=0.0/0.0; }

#define update_hist_ranges(id,x) \
if( (x) < 0.0 ) {m_hist_xl_negx[ id ] = isnan( m_hist_xl_negx[ id ]) ? ( x ) : m_hist_xl_negx[ id ];} \
if( (x) < 0.0 ) {m_hist_xr_negx[ id ] = ( x );} \
if( (x) > 0.0 ) {m_hist_xl_posx[ id ] = isnan( m_hist_xl_posx[ id ]) ? ( x ) : m_hist_xl_posx[ id ];} \
if( (x) > 0.0 ) {m_hist_xr_posx[ id ] = ( x );}

#define ULPS_HIST_RECORD(ulp,x) \
    if      ( (ulp) ==  0.0 )    { /* c0 */ m_hist[  0 ]++; update_hist_ranges(  0, ( x ) ) } \
    else if ( (ulp)  <  1.0 )    { /* c0 */ m_hist[  1 ]++; update_hist_ranges(  1, ( x ) ) } \
    else if ( (ulp)  <  2.0 )    { /* c0 */ m_hist[  2 ]++; update_hist_ranges(  2, ( x ) ) } \
    else if ( (ulp)  <  3.0 )    { /* c0 */ m_hist[  3 ]++; update_hist_ranges(  3, ( x ) ) } \
    else if ( (ulp)  <  4.0 )    { /* c0 */ m_hist[  4 ]++; update_hist_ranges(  4, ( x ) ) } \
    else if ( (ulp)  <  5.0 )    { /* c0 */ m_hist[  5 ]++; update_hist_ranges(  5, ( x ) ) } \
    else if ( (ulp)  <  6.0 )    { /* c0 */ m_hist[  6 ]++; update_hist_ranges(  6, ( x ) ) } \
    else if ( (ulp)  <  7.0 )    { /* c0 */ m_hist[  7 ]++; update_hist_ranges(  7, ( x ) ) } \
    else if ( (ulp)  <  8.0 )    { /* c0 */ m_hist[  8 ]++; update_hist_ranges(  8, ( x ) ) } \
    else if ( (ulp)  <  9.0 )    { /* c0 */ m_hist[  9 ]++; update_hist_ranges(  9, ( x ) ) } \
    else if ( (ulp)  <  10. )    { /* c0 */ m_hist[ 10 ]++; update_hist_ranges( 10, ( x ) ) } \
    else if ( (ulp)  <  11. )    { /* c0 */ m_hist[ 11 ]++; update_hist_ranges( 11, ( x ) ) } \
    else if ( (ulp)  <  12. )    { /* c0 */ m_hist[ 12 ]++; update_hist_ranges( 12, ( x ) ) } \
    else if ( (ulp)  <  13. )    { /* c0 */ m_hist[ 13 ]++; update_hist_ranges( 13, ( x ) ) } \
    else if ( (ulp)  <  14. )    { /* c0 */ m_hist[ 14 ]++; update_hist_ranges( 14, ( x ) ) } \
    else if ( (ulp)  <  15. )    { /* c0 */ m_hist[ 15 ]++; update_hist_ranges( 15, ( x ) ) } \
    else if ( (ulp)  <  16. )    { /* c0 */ m_hist[ 16 ]++; update_hist_ranges( 16, ( x ) ) } \
    else if ( (ulp)  <  17. )    { /* c0 */ m_hist[ 17 ]++; update_hist_ranges( 17, ( x ) ) } \
    else if ( (ulp)  <  18. )    { /* c0 */ m_hist[ 18 ]++; update_hist_ranges( 18, ( x ) ) } \
    else                         { /* c0 */ m_hist[ 19 ]++; update_hist_ranges( 19, ( x ) ) }

#define get_negpos_ranges(id) m_hist_xl_negx[ id ], m_hist_xr_negx[ id ], m_hist_xl_posx[ id ], m_hist_xr_posx[ id ]

#define ULPS_HIST_REPORT \
    printf("\n Histogram\n"LINE_SEPARATOR"\n"); \
    m_hist[  0 ]    &&  printf( "|ulp| == 0.0: %10llu in [ %23.15e, %23.15e ] U [ %23.15e, %23.15e ]\n", ( m_hist[  0 ] ), get_negpos_ranges(  0 ) ); \
    m_hist[  1 ]    &&  printf( "|ulp| <= 1.0: %10llu in [ %23.15e, %23.15e ] U [ %23.15e, %23.15e ]\n", ( m_hist[  1 ] ), get_negpos_ranges(  1 ) ); \
    m_hist[  2 ]    &&  printf( "|ulp| <= 2.0: %10llu in [ %23.15e, %23.15e ] U [ %23.15e, %23.15e ]\n", ( m_hist[  2 ] ), get_negpos_ranges(  2 ) ); \
    m_hist[  3 ]    &&  printf( "|ulp| <= 3.0: %10llu in [ %23.15e, %23.15e ] U [ %23.15e, %23.15e ]\n", ( m_hist[  3 ] ), get_negpos_ranges(  3 ) ); \
    m_hist[  4 ]    &&  printf( "|ulp| <= 4.0: %10llu in [ %23.15e, %23.15e ] U [ %23.15e, %23.15e ]\n", ( m_hist[  4 ] ), get_negpos_ranges(  4 ) ); \
    m_hist[  5 ]    &&  printf( "|ulp| <= 5.0: %10llu in [ %23.15e, %23.15e ] U [ %23.15e, %23.15e ]\n", ( m_hist[  5 ] ), get_negpos_ranges(  5 ) ); \
    m_hist[  6 ]    &&  printf( "|ulp| <= 6.0: %10llu in [ %23.15e, %23.15e ] U [ %23.15e, %23.15e ]\n", ( m_hist[  6 ] ), get_negpos_ranges(  6 ) ); \
    m_hist[  7 ]    &&  printf( "|ulp| <= 7.0: %10llu in [ %23.15e, %23.15e ] U [ %23.15e, %23.15e ]\n", ( m_hist[  7 ] ), get_negpos_ranges(  7 ) ); \
    m_hist[  8 ]    &&  printf( "|ulp| <= 8.0: %10llu in [ %23.15e, %23.15e ] U [ %23.15e, %23.15e ]\n", ( m_hist[  8 ] ), get_negpos_ranges(  8 ) ); \
    m_hist[  9 ]    &&  printf( "|ulp| <= 9.0: %10llu in [ %23.15e, %23.15e ] U [ %23.15e, %23.15e ]\n", ( m_hist[  9 ] ), get_negpos_ranges(  9 ) ); \
    m_hist[ 10 ]    &&  printf( "|ulp| <= 10.: %10llu in [ %23.15e, %23.15e ] U [ %23.15e, %23.15e ]\n", ( m_hist[ 10 ] ), get_negpos_ranges( 10 ) ); \
    m_hist[ 11 ]    &&  printf( "|ulp| <= 11.: %10llu in [ %23.15e, %23.15e ] U [ %23.15e, %23.15e ]\n", ( m_hist[ 11 ] ), get_negpos_ranges( 11 ) ); \
    m_hist[ 12 ]    &&  printf( "|ulp| <= 12.: %10llu in [ %23.15e, %23.15e ] U [ %23.15e, %23.15e ]\n", ( m_hist[ 12 ] ), get_negpos_ranges( 12 ) ); \
    m_hist[ 13 ]    &&  printf( "|ulp| <= 13.: %10llu in [ %23.15e, %23.15e ] U [ %23.15e, %23.15e ]\n", ( m_hist[ 13 ] ), get_negpos_ranges( 13 ) ); \
    m_hist[ 14 ]    &&  printf( "|ulp| <= 14.: %10llu in [ %23.15e, %23.15e ] U [ %23.15e, %23.15e ]\n", ( m_hist[ 14 ] ), get_negpos_ranges( 14 ) ); \
    m_hist[ 15 ]    &&  printf( "|ulp| <= 15.: %10llu in [ %23.15e, %23.15e ] U [ %23.15e, %23.15e ]\n", ( m_hist[ 15 ] ), get_negpos_ranges( 15 ) ); \
    m_hist[ 16 ]    &&  printf( "|ulp| <= 16.: %10llu in [ %23.15e, %23.15e ] U [ %23.15e, %23.15e ]\n", ( m_hist[ 16 ] ), get_negpos_ranges( 16 ) ); \
    m_hist[ 17 ]    &&  printf( "|ulp| <= 17.: %10llu in [ %23.15e, %23.15e ] U [ %23.15e, %23.15e ]\n", ( m_hist[ 17 ] ), get_negpos_ranges( 17 ) ); \
    m_hist[ 18 ]    &&  printf( "|ulp| <= 18.: %10llu in [ %23.15e, %23.15e ] U [ %23.15e, %23.15e ]\n", ( m_hist[ 18 ] ), get_negpos_ranges( 18 ) ); \
    m_hist[ 19 ]    &&  printf( "|ulp| <= inf: %10llu in [ %23.15e, %23.15e ] U [ %23.15e, %23.15e ]\n", ( m_hist[ 19 ] ), get_negpos_ranges( 19 ) );

/* << ULPs HISTOGRAM */

/* >> RANGE FINDER */
/*
 * range [xmin, xmax] assumes xmin <= xmax, this is not validated inside the function
 */

#define MAX_NREL 10

static f64_t 	m_relf_range_left  [ MAX_NREL + 1 ];
static f64_t 	m_relf_range_right [ MAX_NREL + 1 ];
static u32_t	m_relf_range_len;
static fpset_t	m_relf_y_lhs_prev_type;

/*
* Update Range:
* 1. close current range at (m_len) with m_right
* 2. open next range: m_len++
* 3. update type and left value
*/

#define M_CLEAR_PREV_RANGES \
do { \
	for(u32_t i = 0; i < MAX_NREL; i++) \
	{ \
		m_relf_range_left[i]=0; \
		m_relf_range_right[i]=0; \
	} \
	m_relf_range_len = 0; \
} while(0)

#define INIT_RANGE(xl, y_lhs_type) \
do { \
	M_CLEAR_PREV_RANGES; \
	m_relf_y_lhs_prev_type = ( y_lhs_type ); \
	m_relf_range_left[ 0 ] = (xl); /* Open left of first range */ \
} while(0)

#define UPDATE_RANGE(xl, xr, y_lhs, y_lhs_type) \
do{ \
	if( m_relf_y_lhs_prev_type != y_lhs_type ) \
	{ \
		m_relf_range_right[ m_relf_range_len ] = (xl);			/* Close last range */ \
		m_relf_range_len+=(m_relf_range_len<(MAX_NREL-1)?1:0);	/* Count last range closed */ \
		m_relf_range_left[ m_relf_range_len ] = (xr);			/* Open a new range */ \
		m_relf_y_lhs_prev_type = y_lhs_type; \
	} \
} while(0)

#define CLOSE_RANGE(xr) \
do { \
	m_relf_range_right[ m_relf_range_len ] = (xr);			/* Close right of last range */ \
	m_relf_range_len+=m_relf_range_len<(MAX_NREL-1)?1:0;	/* Count last range closed */ \
} while(0)
/* << RANGE FINDER */

/* >> PROGRESS BAR */
static void console_erase_last_n_chars( u8_t n, fpenv_t fpenv )
{
    if( fpenv.log_completation_remove_history )
    {
        printf("\033[%dD", n);
        printf("\033[K");
        fflush(stdout);
    }

    else
    {
        printf("\n");
    }
}

static u32_t fp32_progress_bar_status( f32_t x )
{
    return ( (u32_t) ( 100.0 * ( x < 0.0 ? (39. - (x == 0 ? 0.0 : log10((f64_t) fabs( x )))): (x == 0.0 ? 83. : 83. + 44. + log10((f64_t) fabs( x )))) / 165.0 ) );
}

static void fp32_progress_bar( fpenv_t fpenv, f32_t x_curr, f32_t x_min, f32_t x_max, fp_progress_state_t * pstate )
{
    if( fpenv.log_completation )
    {
		fptriplet_t xtrip = f32_get_triplet( x_curr );

		pstate->count_points_in_the_range++;

        if( xtrip.e != (pstate->xexp) )
        {
            (pstate->xexp) = xtrip.e;

            if( !fp32_equals(x_curr, x_min) && !fp32_equals(x_curr, x_max) )
            {
                console_erase_last_n_chars( 250u, g_fpenv );

                const u32_t progress = fp32_progress_bar_status ( x_curr );

                printf("completed: %3u%% range-bound: %+dE2%+04d in-range-vals: %8llu",	progress,
																				xtrip.s,
                                                                   				pstate->xexp,
                                                                   				pstate->count_points_in_the_range );

                pstate->count_points_in_the_range = 0;
            }
        }

        else
        {
			/* NOP: already counted */
        }
    }
}

static u32_t fp64_progress_bar_status( f64_t x )
{
//    return ( (u32_t) ( 100.0 * ( x < 0.0 ? (39. - (x == 0 ? 0.0 : log10((f64_t) fabs( x )))): (x == 0.0 ? 83. : 83. + 44. + log10((f64_t) fabs( x )))) / 165.0 ) );
    return ( 0 ); /* TODO: adjust progress model */
}

static void fp64_progress_bar( fpenv_t fpenv, f64_t x_curr, f64_t x_min, f64_t x_max, fp_progress_state_t * pstate )
{
    if( fpenv.log_completation )
    {
		fptriplet_t xtrip = f64_get_triplet( x_curr );

		pstate->count_points_in_the_range++;

        if( xtrip.e != (pstate->xexp) )
        {
            (pstate->xexp) = xtrip.e;

            if( !fp64_equals(x_curr, x_min) && !fp64_equals(x_curr, x_max) )
            {
                console_erase_last_n_chars( 250u, g_fpenv );

                const u32_t progress = fp64_progress_bar_status ( x_curr );

                printf("completed: %3u%% range-bound: %+dE2%+05d in-range-vals: %8llu",	progress,
																				xtrip.s,
                                                                   				pstate->xexp,
                                                                   				pstate->count_points_in_the_range );

                pstate->count_points_in_the_range = 0;
            }
        }

        else
        {
			/* NOP: already counted */
        }
    }
}
/* << PROGRESS BAR */

/*
 * f32_ulps
 *
 * compute the number of ULPs between approx and exact.
 *
 * approx and exact are two machine numbers hence this function can't
 * give a resolution less than 1ulp. This is not a problem for the
 * purpose of this function as it is not designed to evaluate
 * rounding errors required by IEEE-754.
 */

f64_t ADDCALL fp32_ulps( f32_t approx, f32_t exact )
{
    f64_t ret;

    fpset_t t_apprx = f32_get_subset( approx );
    fpset_t t_exact = f32_get_subset( exact  );

    if( t_apprx != t_exact )
    {
        if( ( t_apprx == FP_TYPE_NORMAL || t_apprx == FP_TYPE_SUBNORMAL ) && ( t_exact == FP_TYPE_NORMAL || t_exact == FP_TYPE_SUBNORMAL ) )
        {
            ret = f32_ulps_impl( approx, ( (f80_t) approx - (f80_t) exact ) );
        }

        else
        {
            ret = m_inff64.f;
        }
    }

    else if( t_apprx == FP_TYPE_NAN )
    {
        ret = 0.0; /* Controversial: if both are NaN's ULP is zero TODO: research it */
    }

    else if( t_apprx == FP_TYPE_INFINITE )
    {
        if( approx == exact ) /* Same infinite sign */
        {
            ret = 0.0;
        }
        else
        {
            ret = m_inff64.f;
        }
    }

    else if( t_apprx == FP_TYPE_ZERO )
    {
        fw_t u_apprx = { .f = approx };
        fw_t u_exact = { .f = exact  };

        if( u_apprx.u == u_exact.u ) /* is same sign */
        {
            ret = 0.0;
        }

        else
        {
            ret = m_inff64.f;
        }
    }

    else
    {
        ret = f32_ulps_impl( approx, ( (f80_t) approx - (f80_t) exact ) );
    }

    return ( ret );
}

f64_t ADDCALL fp64_ulps( f64_t approx, f64_t exact )
{
    f64_t ret;

    fpset_t t_apprx = f64_get_subset( approx );
    fpset_t t_exact = f64_get_subset( exact  );

    if( t_apprx != t_exact )
    {
        if( ( t_apprx == FP_TYPE_NORMAL || t_apprx == FP_TYPE_SUBNORMAL ) && ( t_exact == FP_TYPE_NORMAL || t_exact == FP_TYPE_SUBNORMAL ) )
        {
            ret = f64_ulps_impl( approx, ( (f80_t) approx - (f80_t) exact ) );
        }

        else
        {
            ret = m_inff64.f;
        }
    }

    else if( t_apprx == FP_TYPE_NAN )
    {
        ret = 0.0; /* Controversial: if both are NaN's ULP is zero TODO: research it, alternative return ( m_qnanf64.f ) */
    }

    else if( t_apprx == FP_TYPE_INFINITE )
    {
        if( approx > 0.0 && exact > 0.0 )
        {
            ret = 0.0;
        }
        else
        {
            ret = m_inff64.f;
        }
    }

    else if( t_apprx == FP_TYPE_ZERO )
    {
        dw_t u_apprx = { .f = approx };
        dw_t u_exact = { .f = exact  };

        if( u_apprx.u == u_exact.u ) /* is same sign */
        {
            ret = 0.0;
        }

        else
        {
            ret = m_inff64.f;
        }
    }

    else
    {
        ret = f64_ulps_impl( approx, ( (f80_t) approx - (f80_t) exact ) );
    }

    return ( ret );
}

static f64_t f64_ulps_impl( f64_t approx, f80_t diff )
{
    dw_t bits = {.f = approx };

    const i16_t e   = ( ( bits.w[ DW0 ] & DW0_EMASK ) >> DOFFS_F64 ) - EMAX_F64;
    const f80_t ulp = ldexpl( 1. , e - ( NDIG_F64 - 1 ) );

    return ( (f64_t) ( diff / ulp ) );
}

static f64_t f32_ulps_impl( f32_t approx, f80_t diff )
{
    fw_t bits = {.f = approx };

    const i16_t e   = ( ( bits.w[ FW0 ] & FW0_EMASK ) >> DOFFS_F32 ) - EMAX_F32;
    const f80_t ulp = ldexpl( 1. , e - ( NDIG_F32 - 1 ) );

    return ( (f64_t) ( diff / ulp ) );
}

f32_t ADDCALL f32_get_named_fp_in_real_line( named_fp_t point )
{
    fw_t ret;

    switch(point)
    {
        case (NAMED_FP_INF )        : { ret.u = m_inff32 .u; break; }
        case (NAMED_FP_MAXNORM )    : { ret.u = m_maxf32 .u; break; }
        case (NAMED_FP_ZERO)        : { ret.u = m_zerof32.u; break; }
        case (NAMED_FP_MINNORM)     : { ret.u = m_minf32.u; break; }
        case (NAMED_FP_MINSUBN)     : { ret.u = m_subnf32.u; break; }
        default                     : { ret.u = m_qnanf32.u; break; }
    }

    return ( ret.f );
}

bool_t ADDCALL fp32_equals( f32_t a, f32_t b )
{
    const fw_t ua = { .f = a };
    const fw_t ub = { .f = b };

    return ( ua.u == ub.u );
}

bool_t ADDCALL fp64_equals( f64_t a, f64_t b )
{
    const dw_t ua = { .f = a };
    const dw_t ub = { .f = b };

    return ( ua.u == ub.u );
}

bool_t ADDCALL fp32_equals_sign( f32_t a, f32_t b )
{
    const fw_t ua = { .f = a };
    const fw_t ub = { .f = b };

    return ( ( ua.w[ FW0 ] & FW0_SMASK ) == ( ub.w[ FW0 ] & FW0_SMASK ) );
}

bool_t ADDCALL fp64_equals_sign( f64_t a, f64_t b )
{
    const dw_t ua = { .f = a };
    const dw_t ub = { .f = b };

    return ( ( ua.w[ FW0 ] & FW0_SMASK ) == ( ub.w[ FW0 ] & FW0_SMASK ) );
}

i16_t ADDCALL fp32_get_exp( f32_t x )
{
    return ( F32_GET_EXP( x ) );
}

i16_t ADDCALL fp64_get_exp( f64_t x )
{
    return ( F64_GET_EXP( x ) );
}

cstr_t ADDCALL fp32_get_subset_name( f32_t x )
{
    cstr_t ret;

    switch ( f32_get_subset( x ) )
    {
        case (FP_TYPE_NAN       )   : { ret = "FP_TYPE_NAN      "; break; }
        case (FP_TYPE_ZERO      )   : { ret = "FP_TYPE_ZERO     "; break; }
        case (FP_TYPE_NORMAL    )   : { ret = "FP_TYPE_NORMAL   "; break; }
        case (FP_TYPE_SUBNORMAL )   : { ret = "FP_TYPE_SUBNORMAL"; break; }
        case (FP_TYPE_INFINITE  )   : { ret = "FP_TYPE_INFINITE "; break; }
        default                     : { ret = "UNKNOWN          "; break; }
    }

    return ( ret );
}

cstr_t ADDCALL fp64_get_subset_name( f64_t x )
{
    cstr_t ret;

    switch ( f64_get_subset( x ) )
    {
        case (FP_TYPE_NAN       )   : { ret = "FP_TYPE_NAN      "; break; }
        case (FP_TYPE_ZERO      )   : { ret = "FP_TYPE_ZERO     "; break; }
        case (FP_TYPE_NORMAL    )   : { ret = "FP_TYPE_NORMAL   "; break; }
        case (FP_TYPE_SUBNORMAL )   : { ret = "FP_TYPE_SUBNORMAL"; break; }
        case (FP_TYPE_INFINITE  )   : { ret = "FP_TYPE_INFINITE "; break; }
        default                     : { ret = "UNKNOWN          "; break; }
    }

    return ( ret );
}

f64_t ADDCALL f64_get_named_fp_in_real_line( named_fp_t point )
{
    dw_t ret;

    switch(point)
    {
        case (NAMED_FP_INF )        : { ret.u = m_inff64 .u; break; }
        case (NAMED_FP_MAXNORM )    : { ret.u = m_maxf64 .u; break; }
        case (NAMED_FP_ZERO)        : { ret.u = m_zerof64.u; break; }
        case (NAMED_FP_MINNORM)     : { ret.u = m_minf64 .u; break; }
        case (NAMED_FP_MINSUBN)     : { ret.u = m_subnf64.u; break; }
        default                     : { ret.u = m_qnanf64.u; break; }
    }

    return ( ret.f );
}

/*
 * get_subset
 *
 * recovers the subset belonging to 'x' according to IEEE 754-2008 standard
 * the algorithm is based on section 3.4 Binary interchange format encodings of cited standard
 *
 */

fpset_t ADDCALL f32_get_subset( f32_t x )
{
    fpset_t ret;

    const fw_t y = { .f = x };

    const u16_t e  = y.w[ FW0 ] & FW0_EMASK;
    const u32_t t  = y.u & 0x7FFFFFu;

    if ( e == FW0_EMASK )
    {
        /* cond. a) E = 2**w - 1 and T != 0 */
        if( t != 0 )
        {
            ret = FP_TYPE_NAN;
        }

        /* cond. b) E = 2**w - 1 and T = 0 */
        else
        {
            ret = FP_TYPE_INFINITE;
        }
    }

    /* cond. c) 0 < E < 2**w - 1 */
    else if ( e != 0 )
    {
        ret = FP_TYPE_NORMAL;
    }

    /* cond. d) E = 0 and T != 0 */
    else if ( t != 0 )
    {
        ret = FP_TYPE_SUBNORMAL;
    }

    /* cond. e) E = 0 and T = 0 */
    else
    {
        ret = FP_TYPE_ZERO;
    }

    return ( ret );
}

fpset_t ADDCALL f64_get_subset( f64_t x )
{
    fpset_t ret;

    const dw_t y = { .f = x };

    const u16_t e  = y.w[ DW0 ] & DW0_EMASK;
    const u64_t t  = y.u & 0x000FFFFFFFFFFFFFuLL;

    if ( e == DW0_EMASK )
    {
        /* cond. a) E = 2**w - 1 and T != 0 */
        if( t != 0 )
        {
            ret = FP_TYPE_NAN;
        }

        /* cond. b) E = 2**w - 1 and T = 0 */
        else
        {
            ret = FP_TYPE_INFINITE;
        }
    }

    /* cond. c) 0 < E < 2**w - 1 */
    else if ( e != 0 )
    {
        ret = FP_TYPE_NORMAL;
    }

    /* cond. d) E = 0 and T != 0 */
    else if ( t != 0 )
    {
        ret = FP_TYPE_SUBNORMAL;
    }

    /* cond. e) E = 0 and T = 0 */
    else
    {
        ret = FP_TYPE_ZERO;
    }

    return ( ret );
}

/*
 * fp64_step_real_line
 *
 * after n calls:
 *
 * x_0 in [ -MAX, +MAX ]
 *
 * x_n = x_n-1 · s**n,
 *
 * with x_n in ( -MAX, -NORM ] U [ +NORM, +MAX ] U { +INF }
 *
 * frac comes from: ratio = 1 + frac, example: ratio = 1.0000001, then frac = 0.0000001
 *
 */

f32_t ADDCALL f32_geom_step_real_line( f32_t x, f32_t frac )
{
    f32_t ret;

    const f32_t right_zero = +f32_get_named_fp_in_real_line( NAMED_FP_MINNORM );
    const f32_t left_zero  = -right_zero;

    if( x < left_zero ) /* [-MAX, -NORM) -> (-MAX, -NORM] U [-NORM, -SUBN] U { -ZERO } */
    {
        ret = x / ( 1. + frac );
    }

    /* else if( x < left_sub )  { TBD } */
    /* else if( x == 0.0 )      { TBD } */
    /* else if( x < right_sub ) { TBD } */

    else if( x < right_zero ) /* x in [-NORM, +NORM] -> +NORM */
    {
        ret = right_zero;
    }

    else
    {
        ret = x * ( 1. + frac );
    }

    return ( ret );
}

f64_t ADDCALL f64_geom_step_real_line( f64_t x, f64_t frac )
{
    f64_t ret;

    const f64_t right_zero = +f64_get_named_fp_in_real_line( NAMED_FP_MINNORM );
    const f64_t left_zero  = -right_zero;

    if( x < left_zero ) /* [-MAX, -NORM) -> (-MAX, -NORM] U [-NORM, -SUBN] U { -ZERO } */
    {
        ret = x / ( 1. + frac ) ;
    }

    /* else if( x < left_sub )  { TBD } */
    /* else if( x == 0.0 )      { TBD } */
    /* else if( x < right_sub ) { TBD } */

    else if( x < right_zero ) /* x in [-NORM, +NORM] -> +NORM */
    {
        ret = right_zero;
    }

    else
    {
        ret = x * ( 1. + frac );
    }

    return ( ret );
}

/**
 * rel_error_for_reals
 *
 * precondition: x, y shall lie in the Real Line
 *
 * ret in [0, 1]
 */

f64_t ADDCALL f64_rel_error_for_reals( f64_t x, f64_t y )
{
    f64_t ret;

    fpset_t tx = f64_get_subset( x );
    fpset_t ty = f64_get_subset( y );

    /* tx or ty inf */
    if( ( tx == FP_TYPE_INFINITE ) || ( ty == FP_TYPE_INFINITE ) )
    {
        if ( x == y ) /* same sign inf */
        {
            ret = 52; /* avoid performing: inf - inf wich evaluates to nan */
        }

        else
        {
            ret = 0;
        }
    }

    else if( ( tx == FP_TYPE_ZERO ) || ( ty == FP_TYPE_ZERO ) )
    {
        const dw_t wx = { .f = x };
        const dw_t wy = { .f = y };

        if ( x == y )
        {
            if( ( wx.w[ FW0 ] & FW0_SMASK ) == ( wy.w[ FW0 ] & FW0_SMASK ) )
            {
                ret = 52;
            }

            else
            {
                ret = 52;
            }
        }

        else
        {
            ret = 0;
        }
    }

    else if( ( tx == FP_TYPE_SUBNORMAL ) || ( ty == FP_TYPE_SUBNORMAL ) )
    {
        ret = 0; /* TODO: renormalize, count the difference between exponent values */
    }

    else
    {
        if( x == y )
        {
            ret = 52;
        }

        else
        {
            ret = fabs( y / x ) - 1.0;
        }
    }

    return ( ret );
}

cstr_t ADDCALL f32_sprint_digits_radix2( char buff [ 32 + 1 + 2 ], char_t sepparator, f32_t x )
{
    static const u16_t index[ 2 ] = { FW0, FW1, };

    fw_t f = { .f = x };

    u16_t w;
    u8_t i = 0;
    u8_t j = 0;
    u8_t k = 0;

    while( k < 2 )
    {
        w = f.w[ index[ k ] ];

        while( j < 16 )
        {
            buff[ i++ ] = w & 0x8000 ? '1' : '0';

            if( ( i == 1 ) || ( i == 10 ) )
            {
                buff[ i++ ] = sepparator;
            }

            w = w << 1;
            j++;
        }

        j = 0;
        k++;
    }

    buff[ i ] = '\0';

    return ( buff );
}

cstr_t ADDCALL f64_sprint_digits_radix2( char buff [ 64 + 1 + 2 ], char_t sepparator, f64_t x )
{
    static const u16_t index[ 4 ] = { DW0, DW1, DW2, DW3, };

    dw_t f = { .f = x };

    u16_t w;
    u8_t i = 0;
    u8_t j = 0;
    u8_t k = 0;

    while( k < 2 )
    {
        w = f.w[ index[ k ] ];

        while( j < 16 )
        {
            buff[ i++ ] = w & 0x8000 ? '1' : '0';

            if( ( i == 1 ) || ( i == 13 ) )
            {
                buff[ i++ ] = sepparator;
            }

            w = w << 1;
            j++;
        }

        j = 0;
        k++;
    }

    buff[ i ] = '\0';

    return ( buff );
}

/**
 * f32_mount_bitfields
 *
 * s: sign bit
 * e: biased exponent bits
 * m: fractional bits
 *
 * return: 32bit floating-point representation according to IEEE-754
 */

f32_t ADDCALL f32_mount_bitfields( u32_t s, u32_t e, u32_t m )
{
    return ( (fw_t) { .u = ( s << 31 ) | ( e << 23 ) | ( m ) } ).f;
}

/**
 * f64_mount_bitfields
 *
 * s: sign bit
 * e: biased exponent bits
 * m: fractional bits
 *
 * return: 64bit floating-point representation according to IEEE-754
 */

f64_t ADDCALL f64_mount_bitfields( u64_t s, u64_t e, u64_t m )
{
    return ( (dw_t) { .u = ( s << 63 ) | ( e << 52 ) | ( m ) }).f;
}

f32_t ADDCALL f32_next_float( f32_t x )
{
    f32_t ret;

    fpset_t type = f32_get_subset( x );

    switch (type)
    {
        case (FP_TYPE_NAN):
        {
            ret = x;
            break;
        }

        case(FP_TYPE_INFINITE):
        {
            if( x < 0.0 )
            {
                ret = -m_maxf32.f;
            }

            else
            {
                ret = m_inff32.f;
            }

            break;
        }

        case(FP_TYPE_ZERO):
        case(FP_TYPE_NORMAL):
        case(FP_TYPE_SUBNORMAL):
        default:
        {
            fw_t this = { .f = x };

            u32_t s = ( this.w[ FW0 ] & 0x8000 ) >> 15;

            i32_t e = ( this.w[ FW0 ] & 0x7F80 ) >> 7;

            u32_t m = ( ( (u32_t) (this.w[ FW0 ] & 0x007F) ) << 16 ) | ( (u32_t) (  this.w[ FW1 ] & 0xFFFF ) );

            if( s == 1 && e == 0 && m == 0 )
            {
                s = 0;
            }

            else if( s == 1 ) /* when x is negative, to move x->x + dx requires to substract values from e, m */
            {
                if( m == 0x000000 ) { e -= 1; m = 0x7FFFFF; } else { m -= 1; }
            }

            else
            {
                if( m >= 0x7FFFFF ) { e += 1; m = 0x000000; } else { m += 1; }
            }

            ret = f32_mount_bitfields( s, e, m );
        }
    }

    return ( ret );
}

f64_t ADDCALL f64_next_float( f64_t x )
{
    f64_t ret;

    fpset_t type = f64_get_subset( x );

    switch (type)
    {
        case (FP_TYPE_NAN):
        {
            ret = x;
            break;
        }

        case(FP_TYPE_INFINITE):
        {
            if( x < 0.0 )
            {
                ret = -m_maxf64.f;
            }

            else
            {
                ret = m_inff64.f;
            }

            break;
        }

        case(FP_TYPE_ZERO):
        case(FP_TYPE_NORMAL):
        case(FP_TYPE_SUBNORMAL):
        default:
        {
            dw_t this = { .f = x };

            u64_t s = ( this.w[ DW0 ] & 0x8000 ) >> 15;
            u64_t e = ( this.w[ DW0 ] & 0x7FF0  ) >> 4;
            u64_t m =   this.u        & 0x000FFFFFFFFFFFFFllu;

            if( s == 1 && e == 0 && m == 0 )
            {
                s = 0;
            }

            else if( s == 1 ) /* when x is negative, to move x->x + dx requires to substract values from e, m */
            {
                if( m == 0x000000 ) { e -= 1; m = 0x000FFFFFFFFFFFFF; } else { m -= 1; }
            }

            else
            {
                if( m >= 0x000FFFFFFFFFFFFF ) { e += 1; m = 0; } else { m += 1; }
            }

            ret = f64_mount_bitfields( s, e, m );
        }
    }

    return ( ret );
}

/*
 * f32_get_triplet
 *
 * The most interesting difference between IEEE-754 enconding is at how NORMAL and SUBNORMAL numbers are identified.
 * Both range of numbers share one value of the exponent range which is emin. In fact all subnormals have e = emin.
 * According to standard, normals has set to 1 the MSB of the fractional part, so this feature is kept here.
 *
 * An alternative approach, not yet considered is to renormalize subnormal numbers moving the fractional part up,
 * which corresponds with a scaling equal to the value of the radix. The scaling can be substracted from the exponent,
 * allowing subnormals to reach values below emin.
 */

fptriplet_t ADDCALL f32_get_triplet( f32_t x )
{
    fptriplet_t ret;

    fpset_t type = f32_get_subset( x );

    fw_t y = { .f = x };

    if( ( y.w[ FW0 ] & FW0_SMASK ) == 0 )
    {
        ret.s = +1;
    }

    else
    {
        ret.s = -1;
    }

    switch (type)
    {
        case (FP_TYPE_NAN):
        {
            ret.e = 1 + ( EMAX_F32 << 1 );
            ret.M = 1;
            break;
        }

        case (FP_TYPE_INFINITE):
        {
            ret.e = 1 + ( EMAX_F32 << 1 );
            ret.M = 0;
            break;
        }

        case (FP_TYPE_ZERO):
        {
            ret.e = 0;
            ret.M = 0;
            break;
        }

        case (FP_TYPE_SUBNORMAL): /* Check that bit 24 is set to 0 */
        {
            ret.e = 1 - EMAX_F32;
            ret.M = ( ( (u32_t) ( y.w[ FW0 ] & FW0_DMASK ) ) << 16 ) | ( (u32_t) ( y.w[ FW1 ] ) );
            break;
        }

        case (FP_TYPE_NORMAL): /* Check that bit 24 is set to 1 */
        default:
        {
            ret.e =   ( (i16_t) ( ( y.w[ FW0 ] & FW0_EMASK ) >> DOFFS_F32 ) ) - EMAX_F32;
            ret.M = ( 1 << 23 ) | ( ( (u32_t) ( y.w[ FW0 ] & FW0_DMASK ) ) << 16 ) | ( (u32_t) ( y.w[ FW1 ] ) );
            break;
        }
    }

    return ( ret );
}

fptriplet_t ADDCALL f64_get_triplet( f64_t x )
{
    fptriplet_t ret;

    fpset_t type = f64_get_subset( x );

    dw_t y = { .f = x };

    if( ( y.w[ DW0 ] & DW0_SMASK ) == 0 )
    {
        ret.s = +1;
    }

    else
    {
        ret.s = -1;
    }

    switch (type)
    {
        case (FP_TYPE_NAN):
        {
            ret.e = 1 + ( EMAX_F64 << 1 );
            ret.M = 1;
            break;
        }

        case (FP_TYPE_INFINITE):
        {
            ret.e = 1 + ( EMAX_F64 << 1 );
            ret.M = 0;
            break;
        }

        case (FP_TYPE_ZERO):
        {
            ret.e = 0;
            ret.M = 0;
            break;
        }

        case (FP_TYPE_SUBNORMAL): /* Check that bit 24 is set to 0 */
        {
            ret.e = 1 - EMAX_F64;
            ret.M = ( ( (u64_t) ( y.w[ DW0 ] & DW0_DMASK ) ) << 16 ) | ( (u64_t) ( y.w[ DW1 ] ) );
            break;
        }

        case (FP_TYPE_NORMAL): /* Check that bit 24 is set to 1 */
        default:
        {
            ret.e =   ( (i16_t) ( ( y.w[ DW0 ] & DW0_EMASK ) >> DOFFS_F64 ) ) - EMAX_F64;

            ret.M = ( 1ull << 52 ) |
                    ( ( (u64_t) ( y.w[ DW0 ] & DW0_DMASK ) ) << ( 16 + 32 ) ) |
                    ( (u64_t) ( y.w[ DW1 ] ) << 32 ) |
                    ( (u64_t) ( y.w[ DW2 ] ) << 16 ) |
                    ( (u64_t) ( y.w[ DW3 ] ) );

            break;
        }
    }

    return ( ret );
}

/**
 * f32_eval_triplet
 *
 * given a floating point triplet x(s,M,e), evaluates the floating-point representation according to the format IEEE-754
 *
 * return the value of evaluate (-1)^s · M · B^{e-p+1}
 */

f32_t ADDCALL f32_eval_triplet( fptriplet_t x )
{
    fw_t ret;

    if( x.e > EMAX_F32 )
    {
        if( x.M != 0 )
        {
            ret.f = m_qnanf32.f;
        }

        else
        {
            ret.f = ( (f32_t) x.s ) * m_inff32.f;
        }
    }

    else
    {
        ret.f = f32_set_exp( 1., x.e );

        if( x.s == -1 )
        {
            ret.w[ FW0 ] = ret.w[ FW0 ] | FW0_SMASK;
        }

        ret.u = ret.u | ( (u32_t) x.M );
    }

    return ( ret.f );
}

/**
 * f64_eval_triplet
 *
 * given a floating point triplet x(s,M,e), evaluates the floating-point representation according to the format IEEE-754
 *
 * return the value of evaluate (-1)^s · M · B^{e-p+1}
 */

f64_t ADDCALL f64_eval_triplet( fptriplet_t x )
{
    dw_t ret;

    if( x.e > EMAX_F64 )
    {
        if( x.M != 0 )
        {
            ret.f = m_qnanf64.f;
        }

        else
        {
            ret.f = ( (f64_t) x.s ) * m_inff64.f;
        }
    }

    else
    {
        ret.f = f64_set_exp( 1., x.e );

        if( x.s == -1 )
        {
            ret.w[ DW0 ] = ret.w[ DW0 ] | DW0_SMASK;
        }

        ret.u = ret.u | x.M;
    }

    return ( ret.f );
}

/**
 * f32_set_exp
 *
 * x: floating-point s, m
 *                                              s       e
 * e: new (unbiased) exponent in the format (-1) · m · B
 *
 * return: a floating-point representation compliant with IEEE-754 | S(1bit) | E(8bit) | M(23bit) |
 *
 * NOTE: if x = 1.0, set_exp( x, e ) is equivalent to y = 2^e, with e in emin <= e <= emax
 */

f32_t ADDCALL f32_set_exp ( f32_t x, i16_t e )
{
	fw_t ret;

    switch ( f32_get_subset( x ) )
    {
        case (FP_TYPE_NORMAL):
        {
            if( e < EMIN_F32 )
            {
                /* TODO: add support to handle return subnormals when EMIN_SUBNORMAL <= e < EMIN_F64 */
                ret.f = 0.0;
            }

            else if( e > EMAX_F32 )
            {
                ret.f = m_inff32.f;
            }

            else
            {
                ret.f = x;
                F32_SET_EXP( ret.w[ FW0 ], e);
            }

            break;
        }

        default:
        case (FP_TYPE_NAN):
        case (FP_TYPE_INFINITE):
        case (FP_TYPE_ZERO):
        case (FP_TYPE_SUBNORMAL):
        {
            ret.f = x;
            break;
        }
    }

    return ret.f;
}

f64_t ADDCALL f64_set_exp ( f64_t x, i16_t e )
{
	dw_t ret;

    switch ( f64_get_subset( x ) )
    {
        case (FP_TYPE_NORMAL):
        {
            if( e < EMIN_F64 )
            {
                /* TODO: add support to handle return subnormals when EMIN_SUBNORMAL <= e < EMIN_F64 */
                ret.f = 0.0;
            }

            else if( e > EMAX_F64 )
            {
                ret.f = m_inff64.f;
            }

            else
            {
                ret.f = x;
                F64_SET_EXP( ret.w[ DW0 ], e);
            }

            break;
        }

        default:
        case (FP_TYPE_NAN):
        case (FP_TYPE_INFINITE):
        case (FP_TYPE_ZERO):
        case (FP_TYPE_SUBNORMAL):
        {
            ret.f = x;
            break;
        }
    }

    return ret.f;
}

fp32_vec2_t ADDCALL fp32_find_control_boundaries( f32_t at_x, const f32_t * control_points,  f64_t boundary_semi_length_ulps, i32_t n )
{
    f32_t local_epsl;
    f32_t local_epsr;

    f32_t control_left  = -m_inff32.f;
    f32_t control_right = +m_inff32.f;

    /*
     * As a protection mechanism, in case of passing null-pointer in control_points,
     * It is considered there is one single range, with left and right points at -inf and +inf respectively
     */

    if( control_points == NULLPTR )
    {
        /* NOP: already set to inf */
    }

    else
    {
        for( i32_t i = 0; i < n; i++ )
        {
            const f32_t control_center = control_points[ i ];

            if( control_center == 1.0 || control_center == -1.0 )
            {
                local_epsl      = 0.5f * m_epsf32.f;
                local_epsr      =        m_epsf32.f;
            }

            else if( control_center == 0.0 )
            {
                control_left    = -0.0;
                control_right   = +0.0;

                local_epsl      = 0.0f;
                local_epsr      = 0.0f;
            }

            else /* general method */
            {
                local_epsl      = ( f32_next_float( control_center ) - control_center ) / control_center;
                local_epsr      = local_epsl;
            }

            control_left  = control_center * ( 1.f - boundary_semi_length_ulps * local_epsl );
            control_right = control_center * ( 1.f + boundary_semi_length_ulps * local_epsr );

            if ( control_right > at_x )
            {
                break;
            }
        }
    }

    return ( (fp32_vec2_t) { .x0 = control_left, .x1 = control_right, } );
}

fp64_vec2_t ADDCALL fp64_find_control_boundaries( f64_t at_x, const f64_t * control_points,  f64_t boundary_semi_length_ulps, i32_t n )
{
    f64_t local_epsl;
    f64_t local_epsr;

    f64_t control_left  = -m_inff64.f;
    f64_t control_right = +m_inff64.f;

    /*
     * As a protection mechanism, in case of passing null-pointer in control_points,
     * It is considered there is one single range, with left and right points at -inf and +inf respectively
     */

    if( control_points == NULLPTR )
    {
        /* NOP: already set to inf */
    }

    else
    {
        for( i32_t i = 0; i < n; i++ )
        {
            const f64_t control_center = control_points[ i ];

            if( control_center == 1.0 || control_center == -1.0 )
            {
                local_epsl      = 0.5 * m_epsf64.f;
                local_epsr      =       m_epsf64.f;
            }

            else if( control_center == 0.0 )
            {
                control_left    = -0.0;
                control_right   = +0.0;

                local_epsl      = 0.0;
                local_epsr      = 0.0;
            }

            else /* general method */
            {
                local_epsl      = ( f64_next_float( control_center ) - control_center ) / control_center;
                local_epsr      = local_epsl;
            }

            control_left  = control_center * ( 1. - boundary_semi_length_ulps * local_epsl );
            control_right = control_center * ( 1. + boundary_semi_length_ulps * local_epsr );

            if ( control_right > at_x )
            {
                break;
            }
        }
    }

    return ( (fp64_vec2_t) { .x0 = control_left, .x1 = control_right, } );
}

/*
 * fp32_next_x
 *
 * compute the next value of x, considering certain points of control,
 * in which the next value is computed step by step to ensure the point of control is not skipped.
 *
 * Cases description:
 *
 * - CASE 1: x left   control -> next left control            : geom grow
 * - CASE 2: x left   control -> next passes left control     : step grow
 * - CASE 3: x inside control                                 : step grow
 * - CASE 4: x inside control zero (child of CASE 3)          : step grow
 * - CASE 5: x is one of the two infinites                    : step grow (-inf), clamp (+inf)
 */

f32_t ADDCALL fp32_next_x( f32_t x, f32_t frac, const f32_t * control_points, f64_t round_ulps, i32_t n )
{
    const fp32_vec2_t   control_point = fp32_find_control_boundaries( x, control_points, round_ulps, n );
    const f32_t         control_left  = control_point.x0;
    const f32_t         control_right = control_point.x1;

    f32_t next = f32_geom_step_real_line( x, frac );

    /* CASE 2: x left   control -> next passes left control     : step grow */
    if( ( x < control_left ) && ( next >= control_left ) )
    {
        if( control_left == 0.0 )
        {
            next = -0.0; /* This helps to the user to avoid set the control point 0.0 as { -0.0 }*/
            g_fpenv.debug_next_x_inside_boundaries && printf("CASE 2: %+.1f, [%f -> %f], %+.1f\n", control_left, x, next, control_right);
        }

        else
        {
            next = control_left;
        }
    }

    /* CASE 3: x inside control                                 : step grow */
    else if( ( control_left <= x ) && ( x <= control_right ) )
    {
        next = f32_next_float( x );
        g_fpenv.debug_next_x_inside_boundaries && printf("CASE 3: %f, [%f -> %f], %f\n", control_left, x, next, control_right);
    }

    /* CASE 4: x inside control zero (child of CASE 3)          : step grow */
    else if ( f32_get_subset( x ) == FP_TYPE_ZERO )
    {
        next = f32_next_float( x );
        g_fpenv.debug_next_x_inside_boundaries && printf("CASE 4: %+.1f, [%f -> %f], %+.1f\n", control_left, x, next, control_right);
    }

    /* CASE 5: x is one of the two infinites                    : step grow (-inf), clamp (+inf) */
    else if ( f32_get_subset( x ) == FP_TYPE_INFINITE )
    {
        if( x == (-m_inff32.f) )
        {
            next = f32_next_float( x );
            g_fpenv.debug_next_x_inside_boundaries && printf("CASE 5: %+.1f, [%f -> %f], %+.1f\n", control_left, x, next, control_right);
        }

        else
        {
            next = x;
        }
    }

    /* CASE 1: x left   control -> next left control            : geom grow */
    else
    {
        /* NOP: geometrical grow already applied */
    }

    return ( next );
}

f64_t ADDCALL fp64_next_x( f64_t x, f64_t frac, const f64_t * control_points, f64_t round_ulps, i32_t n )
{
    const fp64_vec2_t   control_point = fp64_find_control_boundaries( x, control_points, round_ulps, n );
    const f64_t         control_left  = control_point.x0;
    const f64_t         control_right = control_point.x1;

    f64_t next = f64_geom_step_real_line( x, frac );

    /* CASE 2: x left   control -> next passes left control     : step grow */
    if( ( x < control_left ) && ( next >= control_left ) )
    {
        if( control_left == 0.0 )
        {
            next = -0.0; /* This helps to the user to avoid set the control point 0.0 as { -0.0 }*/
            g_fpenv.debug_next_x_inside_boundaries && printf("CASE 2: %+.1f, [%f -> %f], %+.1f\n", control_left, x, next, control_right);
        }

        else
        {
            next = control_left;
        }
    }

    /* CASE 3: x inside control                                 : step grow */
    else if( ( control_left <= x ) && ( x <= control_right ) )
    {
        next = f64_next_float( x );
        g_fpenv.debug_next_x_inside_boundaries && printf("CASE 3: %f, [%f -> %f], %f\n", control_left, x, next, control_right);
    }

    /* CASE 4: x inside control zero (child of CASE 3)          : step grow */
    else if ( f64_get_subset( x ) == FP_TYPE_ZERO )
    {
        next = f64_next_float( x );
        g_fpenv.debug_next_x_inside_boundaries && printf("CASE 4: %+.1f, [%f -> %f], %+.1f\n", control_left, x, next, control_right);
    }

    /* CASE 5: x is one of the two infinites                    : step grow (-inf), clamp (+inf) */
    else if ( f64_get_subset( x ) == FP_TYPE_INFINITE )
    {
        if( x == (-m_inff64.f) )
        {
            next = f64_next_float( x );
            g_fpenv.debug_next_x_inside_boundaries && printf("CASE 5: %+.1f, [%f -> %f], %+.1f\n", control_left, x, next, control_right);
        }

        else
        {
            next = x;
        }
    }

    /* CASE 1: x left   control -> next left control            : geom grow */
    else
    {
        /* NOP: geometrical grow already applied */
    }

    return ( next );
}

/*
 * test_range_analyzer ( desc, approx, expect, accept, reject )
 *
 * desc:    name of the function under test plus the test description, that will appear in the report
 * lhs:     function pointer to the left-hand-side value to compare with the rhs
 * rhs:     function pointer to the right-hand-side value to compare with the lhs
 * accept:  max difference in ulps between lhs and rhs to consider the test is passed
 * reject:  max difference in ulps between lhs and rhs to consider the test is reject
 *
 * state:   the state of the function when a test is rejected, useful to fast iterate a solution.
 *          If a not NULL state is given, the function will continue from that point, so in case of
 *          early rejection, the programmer can fix-up the code and continue the test from that, in order
 *          to know if the error is fixed or not. After this is done, the test can be reset by just
 *          passing a NULL pointer into state.
 *
 *          The state is passed as an ADT, but no ctor or dtor are given, so instead of a void *
 *          the programmer is requested to hold a minimum array size in bytes and pass its address.
 *          The state is expected to be a tiny structure, so there is no problem to hold it in the stack.
 */

void ADDCALL fp32_range_analyzer( cstr_t desc, f3232_t lhs, f3232_t rhs, f3232_t next_x, f64_t accept, f64_t reject, u8_t state[ /* TBD: Size in bytes of the state */ ] )
{
    printf("Analyzing %s\n"LINE_SEPARATOR"\n", desc);

    /* >> MAIN ANALYSIS */
    const f32_t x_min   = -f32_get_named_fp_in_real_line( NAMED_FP_INF );
    const f32_t x_max   = +f32_get_named_fp_in_real_line( NAMED_FP_INF );

    f32_t x_curr        = x_min;
    /* << MAIN ANALYSIS */

    /* >> SERVICE: PROGRESS BAR */
    fp_progress_state_t pstate =
    {
        .xexp = EMIN_F32,
        .count_points_in_the_range = 0,
    };
    /* << SERVICE: PROGRESS BAR */

    /* >> SERVICE: POINTS COUNTER */
    u64_t count_hi = 0;
    u64_t total_fp32 = (~0u   - (2u   * 8388607u  ) + 1u   ) ; /* 2^32 - all f32 NaNs */
    /* << SERVICE: POINTS COUNTER */

    /* >> SERVICE: ULP */
    ULPS_HIST_INIT;
    /* << SERVICE: ULP */

    /* >> SERVICE: RANGE FINDER */
    f32_t x_prev = x_curr;
    INIT_RANGE( x_prev, f32_get_subset( lhs ( x_prev ) ) );
    /* << SERVICE: RANGE FINDER */

    /* >> SERVICE: ERROR LOOP TERMINATION */
    i32_t x_curr_ntimes_analyzed = 0;
    /* << SERVICE: ERROR LOOP TERMINATION */

    /* >> MAIN ANALYSIS */
    while( ( x_curr <= x_max ) && ( x_curr_ntimes_analyzed == 0 ) )
    {
        /* >> MAIN ANALYSIS */
        const f32_t     y_lhs       = lhs ( x_curr );
        const f32_t     y_rhs       = rhs ( x_curr );
        const fpset_t   y_lhs_type  = f32_get_subset( y_lhs );
        /* << MAIN ANALYSIS */

        /* >> SERVICE: ULP */
        const f64_t curr_ulp = fabs( fp32_ulps( y_lhs, y_rhs ) );
        ULPS_HIST_RECORD( curr_ulp, x_curr );
        /* << SERVICE: ULP */

        {/* >> SERVICE: PROGRESS BAR */
            fp32_progress_bar( g_fpenv, x_curr, x_min, x_max, &pstate );
        }/* << SERVICE: PROGRESS BAR */

        {/* >> SERVICE: POINTS COUNTER */
            count_hi++;
        }/* << SERVICE: POINTS COUNTER */

        {/* >> SERVICE: RANGE FINDER */
            UPDATE_RANGE(x_prev, x_curr, y_lhs, y_lhs_type);
            x_prev = x_curr;
        }/* << SERVICE: RANGE FINDER */

        {/* >> MAIN ANALYSIS */
            x_curr = next_x( x_curr ); /* Update x */
        }/* >> MAIN ANALYSIS */

        {/* >> SERVICE: ERROR LOOP TERMINATION */
            x_curr_ntimes_analyzed = ( fp32_equals( x_curr, x_prev ) ? ( x_curr_ntimes_analyzed + 1 ) : 0 );
        }/* << SERVICE: ERROR LOOP TERMINATION */
    }
    /* << MAIN ANALYSIS */

    {/* >> SERVICE: PROGRESS BAR */
        console_erase_last_n_chars( 255u, g_fpenv );
    }/* << SERVICE: PROGRESS BAR */

    {/* >> SERVICE: ERROR LOOP TERMINATION */
        if( x_curr_ntimes_analyzed != 0 )
        {
            printf("[INFO]: loop stopped to prevent double analysis of point %f\n", x_curr);
        }
    }/* << SERVICE: ERROR LOOP TERMINATION */

    {/* >> SERVICE: RANGE FINDER */
    CLOSE_RANGE(x_prev);
    }/* << SERVICE: RANGE FINDER */

    /*
     * REPORTING SECTION
     * TODO: define the state of this function and send all of this to another function is possible
     */

    {/* >> SERVICE: RANGE FINDER */
        printf("\n Range Len:\n"LINE_SEPARATOR"\n");
        for(u32_t i = 0; i < m_relf_range_len; i++)
        {
            const f32_t yl = lhs( m_relf_range_left[i] );
            const f32_t yr = lhs( m_relf_range_right[i] );

            fpset_t yl_type = f32_get_subset( yl );
            cstr_t  yl_name = fp32_get_subset_name( yl );

            if( ( yl_type == FP_TYPE_NAN ) || ( yl_type == FP_TYPE_INFINITE ) || ( yl_type == FP_TYPE_ZERO ) )
            {
                printf("[ %+12.2e, %+12.2e ] -> %s\n", m_relf_range_left[i], m_relf_range_right[i], yl_name );
            }
            else
            {
                printf("[ %+12.2e, %+12.2e ] -> %s [ %+.2e, %+.2e ]\n", m_relf_range_left[i], m_relf_range_right[i], yl_name, yl, yr );
            }
        }
    }/* << SERVICE: RANGE FINDER */

    {/* >> SERVICE: PROGRESS BAR */
        f64_t count_exp  = log10( count_hi );
        i64_t count_iexp = (i64_t) count_exp;
        f64_t count_frac = ( count_exp > 1.0 ) ? ( count_exp - count_iexp ) : count_exp;

        printf("\n Analyzed points:\n"LINE_SEPARATOR"\n%.1fe+%lld (%.1f%%) in [ %e, %e ]\n", pow(10.0, count_frac), count_iexp, 100. * ( (f64_t) count_hi / ((f64_t) total_fp32) ), x_min, x_curr );
    }/* >> SERVICE: PROGRESS BAR */

    {/* >> SERVICE: ULP */
        ULPS_HIST_REPORT;
    }/* << SERVICE: ULP */
}

void ADDCALL fp64_range_analyzer( cstr_t desc, f6464_t lhs, f6464_t rhs, f6464_t next_x, f64_t accept, f64_t reject, u8_t state[ /* TBD: Size in bytes of the state */ ] )
{
    printf("Analyzing %s\n"LINE_SEPARATOR"\n", desc);

    /* >> MAIN ANALYSIS */
    const f64_t x_min   = -f64_get_named_fp_in_real_line( NAMED_FP_INF );
    const f64_t x_max   = +f64_get_named_fp_in_real_line( NAMED_FP_INF );

    f64_t x_curr        = x_min;
    /* << MAIN ANALYSIS */

    /* >> SERVICE: PROGRESS BAR */
    fp_progress_state_t pstate =
    {
        .xexp = EMIN_F64,
        .count_points_in_the_range = 0,
    };
    /* << SERVICE: PROGRESS BAR */

    /* >> SERVICE: POINTS COUNTER */
    u64_t count_hi = 0;
    u64_t total_fp64 = (~0llu - (2llu * 8388607llu) + 1llu ) ; /* 2^64 - all f64 NaNs */

    /* << SERVICE: POINTS COUNTER */

    /* >> SERVICE: ULP */
    ULPS_HIST_INIT;
    /* << SERVICE: ULP */

    /* >> SERVICE: RANGE FINDER */
    f64_t x_prev = x_curr;
    INIT_RANGE( x_prev, f64_get_subset( lhs ( x_prev ) ) );
    /* << SERVICE: RANGE FINDER */

    /* >> SERVICE: ERROR LOOP TERMINATION */
    i32_t x_curr_ntimes_analyzed = 0;
    /* << SERVICE: ERROR LOOP TERMINATION */

    /* >> MAIN ANALYSIS */
    while( ( x_curr <= x_max ) && ( x_curr_ntimes_analyzed == 0 ) )
    {
        /* >> MAIN ANALYSIS */
        const f64_t     y_lhs       = lhs ( x_curr );
        const f64_t     y_rhs       = rhs ( x_curr );
        const fpset_t   y_lhs_type  = f64_get_subset( y_lhs );
        /* << MAIN ANALYSIS */

        /* >> SERVICE: ULP */
        const f64_t curr_ulp = fabs( fp64_ulps( y_lhs, y_rhs ) );
        ULPS_HIST_RECORD( curr_ulp, x_curr );
        /* << SERVICE: ULP */

        {/* >> SERVICE: PROGRESS BAR */
            fp64_progress_bar( g_fpenv, x_curr, x_min, x_max, &pstate );
        }/* << SERVICE: PROGRESS BAR */

        {/* >> SERVICE: POINTS COUNTER */
            count_hi++;
        }/* << SERVICE: POINTS COUNTER */

        {/* >> SERVICE: RANGE FINDER */
            UPDATE_RANGE(x_prev, x_curr, y_lhs, y_lhs_type);
            x_prev = x_curr;
        }/* << SERVICE: RANGE FINDER */

        {/* >> MAIN ANALYSIS */
            x_curr = next_x( x_curr ); /* Update x */
        }/* >> MAIN ANALYSIS */

        {/* >> SERVICE: ERROR LOOP TERMINATION */
            x_curr_ntimes_analyzed = ( fp64_equals( x_curr, x_prev ) ? ( x_curr_ntimes_analyzed + 1 ) : 0 );
        }/* << SERVICE: ERROR LOOP TERMINATION */
    }
    /* << MAIN ANALYSIS */

    {/* >> SERVICE: PROGRESS BAR */
        console_erase_last_n_chars( 255u, g_fpenv );
    }/* << SERVICE: PROGRESS BAR */

    {/* >> SERVICE: ERROR LOOP TERMINATION */
        if( x_curr_ntimes_analyzed != 0 )
        {
            printf("[INFO]: loop stopped to prevent double analysis of point %f\n", x_curr);
        }
    }/* << SERVICE: ERROR LOOP TERMINATION */

    {/* >> SERVICE: RANGE FINDER */
    CLOSE_RANGE(x_prev);
    }/* << SERVICE: RANGE FINDER */

    /*
     * REPORTING SECTION
     * TODO: define the state of this function and send all of this to another function is possible
     */

    {/* >> SERVICE: RANGE FINDER */
        printf("\n Range Len:\n"LINE_SEPARATOR"\n");
        for(u32_t i = 0; i < m_relf_range_len; i++)
        {
            const f64_t yl = lhs( m_relf_range_left[i] );
            const f64_t yr = lhs( m_relf_range_right[i] );

            fpset_t yl_type = f64_get_subset( yl );
            cstr_t  yl_name = fp64_get_subset_name( yl );

            if( ( yl_type == FP_TYPE_NAN ) || ( yl_type == FP_TYPE_INFINITE ) || ( yl_type == FP_TYPE_ZERO ) )
            {
                printf("[ %+12.2e, %+12.2e ] -> %s\n", m_relf_range_left[i], m_relf_range_right[i], yl_name );
            }
            else
            {
                printf("[ %+12.2e, %+12.2e ] -> %s [ %+.2e, %+.2e ]\n", m_relf_range_left[i], m_relf_range_right[i], yl_name, yl, yr );
            }
        }
    }/* << SERVICE: RANGE FINDER */

    {/* >> SERVICE: PROGRESS BAR */
        f64_t count_exp  = log10( count_hi );
        i64_t count_iexp = (i64_t) count_exp;
        f64_t count_frac = ( count_exp > 1.0 ) ? ( count_exp - count_iexp ) : count_exp;

        printf("\n Analyzed points:\n"LINE_SEPARATOR"\n%.1fe+%lld (%.1f%%) in [ %e, %e ]\n", pow(10.0, count_frac), count_iexp, 100. * ( (f64_t) count_hi / ((f64_t) total_fp64) ), x_min, x_curr );
    }/* >> SERVICE: PROGRESS BAR */

    {/* >> SERVICE: ULP */
        ULPS_HIST_REPORT;
    }/* << SERVICE: ULP */
}
