#include <math.h>
#include "public_math.h"
#include "private_math.h"

static i16_t    math_pow_safe_exp   ( f64_t y                           );
static u32_t    math_pow_iabs       ( i32_t n                           );
static f64_t    math_pow_impl2       ( f64_t x,  f64_t y                 );
static f64_t    math_pow_i16        ( f64_t x,  i16_t n                 );
static f64_t    math_pow_ydexp      ( i64_t zl, f64_t y,    f64_t dexp  );
static i64_t    math_pow_clamp      ( i64_t x,  i64_t xmin, i64_t xmax  );

extern f64_t math_pow(f64_t x, f64_t y)
{
    dnorm_t result      = { 0 };

    result.type         = math_type( x );
    result.f.w[ W0 ]    = result.type;

    switch( result.type )
    {
        case (FINITE):
        {
            if( x > 0.0 )
            {
                if( x == 1.0 )
                {
                    result.f.d = 1.0;
                }

                else
                {
                    result.f.d = math_pow_impl2( x, y );
                }
            }

            else
            {
                if( math_intrnd( y ) == y )
                {
                   if( (__int128_t)y%2==0 ) /* TODO: we need a typedef of i128_t, to detect when y is even safely */
                    {
                        result.f.d = math_pow_impl2( -x, y );
                    }

                    else
                    {
                        result.f.d = -math_pow_impl2( -x, y );
                    }
                }

                else
                {
                    result.f.w[ W0 ] = NIL;
                }
            }

            break;
        }

        case (INF):
        {
            result.f.w[ W0 ] = INF;
            break;
        }

        case (NIL):
        case (GRADZ):
        {
            result.f.w[ W0 ] = NIL;
            break;
        }

        case (ZERO): /* x^0 = 1, for all x >= 0, even 0^0 */
        {
            if( math_type( y ) == ZERO )
            {
                result.f.d = 1.0;
            }

            else
            {
                result.f.d = 0.0;
            }

            break;
        }
    }

    return result.f.d;
}

/*
 * math_pow_safe_exp
 *
 * cast exponent y of function pow( x, y ) in a safe range of |x| < 1024
 * which corresponds with the range of 11bits of e in positive and negative direction: 2**{11-1} = 1023
 */

static i16_t math_pow_safe_exp( f64_t y )
{
    i16_t result;

    static const f64_t safe_exp = (f64_t) (DMAX>>1); /* 1023 */

    if( ( -safe_exp <= y ) && ( y <= +safe_exp ) )
    {
        result = (i16_t) y;
    }
    else
    {
        result = 0.0;
    }

    return ( result );
}

static f64_t math_pow_ydexp( i64_t zl, f64_t y, f64_t dexp )
{
    static const f64_t ln2 = 0.69314718055994530942;

    f64_t result;

    if ( zl == 0 )
    {
        result = y * dexp;
    }

    else
    {
        const f64_t yrnd = math_intrnd( y );

        result = ( ( yrnd * dexp ) - (f64_t) zl ) + ( y - yrnd ) * dexp;
    }

    return ( ln2 * result );
}

static f64_t math_pow_impl2( f64_t x, f64_t y )
{
/*  3. determine m */
/*  5. determine r, g */
/*  6. determine p */
/*  7. determine z */
/*  8. determine u2 = R(z) */
/*  9. determine w = y*u */
/* 10. determine w > BIGX */
/* 11. determine w < SMALLX */
/* 13. determine p', r', m' */
/* 14. determine z = n**w2 */
/* 15. determine result = z * B**m' * n**(-r'-p'/C) */
    return ( x*y );
}

static u32_t math_pow_iabs( i32_t n )
{
    u32_t result;

    if ( n > 0 )
    {
        result  = ( u32_t ) n;
    }

    else
    {
        result  = ( u32_t ) -n;
    }

    return ( result );
}

static f64_t math_pow_i16( f64_t x, i16_t n )
{
    f64_t result = 1.0; /* exact result when n = 0 */

    u32_t k = math_pow_iabs( n );

    f64_t x2 = x;

    while ( k > 0 )
    {
        if ( k & 1 )
        {
            result = result * x2;
        }

        k = k >> 1;

        if ( k > 0 )
        {
            x2 = x2 * x2;
        }
    }

    if ( n < 0 )
    {
        result = 1.0 / result;
    }

    return ( result );
}

static i64_t math_pow_clamp( i64_t x, i64_t xmin, i64_t xmax )
{
    i64_t result;

    if( x <= xmin )
    {
        result = xmin;
    }

    else if ( x < xmax )
    {
        result = x;
    }

    else
    {
        result = xmax;
    }

    return ( result );
}
