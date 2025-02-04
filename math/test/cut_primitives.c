#include "../vendor/Unity-v2.6.0/unity.h"
#include "../src/public_math.h"
#include "../src/math_primitives.c"
#include "../src/math_pow.c"

void setUp    ( void ) { /* NOP */ }
void tearDown ( void ) { /* NOP */ }

/*
 * CUT: math_type
 *
 * each test function is testing extremal values of a NORMAL RANGE,
 * the combination of all tests gives a conclusion about BVA
 */

#define INIT_WDOUBLE( w0, w1, w2, w3) .w[ W0 ] = (w0), .w[ W1 ] = (w1), .w[ W2 ] = (w2), .w[ W3 ] = (w3),

void test_math_type_zero( void )
{
    TEST_ASSERT_EQUAL_UINT16(ZERO, math_type(-0.0));
    TEST_ASSERT_EQUAL_UINT16(ZERO, math_type(+0.0));
}

void test_math_type_negative_underflow( void )
{
    TEST_ASSERT_EQUAL_UINT16(GRADZ, math_type( (dw_t) { INIT_WDOUBLE( 0x800F, 0xFFFF, 0xFFFF, 0xFFFF ) }.d )); /* -2.22507385850720088902459e-308 */
    TEST_ASSERT_EQUAL_UINT16(GRADZ, math_type( (dw_t) { INIT_WDOUBLE( 0x8000, 0x0000, 0x0000, 0x0001 ) }.d )); /* -4.94065645841246544176569e-324 */
}

void test_math_type_positive_underflow( void )
{
    TEST_ASSERT_EQUAL_UINT16(GRADZ, math_type( (dw_t) { INIT_WDOUBLE( 0x0000, 0x0000, 0x0000, 0x0001 ) }.d )); /* +4.94065645841246544176569e-324 */
    TEST_ASSERT_EQUAL_UINT16(GRADZ, math_type( (dw_t) { INIT_WDOUBLE( 0x000F, 0xFFFF, 0xFFFF, 0xFFFF ) }.d )); /* +2.22507385850720088902459e-308 */
}

void test_math_type_negative_finite( void )
{
    TEST_ASSERT_EQUAL_UINT16(FINITE, math_type( (dw_t) { INIT_WDOUBLE( 0x8010, 0x0000, 0x0000, 0x0001 ) }.d )); /* -2.22507385850720138309023e-308 */
    TEST_ASSERT_EQUAL_UINT16(FINITE, math_type( (dw_t) { INIT_WDOUBLE( 0xFFEF, 0xFFFF, 0xFFFF, 0xFFFF ) }.d )); /* -1.79769313486231570814527e+308 */
}

void test_math_type_positive_finite( void )
{
    TEST_ASSERT_EQUAL_UINT16(FINITE, math_type( (dw_t) { INIT_WDOUBLE( 0x0010, 0x0000, 0x0000, 0x0001 ) }.d )); /* +2.22507385850720138309023e-308 */
    TEST_ASSERT_EQUAL_UINT16(FINITE, math_type( (dw_t) { INIT_WDOUBLE( 0x7FEF, 0xFFFF, 0xFFFF, 0xFFFF ) }.d )); /* +1.79769313486231570814527e+308 */
}

void test_math_type_minus_inf( void )
{
    TEST_ASSERT_EQUAL_UINT16( INF, math_type(-1./0.) );
}

void test_math_type_plus_inf( void )
{
    TEST_ASSERT_EQUAL_UINT16( INF, math_type(+1./0.) );
}

void test_math_type_nan( void )
{
    f64_t plus_inf = (dw_t) { INIT_WDOUBLE( 0x7FF0, 0x0000, 0x0000, 0x0000 ) }.d;

    TEST_ASSERT_EQUAL_UINT16( NIL, math_type(0./0.) );
    TEST_ASSERT_EQUAL_UINT16( NIL, math_type(0. * plus_inf) );
    TEST_ASSERT_EQUAL_UINT16( NIL, math_type(plus_inf - plus_inf) );

    /* There are 1+52 bits of combinations to create different values for the NaN symbol, here are a few examples */

    TEST_ASSERT_EQUAL_UINT16( NIL, math_type( (dw_t) { INIT_WDOUBLE( 0x7FF0, 0x0000, 0x0000, 0x0001 ) }.d ) );
    TEST_ASSERT_EQUAL_UINT16( NIL, math_type( (dw_t) { INIT_WDOUBLE( 0x7FF0, 0x0000, 0x0000, 0x0010 ) }.d ) );
    TEST_ASSERT_EQUAL_UINT16( NIL, math_type( (dw_t) { INIT_WDOUBLE( 0x7FF0, 0x0000, 0x0000, 0x0100 ) }.d ) );
    TEST_ASSERT_EQUAL_UINT16( NIL, math_type( (dw_t) { INIT_WDOUBLE( 0x7FF0, 0x0000, 0x0000, 0x1000 ) }.d ) );

    TEST_ASSERT_EQUAL_UINT16( NIL, math_type( (dw_t) { INIT_WDOUBLE( 0x7FF0, 0x0000, 0x0001, 0x0000 ) }.d ) );
    TEST_ASSERT_EQUAL_UINT16( NIL, math_type( (dw_t) { INIT_WDOUBLE( 0x7FF0, 0x0000, 0x0010, 0x0000 ) }.d ) );
    TEST_ASSERT_EQUAL_UINT16( NIL, math_type( (dw_t) { INIT_WDOUBLE( 0x7FF0, 0x0000, 0x0100, 0x0000 ) }.d ) );
    TEST_ASSERT_EQUAL_UINT16( NIL, math_type( (dw_t) { INIT_WDOUBLE( 0x7FF0, 0x0000, 0x1000, 0x0000 ) }.d ) );

    TEST_ASSERT_EQUAL_UINT16( NIL, math_type( (dw_t) { INIT_WDOUBLE( 0x7FF0, 0x0001, 0x0000, 0x0000 ) }.d ) );
    TEST_ASSERT_EQUAL_UINT16( NIL, math_type( (dw_t) { INIT_WDOUBLE( 0x7FF0, 0x0010, 0x0000, 0x0000 ) }.d ) );
    TEST_ASSERT_EQUAL_UINT16( NIL, math_type( (dw_t) { INIT_WDOUBLE( 0x7FF0, 0x0100, 0x0000, 0x0000 ) }.d ) );
    TEST_ASSERT_EQUAL_UINT16( NIL, math_type( (dw_t) { INIT_WDOUBLE( 0x7FF0, 0x1000, 0x0000, 0x0000 ) }.d ) );

    TEST_ASSERT_EQUAL_UINT16( NIL, math_type( (dw_t) { INIT_WDOUBLE( 0xFFF0, 0x0000, 0x0000, 0x0001 ) }.d ) );
    TEST_ASSERT_EQUAL_UINT16( NIL, math_type( (dw_t) { INIT_WDOUBLE( 0xFFF0, 0x0000, 0x0000, 0x0010 ) }.d ) );
    TEST_ASSERT_EQUAL_UINT16( NIL, math_type( (dw_t) { INIT_WDOUBLE( 0xFFF0, 0x0000, 0x0000, 0x0100 ) }.d ) );
    TEST_ASSERT_EQUAL_UINT16( NIL, math_type( (dw_t) { INIT_WDOUBLE( 0xFFF0, 0x0000, 0x0000, 0x1000 ) }.d ) );

    TEST_ASSERT_EQUAL_UINT16( NIL, math_type( (dw_t) { INIT_WDOUBLE( 0xFFF0, 0x0000, 0x0001, 0x0000 ) }.d ) );
    TEST_ASSERT_EQUAL_UINT16( NIL, math_type( (dw_t) { INIT_WDOUBLE( 0xFFF0, 0x0000, 0x0010, 0x0000 ) }.d ) );
    TEST_ASSERT_EQUAL_UINT16( NIL, math_type( (dw_t) { INIT_WDOUBLE( 0xFFF0, 0x0000, 0x0100, 0x0000 ) }.d ) );
    TEST_ASSERT_EQUAL_UINT16( NIL, math_type( (dw_t) { INIT_WDOUBLE( 0xFFF0, 0x0000, 0x1000, 0x0000 ) }.d ) );

    TEST_ASSERT_EQUAL_UINT16( NIL, math_type( (dw_t) { INIT_WDOUBLE( 0xFFF0, 0x0001, 0x0000, 0x0000 ) }.d ) );
    TEST_ASSERT_EQUAL_UINT16( NIL, math_type( (dw_t) { INIT_WDOUBLE( 0xFFF0, 0x0010, 0x0000, 0x0000 ) }.d ) );
    TEST_ASSERT_EQUAL_UINT16( NIL, math_type( (dw_t) { INIT_WDOUBLE( 0xFFF0, 0x0100, 0x0000, 0x0000 ) }.d ) );
    TEST_ASSERT_EQUAL_UINT16( NIL, math_type( (dw_t) { INIT_WDOUBLE( 0xFFF0, 0x1000, 0x0000, 0x0000 ) }.d ) );
}

void test_math_significand_shall_set_to_zero_all_non_significand_bits_from_input( void )
{
    dw_t given      = { .w[ W0 ] = 0xffff, .w[ W1 ] = 0xffff, .w[ W2 ] = 0xffff, .w[ W3 ] = 0xffff, };
    dw_t expected   = { .w[ W0 ] = 0x800f, .w[ W1 ] = 0xffff, .w[ W2 ] = 0xffff, .w[ W3 ] = 0xffff, };

    TEST_ASSERT_EQUAL_DOUBLE( expected.d, math_significand( given.d ) );
}

void test_math_cwnormalize( void )
{
    /*
     * ===================================================================
     *  alternative method to compute normalized values for |x| >= 1/B
     * ===================================================================
     *
     * OBJECTIVE:
     *
     * determine the expected result of math_cwnormalize( x )
     *
     * DESCRIPTION:
     *
     * expressing 'x' in the normalized format
     *
     * -------------------------------------------------------------------
     *                  x = f·B**e
     *                  f in [1/B, 1)       as fp               (eq.0)
     *                  e in [emin, emax]   as integer
     * -------------------------------------------------------------------
     *
     * solving (eq.0) the significand 'f' is uniquely determined by (eq.1)
     *
     * -------------------------------------------------------------------
     *                  f = x·B**-e
     *                  f in [1/B, 1)       as fp               (eq.1)
     *                  e in [emin, emax]   as integer
     * -------------------------------------------------------------------
     *
     * To determine the value of 'f', (eq.1) is transformed according to
     * o* := logB( o ), where B is the radix-B of the fp representation
     *
     * -------------------------------------------------------------------
     *                  f* = x* - e
     *                  f* in [-1, 0)       as fp               (eq.2)
     *                  e in [emin, emax]   as integer
     * -------------------------------------------------------------------
     *
     * rhs of (eq.2) is constrained to f* range [-1, 0), thus e satisfies:
     *
     * -------------------------------------------------------------------
     *                  -1-x* <= e < -x*                        (eq.4)
     * -------------------------------------------------------------------
     *
     * according to (eq.4) the solution of e must be a value of -1-x*, or
     * one unit less than that, which corresponds to (eq.5), where
     * operator [o] corresponds to rounding towards zero.
     *
     * -------------------------------------------------------------------
     *                  e = [-1-x*]                             (eq.5)
     * -------------------------------------------------------------------
     *
     * finally the solution to normalized significand f is given by (eq.6)
     *
     * -------------------------------------------------------------------
     *                  f = x·B**-[1+x*], with |x| >= 1/B       (eq.6)
     * -------------------------------------------------------------------
     */

    #define expect_cwnormalize(_x_) (_x_)*pow(2.,-(long long)(1+log2(_x_))) /* valid for x >= |1/B|, and arbitrarily for x=0 */

    /* POSITIVE BRANCH */
    for( double x = +0.5; x < +1000.; x += 0.001 )
    TEST_ASSERT_EQUAL_DOUBLE( +expect_cwnormalize(+x), math_cwnormalize(+x).f.d );

    /* NEGATIVE BRANCH */
    for( double x = -0.5; x > -1000.; x -= 0.001 )
    TEST_ASSERT_EQUAL_DOUBLE( -expect_cwnormalize(-x), math_cwnormalize(+x).f.d );

    /* POSITIVE BRANCH with 0 <= x < |1/B| */
    for( double x = +0.0; x < +0.5; x += 0.001 )
    TEST_ASSERT_EQUAL_DOUBLE( +2.*expect_cwnormalize(+x), math_cwnormalize(+x).f.d );

    /* NEGATIVE BRANCH with -|1/B| < x <= 0 */
    for( double x = -0.0; x > -0.5; x -= 0.001 )
    TEST_ASSERT_EQUAL_DOUBLE( -2.*expect_cwnormalize(-x), math_cwnormalize(+x).f.d );

    #undef expect_cwnormalize
}

void test_math_cwnormalize_dunderflow_left( void )
{
    /* GRADUAL UNDEFLOW x=2.225074e-308 */
    dw_t x = { .w[ W0 ] = 0x000F, .w[ W1 ] = 0xFFFF, .w[ W2 ] = 0xFFFF, .w[ W3 ] = 0xFFFF };
    TEST_ASSERT_EQUAL_DOUBLE( 0.50000003179507908825935950226121, math_cwnormalize(+x.d).f.d );
}

void test_math_cwnormalize_dunderflow_right( void )
{
    /* GRADUAL UNDEFLOW x=2.225042e-308 */
    dw_t x = { .w[ W0 ] = 0x0000, .w[ W1 ] = 0x0000, .w[ W2 ] = 0x0000, .w[ W3 ] = 0x0001 };
    TEST_ASSERT_EQUAL_DOUBLE( 0.9999856820450792779914571700809, math_cwnormalize(+x.d).f.d );
}

void test_math_cwnormalize_double_underflow_left( void )
{
    /* GRADUAL UNDEFLOW x=2.225074e-308 */
    dw_t x = { .w[ W0 ] = 0x000F, .w[ W1 ] = 0xFFFF, .w[ W2 ] = 0xFFFF, .w[ W3 ] = 0xFFFF };
    TEST_IGNORE_MESSAGE("Peding to implement Normalization Underflow");
    TEST_ASSERT_EQUAL_DOUBLE( 0.50000003179507908825935950226121, math_cwnormalize(+x.d).f.d );
}

void test_math_cwnormalize_double_underflow_right( void )
{
    /* GRADUAL UNDEFLOW x=2.225042e-308 */
    dw_t x = { .w[ W0 ] = 0x0000, .w[ W1 ] = 0x0000, .w[ W2 ] = 0x0000, .w[ W3 ] = 0x0001 };
    TEST_IGNORE_MESSAGE("Peding to implement Normalization Underflow");
    TEST_ASSERT_EQUAL_DOUBLE( 0.9999856820450792779914571700809, math_cwnormalize(+x.d).f.d );
}

void test_math_cwnormalize_float_underflow_left( void )
{
    /* GRADUAL UNDEFLOW x=1.401298464324817070923730e-45 */
    f32_t x = 1.401298464324817070923730e-45;
    TEST_ASSERT_EQUAL_DOUBLE( 0.5, math_cwnormalize(+x).f.d );
}

void test_math_cwnormalize_float_underflow_right( void )
{
    /* GRADUAL UNDEFLOW x=1.175494210692441075487029-38 */
    f32_t x = 1.175494210692441075487029e-38;
    TEST_ASSERT_EQUAL_DOUBLE( 0.99999988079071044921874962156408, math_cwnormalize(+x).f.d );
}

void test_math_cwsetexp( void )
{
    #define expect_cwsetexp(_f_,_n_) math_cwnormalize((_f_)).f.d * pow(2.0, (_n_))

    /* test math_cwsetexp( f, n ), with 0.5 <= f < 1, and -1021 <= n <= +1023 */
    for( int n = -1021; n <= 1023; n++ )
    for( double f =-10.; f < 10.; f += 0.1 )
    TEST_ASSERT_EQUAL_DOUBLE( expect_cwsetexp(f, n), math_cwsetexp(f, n) );

    /* test some specific points */
    TEST_ASSERT_EQUAL_DOUBLE( expect_cwsetexp(+0.5, 2), math_cwsetexp(+0.5, 2) );
    TEST_ASSERT_EQUAL_DOUBLE( expect_cwsetexp(-0.5, 2), math_cwsetexp(-0.5, 2) );

    TEST_IGNORE_MESSAGE("Skip til add test to all TYPES");
    TEST_ASSERT_EQUAL_DOUBLE( expect_cwsetexp(+0.0, 2), math_cwsetexp(+0.0, 2) );
    TEST_ASSERT_EQUAL_DOUBLE( expect_cwsetexp(-0.0, 2), math_cwsetexp(-0.0, 2) );

    #undef expect_cwsetexp
}

void test_math_intrnd( void )
{
    /* NORMAL */
    TEST_ASSERT_EQUAL_DOUBLE( +1.00000e00, math_intrnd( +1.10000e00 ) );
    TEST_ASSERT_EQUAL_DOUBLE( +1.00000e01, math_intrnd( +1.01000e01 ) );
    TEST_ASSERT_EQUAL_DOUBLE( +1.00000e02, math_intrnd( +1.00100e02 ) );
    TEST_ASSERT_EQUAL_DOUBLE( +1.00000e03, math_intrnd( +1.00010e03 ) );
    TEST_ASSERT_EQUAL_DOUBLE( +1.00000e04, math_intrnd( +1.00001e04 ) );
    TEST_ASSERT_EQUAL_DOUBLE( -1.00000e00, math_intrnd( -1.10000e00 ) );
    TEST_ASSERT_EQUAL_DOUBLE( -1.00000e01, math_intrnd( -1.01000e01 ) );
    TEST_ASSERT_EQUAL_DOUBLE( -1.00000e02, math_intrnd( -1.00100e02 ) );
    TEST_ASSERT_EQUAL_DOUBLE( -1.00000e03, math_intrnd( -1.00010e03 ) );
    TEST_ASSERT_EQUAL_DOUBLE( -1.00000e04, math_intrnd( -1.00001e04 ) );

    /* BVA */
    TEST_ASSERT_EQUAL_DOUBLE( 0.0, math_intrnd( -0.00000001 ) );
    TEST_ASSERT_EQUAL_DOUBLE( 0.0, math_intrnd( +0.00000001 ) );

    TEST_ASSERT_EQUAL_DOUBLE(  0.0  , math_intrnd( +0.499999    ) );
    TEST_ASSERT_EQUAL_DOUBLE( +1.0  , math_intrnd( +0.500001    ) );

    TEST_ASSERT_EQUAL_DOUBLE(  0.0  , math_intrnd( -0.499999    ) );
    TEST_ASSERT_EQUAL_DOUBLE( -1.0  , math_intrnd( -0.500001    ) );

    TEST_ASSERT_EQUAL_DOUBLE(  5.0  , math_intrnd( +5.499999    ) );
    TEST_ASSERT_EQUAL_DOUBLE( +6.0  , math_intrnd( +5.500001    ) );

    TEST_ASSERT_EQUAL_DOUBLE( -5.0  , math_intrnd( -5.499999    ) );
    TEST_ASSERT_EQUAL_DOUBLE( -6.0  , math_intrnd( -5.500001    ) );

    /* ROBUSTNESS */
    TEST_ASSERT_EQUAL_DOUBLE( 0.0   , math_intrnd( (dw_t){ .w[ W3 ] = 1 }.d )   ); /* UNDERFLOW */
    TEST_ASSERT_EQUAL_DOUBLE( 0.0   , math_intrnd( 0.0   )                      ); /* ZERO returns ZERO */
    TEST_ASSERT_EQUAL_DOUBLE( 1./0. , math_intrnd( 1./0. )                      ); /* INF  returns INF  */
    TEST_ASSERT_EQUAL_DOUBLE( 0./0. , math_intrnd( 0./0. )                      ); /* NAN  returns NAN  */

    /* RANGE, I'm not sure about how to test it */
    TEST_ASSERT_EQUAL_DOUBLE( 3.402823466385288598117042E+38  , math_intrnd( 3.402823466385288598117042E+38 ));
}

void test_math_cwnormalized_exponent( void )
{
    TEST_ASSERT_EQUAL_INT16( -2, math_cwnormalize( 0.125 ).e );
    TEST_ASSERT_EQUAL_INT16( -1, math_cwnormalize( 0.25  ).e );
    TEST_ASSERT_EQUAL_INT16(  0, math_cwnormalize( 0.50  ).e );
    TEST_ASSERT_EQUAL_INT16( +1, math_cwnormalize( 1.0   ).e );
    TEST_ASSERT_EQUAL_INT16( +2, math_cwnormalize( 2.0   ).e );
    TEST_ASSERT_EQUAL_INT16( +3, math_cwnormalize( 4.0   ).e );
}

void test_math_cwnormalized_exponent_range( void )
{
    #define set_cwnormalized(_f_,_n_) (_f_)*pow(2., (_n_))

    #define emin 1022
    #define emax 1023

    for(int i = -emin + 1; i < emax; i++)
    TEST_ASSERT_EQUAL_INT16( i, math_cwnormalize( set_cwnormalized(0.5, i) ).e );

    for(int i = -emin + 1; i < emax; i++)
    TEST_ASSERT_EQUAL_INT16( i, math_cwnormalize( set_cwnormalized(0.9, i) ).e );

    #undef set_cwnormalized
}

void test_math_scale( void )
{

}

int main( void )
{
    UNITY_BEGIN();

    RUN_TEST( test_math_type_zero );
    RUN_TEST( test_math_type_negative_underflow );
    RUN_TEST( test_math_type_positive_underflow );
    RUN_TEST( test_math_type_negative_finite );
    RUN_TEST( test_math_type_positive_finite );

    RUN_TEST( test_math_type_plus_inf );
    RUN_TEST( test_math_type_minus_inf );
    RUN_TEST( test_math_type_nan );

    RUN_TEST( test_math_significand_shall_set_to_zero_all_non_significand_bits_from_input );

    RUN_TEST( test_math_cwnormalize );
    RUN_TEST( test_math_cwnormalize_double_underflow_left );
    RUN_TEST( test_math_cwnormalize_double_underflow_right );
    RUN_TEST( test_math_cwnormalize_float_underflow_left );
    RUN_TEST( test_math_cwnormalize_float_underflow_right );
    RUN_TEST( test_math_cwnormalized_exponent );
    RUN_TEST( test_math_cwnormalized_exponent_range );

    RUN_TEST( test_math_cwsetexp );

    RUN_TEST( test_math_intrnd );

    RUN_TEST( test_math_scale );

    return UNITY_END();
}
