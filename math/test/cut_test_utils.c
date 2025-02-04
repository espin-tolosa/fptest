#include "math.h"

#include "../vendor/Unity-v2.6.0/unity.h"

#include "../test/test_utils.h"

void setUp    ( void ) { /* NOP */ }
void tearDown ( void ) { /* NOP */ }

double sinx( double x ) { return ( ( -HUGE_SINX < x && x < +HUGE_SINX ) ? sin(x) : tmath_range_init( TD_XINIT_NAN ) );  }

/*
 * CUT: test_utils
 */

/* Requirement Based Tests */

void test_math_sin_plus_inf_is_nan( void )
{
    TEST_ASSERT_EQUAL_UINT16( NIL, math_type( math_sin( tmath_range_init( TD_PLUS_INF ) ) ) );
}

void test_math_sin_minus_inf_is_nan( void )
{
    TEST_ASSERT_EQUAL_UINT16( NIL, math_type( math_sin( tmath_range_init( TD_MINUS_INF ) ) ) );
}

/* Range Tests */
void test_math_sin_matches_double_precision( void )
{
    vector_t vtable = kv_new( 0 );
    (void) kv_push( &vtable, math_sin, "math_sin" );

    tmath_maxerr_t status = tmath_test_maxerr_f32_range(&vtable, math_sin, sinx, 0.0, RESULT_FILENAME );

    TEST_ASSERT_EQUAL_DOUBLE( 0.0, status.maxerr );

    system("cat " RESULT_FILENAME);
}

int main( void )
{
    UNITY_BEGIN();

    RUN_TEST( test_math_sin_plus_inf_is_nan );
    RUN_TEST( test_math_sin_plus_inf_is_nan );
    RUN_TEST( test_math_sin_matches_double_precision );

    return UNITY_END();
}
