#include <stdio.h>
#include "../vendor/Unity-v2.6.0/unity.h"
#include "../src/math_primitives.c"
#include "../src/math_arithmetic.c"

void setUp    ( void ) { /* NOP */ }
void tearDown ( void ) { /* NOP */ }

/*
 * CUT: math_add_u64
 */

void math_add_u64_base_cases( void )
{
    TEST_ASSERT_EQUAL_UINT64( 2  , math_add_u64( 1, 1, NULL ) );
    TEST_ASSERT_EQUAL_UINT64( 3  , math_add_u64( 1, 2, NULL ) );
    TEST_ASSERT_EQUAL_UINT64( 101, math_add_u64( 1, 100, NULL ) );
    TEST_ASSERT_EQUAL_UINT64( 1024, math_add_u64( 1000, 24, NULL ) );
    TEST_ASSERT_EQUAL_UINT64( 5032, math_add_u64(32, 5000, NULL ) );
}

void math_add_u64_with_prev_carry( void )
{
    bool_t carry;

    carry = 1;
    TEST_ASSERT_EQUAL_UINT64( 3  , math_add_u64( 1, 1, &carry ) );

    carry = 1;
    TEST_ASSERT_EQUAL_UINT64( 4  , math_add_u64( 1, 2, &carry ) );

    carry = 1;
    TEST_ASSERT_EQUAL_UINT64( 102, math_add_u64( 1, 100, &carry ) );

    carry = 1;
    TEST_ASSERT_EQUAL_UINT64( 1025, math_add_u64( 1000, 24, &carry ) );

    carry = 1;
    TEST_ASSERT_EQUAL_UINT64( 5033, math_add_u64(32, 5000, &carry ) );
}

void math_add_u64_with_last_carry( void )
{
    bool_t carry = 0;

    TEST_ASSERT_EQUAL_UINT64( (1ULL<<33)  , math_add_u64( (1ULL<<32), (1ULL<<32), &carry ) );
    TEST_ASSERT_EQUAL_INT( 0, carry );

    TEST_ASSERT_EQUAL_UINT64( (1ULL<<63)  , math_add_u64( (1ULL<<62), (1ULL<<62), &carry ) );
    TEST_ASSERT_EQUAL_INT( 0, carry );

    TEST_ASSERT_EQUAL_UINT64( 0  , math_add_u64( (1ULL<<63), (1ULL<<63), &carry ) );
    TEST_ASSERT_EQUAL_INT( 1, carry );
}

/*
 * CUT: math_add
 */

void math_add_check_trivial_cases( void )
{
    const f64_t pinf = +1./0.;
    const f64_t ninf = -1./0.;
    const f64_t nan  =  0./0.;

    TEST_ASSERT_EQUAL_UINT16( INF, math_type( math_add( pinf, pinf, NULL ) ) );
    TEST_ASSERT_EQUAL_UINT16( INF, math_type( math_add( ninf, ninf, NULL ) ) );

    TEST_ASSERT_EQUAL_UINT16( NIL, math_type( math_add( pinf, ninf, NULL ) ) );
    TEST_ASSERT_EQUAL_UINT16( NIL, math_type( math_add( ninf, pinf, NULL ) ) );

    TEST_ASSERT_EQUAL_UINT16( NIL, math_type( math_add( nan , 1.0 , NULL ) ) );
    TEST_ASSERT_EQUAL_UINT16( NIL, math_type( math_add( 1.0 , nan , NULL ) ) );

    TEST_ASSERT_EQUAL_DOUBLE( 1.0 , math_add( 0.0 , 1.0 , NULL ) );
    TEST_ASSERT_EQUAL_DOUBLE( 1.0 , math_add( 1.0 , 0.0 , NULL ) );
}

void math_add_check_impl( void )
{
    TEST_ASSERT_EQUAL_DOUBLE( 2.000000 , math_add( 1.000000 , 1.000000 , NULL ) );
    TEST_ASSERT_EQUAL_DOUBLE( 4.000003 , math_add( 2.000001 , 2.000002 , NULL ) );
    TEST_ASSERT_EQUAL_DOUBLE( 4.000003 , math_add( 2.000001 , 2.000002 , NULL ) );
    TEST_ASSERT_EQUAL_DOUBLE( 310.0    , math_add( 155.0000 , 155.0000 , NULL ) );
    TEST_ASSERT_EQUAL_DOUBLE( 2.50000 , math_add( 1.1, 1.4, NULL ) );

}

void math_add_with_diff_exp( void )
{
    /* TODO: math_add_impl does not work completely because, renormalization and packing is missing */
    TEST_ASSERT_EQUAL_DOUBLE( 3.0       , math_add( 1.0 , 2.0 , NULL ) );
}

int main( void )
{
    UNITY_BEGIN();

    RUN_TEST( math_add_u64_base_cases );
    RUN_TEST( math_add_u64_with_prev_carry );
    RUN_TEST( math_add_u64_with_last_carry );

    RUN_TEST( math_add_check_trivial_cases );
    RUN_TEST( math_add_check_impl );
    RUN_TEST( math_add_with_diff_exp );

    return ( UNITY_END() );
}
