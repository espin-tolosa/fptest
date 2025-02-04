#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "../vendor/fptest/fptest.h"
#include "../src/math_primitives.c"
#include "../src/math_sqrt.c"
#include "../src/math_exp.c"
#include "../src/math_log.c"
#include "../src/math_sin.c"
#include "../src/math_asin.c"
#include "../src/math_pow.c"
#include "../src/math_correlation.c"
#include "../src/math_surface.c"

/* CLEAR TEST TYPE SELECTOR */

/* RUN TEST TYPE SELECTOR */
#undef BENCHMARK_ANALISYS
#undef RANGE_ANALISYS
#undef ANALYZE_SQRT
#undef ANALYZE_SQRT
#undef ANALYZE_EXP
#undef ANALYZE_LOG
#undef ANALYZE_SIN
#undef ANALYZE_COS
#undef ANALYZE_ASIN
#undef ANALYZE_POW

#undef ANALYZE_CORRELATION
#undef ANALYZE_SURFACE

/* TOOLS */
static int load_array_from_file(const char *filename, double *array, int size);
static double read_pearson_correlation(const char *filename);
static void write_xy(const char *filename, const double * x, const double * y, unsigned n);
#define TEST_ASSERT_EQUAL_DOUBLE(x,y) do{ f64_t expect = (x) ; f64_t given = (y); if( !fp64_equals ( given, expect ) ){ fprintf(stderr, "[FAILED] %s, expect: %e, given: %e, lc:%d\n",__func__, expect, given, __LINE__ ); } } while(0)

/* TEST FUNCTIONS */
void test_math_cwsetexp( f64_t (*tested_cwsetexp)(f64_t, i16_t), bool_t active );

/* WRAPPERS */
static f64_t x, y, m_y; f64_t math_powy( f64_t x ); f64_t powy( f64_t x ); f64_t sinx( f64_t x ); f64_t cosx( f64_t x );
#define PRINT(x,fun,y) printf( #fun"(%.10e) = %.10e, expect: %.10e, ULPs = %f\n\n", x, fun(x), y, fp64_ulps(fun(x),y) ) /* example: PRINT(x,math_sqrt,sqrt(x));*/

int main( void )
{
    srand(time(NULL));

    g_fpenv.log_completation                    = TRUE;
    g_fpenv.log_completation_remove_history     = TRUE;
    g_fpenv.fp64_grow_frac                      = 1.e-7;
    g_fpenv.fp64_ctrl_ulps                      = 1000.;

    g_fpenv.fp64_range_sin_min = -1.e3;
    g_fpenv.fp64_range_sin_max = +1.e3;

    g_fpenv.fp64_range_cos_min = -1.e3;
    g_fpenv.fp64_range_cos_max = +1.e3;

    g_fpenv.fp64_sqrt_rejected  = 2.0;
    g_fpenv.fp64_exp_rejected   = 4.0;
    g_fpenv.fp64_log_rejected   = 2.0;
    g_fpenv.fp64_sin_rejected   = 2.0;
    g_fpenv.fp64_cos_rejected   = 2.0;
    g_fpenv.fp64_asin_rejected  = 3.0;
    g_fpenv.fp64_powy_rejected  = 10.0;

    fp64_test_sqrt( math_sqrt, TRUE );
    fp64_test_exp ( math_exp , TRUE );
    fp64_test_log ( math_log , TRUE );
    fp64_test_sin ( math_sin , TRUE );
    fp64_test_cos ( math_cos , TRUE );
    fp64_test_asin( math_asin, TRUE );
    fp64_test_pow ( math_pow , FALSE );

    test_math_cwsetexp( math_cwsetexp, TRUE );
    return 0;

#ifdef ANALYZE_POW
    x = -1./0.;

    y = -1.0 ;  assert( fp64_ulps( math_pow(x, y), pow(x, y) ) == 0.0 );
    y = -3.0 ;  assert( fp64_ulps( math_pow(x, y), pow(x, y) ) == 0.0 );
    y = -5.0 ;  assert( fp64_ulps( math_pow(x, y), pow(x, y) ) == 0.0 );

    x = +1./0.;

    y = -1.0 ;  assert( fp64_ulps( math_pow(x, y), pow(x, y) ) == 0.0 );
    y = -2.0 ;  assert( fp64_ulps( math_pow(x, y), pow(x, y) ) == 0.0 );
    y = -3.0 ;  assert( fp64_ulps( math_pow(x, y), pow(x, y) ) == 0.0 );
    y = -4.0 ;  assert( fp64_ulps( math_pow(x, y), pow(x, y) ) == 0.0 );
    y = -4.5 ;  assert( fp64_ulps( math_pow(x, y), pow(x, y) ) == 0.0 );
    y = -5.3 ;  assert( fp64_ulps( math_pow(x, y), pow(x, y) ) == 0.0 );

    y = -0.0;

    x = -1./0.; assert( fp64_ulps( math_pow(x, y), pow(x, y) ) == 0.0 );
    x = +1./0.; assert( fp64_ulps( math_pow(x, y), pow(x, y) ) == 0.0 );

    x = -1./0.;

    y = +1.0 ;  assert( fp64_ulps( math_pow(x, y), pow(x, y) ) == 0.0 );
    y = +3.0 ;  assert( fp64_ulps( math_pow(x, y), pow(x, y) ) == 0.0 );
    y = +5.0 ;  assert( fp64_ulps( math_pow(x, y), pow(x, y) ) == 0.0 );

    x = -1./0.;

    y = +1.0 ;  assert( fp64_ulps( math_pow(x, y), pow(x, y) ) == 0.0 );
    y = +2.0 ;  assert( fp64_ulps( math_pow(x, y), pow(x, y) ) == 0.0 );
    y = +3.0 ;  assert( fp64_ulps( math_pow(x, y), pow(x, y) ) == 0.0 );
    y = +4.0 ;  assert( fp64_ulps( math_pow(x, y), pow(x, y) ) == 0.0 );
    y = +4.5 ;  assert( fp64_ulps( math_pow(x, y), pow(x, y) ) == 0.0 );
    y = +5.3 ;  assert( fp64_ulps( math_pow(x, y), pow(x, y) ) == 0.0 );

    x = +1./0.;

    y = -6.0 ;  assert( fp64_ulps( math_pow(x, y), pow(x, y) ) == 0.0 );
    y = -5.3 ;  assert( fp64_ulps( math_pow(x, y), pow(x, y) ) == 0.0 );
    y = -4.5 ;  assert( fp64_ulps( math_pow(x, y), pow(x, y) ) == 0.0 );
    y = -4.0 ;  assert( fp64_ulps( math_pow(x, y), pow(x, y) ) == 0.0 );
    y = -3.75;  assert( fp64_ulps( math_pow(x, y), pow(x, y) ) == 0.0 );
    y = -3.0 ;  assert( fp64_ulps( math_pow(x, y), pow(x, y) ) == 0.0 );
    y = -2.0 ;  assert( fp64_ulps( math_pow(x, y), pow(x, y) ) == 0.0 );
    y = -1.0 ;  assert( fp64_ulps( math_pow(x, y), pow(x, y) ) == 0.0 );
    y = -0.0 ;  assert( fp64_ulps( math_pow(x, y), pow(x, y) ) == 0.0 );

    y = +6.0 ;  assert( fp64_ulps( math_pow(x, y), pow(x, y) ) == 0.0 );
    y = +5.3 ;  assert( fp64_ulps( math_pow(x, y), pow(x, y) ) == 0.0 );
    y = +4.5 ;  assert( fp64_ulps( math_pow(x, y), pow(x, y) ) == 0.0 );
    y = +4.0 ;  assert( fp64_ulps( math_pow(x, y), pow(x, y) ) == 0.0 );
    y = +3.75;  assert( fp64_ulps( math_pow(x, y), pow(x, y) ) == 0.0 );
    y = +3.0 ;  assert( fp64_ulps( math_pow(x, y), pow(x, y) ) == 0.0 );
    y = +2.0 ;  assert( fp64_ulps( math_pow(x, y), pow(x, y) ) == 0.0 );
    y = +1.0 ;  assert( fp64_ulps( math_pow(x, y), pow(x, y) ) == 0.0 );
    y = +0.0 ;  assert( fp64_ulps( math_pow(x, y), pow(x, y) ) == 0.0 );

    m_y = -1.e4;
    fp64_range_analyzer("powy -10E4", math_powy, powy , fp64_geometric_grow, 1.0, 8.0, NULL);

    m_y = -2.0;
    fp64_range_analyzer("powy -2", math_powy, powy , fp64_geometric_grow, 1.0, 8.0, NULL);

    m_y = -1.3;
    fp64_range_analyzer("powy -1.3", math_powy, powy , fp64_geometric_grow, 1.0, 8.0, NULL);

    m_y = -1.0;
    fp64_range_analyzer("powy -1", math_powy, powy , fp64_geometric_grow, 1.0, 8.0, NULL);

    m_y = -0.0;
    fp64_range_analyzer("powy -0", math_powy, powy , fp64_geometric_grow, 1.0, 8.0, NULL);

    m_y = +0.0;
    fp64_range_analyzer("powy +0", math_powy, powy , fp64_geometric_grow, 1.0, 8.0, NULL);

    m_y = -1.e-309;
    fp64_range_analyzer("powy -sub", math_powy, powy , fp64_geometric_grow, 1.0, 8.0, NULL);

    m_y = +1.e-309;
    fp64_range_analyzer("powy +sub", math_powy, powy , fp64_geometric_grow, 1.0, 8.0, NULL);

    m_y = +1.0;
    fp64_range_analyzer("powy +1", math_powy, powy , fp64_geometric_grow, 1.0, 8.0, NULL);

    m_y = 2.00;
    fp64_range_analyzer("powy 2", math_powy, powy , fp64_geometric_grow, 1.0, 8.0, NULL);

    m_y = 3.75;
    fp64_range_analyzer("powy 3.75", math_powy, powy , fp64_geometric_grow, 1.0, 8.0, NULL);

    m_y = 1.e4;
    fp64_range_analyzer("powy 10E4", math_powy, powy , fp64_geometric_grow, 1.0, 8.0, NULL);
#   endif

#   ifdef ANALYZE_CORRELATION
    {
        f64_t k    = 1./128.;
    do
    {
        #define ARRAY_SIZE (64/*1*60*60*32*/)

        f64_t N = ARRAY_SIZE;

        // Function to load an array from a file
        static double x[ARRAY_SIZE];
        static double y[ARRAY_SIZE];


        for( u32_t i = 0; i < (ARRAY_SIZE -1); i+=2 )
        {
            x[i  ] = (f64_t)i/2.;
            x[i+1] = x[i];

            y[i  ] = x[i];
            y[i+1] = y[i] + k;

//            y[i] = k*x[i] + fp64_rand_in_range(-1., +1.);
        }

        write_xy("xy.csv", x,y,SIZEOF(x));
//        system("python .\\test\\pearson_correlation.py xy.csv ");

        f64_t report_elapsed;

        f64_t rxy_v1;
        f64_t rxy_v2;
        {
            struct timeval t1, t2;
            double elapsedTime;

            // start timer
            mingw_gettimeofday(&t1, NULL);

            // do something
            // ...
            rxy_v1 = math_correlation( x, y, SIZEOF( x ) );
//            rxy_v2 = math_pearson_square( x, y, SIZEOF( x ) );

            // stop timer
            mingw_gettimeofday(&t2, NULL);

            // compute and print the elapsed time in millisec
            elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0;      // sec to ms
            elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;   // us to ms

            report_elapsed = elapsedTime;
        }

        f64_t expect = 1./( sqrt(1.+(0.1876*k*k*N/(N*N-1.))));

//        f64_t python = read_pearson_correlation("pearson_correlation.txt");
        f64_t ulp_v1 = ( fp64_ulps( rxy_v1, expect ) );

        printf("%+5.f, %e, %.15e, %.15e\n",report_elapsed, k, expect, rxy_v1);
//        printf("%+5.f, %1.e, %+.20e, %+2.f,\n", report_elapsed, k, rxy_v1, ulp_v1 );

        k = 2.*k;
    }
    while(k < 1.e+20);

    }
#   endif

#   ifdef ANALYZE_SURFACE
    {
        f64_t w = 0.0;
        f64_t dw = FP64_PI / 16.0;
        i32_t niter  = 0;

        do
        {
            #define ARRAY_SIZE (23*60*60*32)

            f64_t N = ARRAY_SIZE;

            // Function to load an array from a file
            static complex64_t x[ARRAY_SIZE];
            static complex64_t y[ARRAY_SIZE];

            for( u32_t i = 0; i < ARRAY_SIZE; i++ )
            {
                f64_t t = FP64_PI * ( i / (ARRAY_SIZE - 1.) );

                x[ i ]  = (complex64_t) { .re = cos( t     ), .im = 0.0, };
                y[ i ]  = (complex64_t) { .re = cos( t + w ), .im = 0.0, };
            }

            f64_t surface = math_surface( x, y, ARRAY_SIZE );

            printf("%d, %f, %.15e\n", niter, w, surface);
            niter++;
            w = niter * dw;

        }
        while( w <= 2.*FP64_PI );
    }

    {
        f64_t k = 0.0;
        f64_t dk = 0.03125;
        i32_t niter  = 0;

        do
        {
//            #define ARRAY_SIZE (1*60*60*32)

            f64_t N = ARRAY_SIZE;

            // Function to load an array from a file
            static complex64_t x[ARRAY_SIZE];
            static complex64_t y[ARRAY_SIZE];

            for( u32_t i = 0; i < ARRAY_SIZE; i++ )
            {
                f64_t t = FP64_PI * ( i / (ARRAY_SIZE - 1.) );

                x[ i ]  = (complex64_t) { .re =     cos( t ), .im =     sin( t ), };
                y[ i ]  = (complex64_t) { .re = k * cos( t ), .im = k * sin( t ), };
            }

            f64_t surface, report_elapsed;
            {
                struct timeval t1, t2;
                double elapsedTime;

                // start timer
                mingw_gettimeofday(&t1, NULL);

                surface = math_surface( x, y, ARRAY_SIZE );

                // stop timer
                mingw_gettimeofday(&t2, NULL);

                // compute and print the elapsed time in millisec
                elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0;      // sec to ms
                elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;   // us to ms

                report_elapsed = elapsedTime;
            }

            f64_t expect = 1.-fabs(1.-k)/fabs(1.+k);

            printf("%+5.f, %.15e, %.15e, %f\n", report_elapsed, expect, surface, fp64_ulps(surface, expect));
            niter++;
            k += dk;
        }
        while( k <= 5. );
    }

    {
        do
        {
//            #define ARRAY_SIZE (64 /*1*60*60*32*/)

            f64_t N = ARRAY_SIZE;

            // Function to load an array from a file
            static complex64_t x[ARRAY_SIZE];
            static complex64_t y[ARRAY_SIZE];

            for( u32_t i = 0; i < ARRAY_SIZE; i++ )
            {
                f64_t t = FP64_PI * ( i / (ARRAY_SIZE - 1.) );

                x[ i ]  = (complex64_t) { .re = 1.0, .im = 0.0, };
                y[ i ]  = (complex64_t) { .re = 0.0, .im = 1.0, };
            }

            f64_t surface, report_elapsed;
            {
                struct timeval t1, t2;
                double elapsedTime;

                // start timer
                mingw_gettimeofday(&t1, NULL);

                surface = math_surface( x, y, ARRAY_SIZE );

                // stop timer
                mingw_gettimeofday(&t2, NULL);

                // compute and print the elapsed time in millisec
                elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0;      // sec to ms
                elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;   // us to ms

                report_elapsed = elapsedTime;
            }

            f64_t expect = 0.29289321881345247559915563789515;

            printf("%+5.f, %.15e, %.15e, %f\n", report_elapsed, expect, surface, fp64_ulps(surface, expect));
        }
        while( 0 );
    }

    {
        do
        {
//            #define ARRAY_SIZE (23*60*60*32)

            f64_t N = ARRAY_SIZE;

            // Function to load an array from a file
            static complex64_t x[ARRAY_SIZE];
            static complex64_t y[ARRAY_SIZE];

            for( u32_t i = 0; i < ARRAY_SIZE; i++ )
            {
                x[ i ]  = (complex64_t) { .re = 1.0, .im = 1.0, };
                y[ i ]  = (complex64_t) { .re = 1.0, .im = 1.0, };
            }

            f64_t surface, report_elapsed;
            {
                struct timeval t1, t2;
                double elapsedTime;

                // start timer
                mingw_gettimeofday(&t1, NULL);

                surface = math_surface( x, y, ARRAY_SIZE );

                // stop timer
                mingw_gettimeofday(&t2, NULL);

                // compute and print the elapsed time in millisec
                elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0;      // sec to ms
                elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;   // us to ms

                report_elapsed = elapsedTime;
            }

            f64_t expect = 1.0;

            printf("%+5.f, %.15e, %.15e, %f\n", report_elapsed, expect, surface, fp64_ulps(surface, expect));
        }
        while( 0 );
    }
#   endif

#   ifdef ANALYZE_COMPLEX
{
    const f64_t eps = fp64_get_named_fp_in_real_line( NAMED_FP_MINNORM );

    complex64_t x = { .re = 3./4., .im = 3./4.*(1. - 4.* eps ) };
    complex64_t y = { .re = 2./3.*(1.+11.*eps ), .im = 2./3.*(1. + 5.*eps) };

    printf("X: %.18e + %.18e\n", x.re, x.im );
    printf("Y: %.18e + %.18e\n", y.re, y.im );

    complex64_t z = { .re = x.re*y.re - x.im*y.im, .im =x.re*y.im + x.im*y.re };

    printf("Z: %.15e + %.15e\n", z.re, z.im );


    complex64_t x0 = { .re = 3./4., .im = 3./4. };
    complex64_t y0 = { .re = 2./3., .im = 2./3. };

    printf("X: %.18e + %.18e\n", x0.re, x0.im );
    printf("Y: %.18e + %.18e\n", y0.re, y0.im );

    complex64_t z0 = { .re = x0.re*y0.re - x0.im*y0.im, .im =x0.re*y0.im + x0.im*y0.re };

    printf("Z: %.15e + %.15e\n", z0.re, z0.im );

}
#   endif

    return 0;
}

#define HUGE_EXP 1.e100
static f64_t m_y = 0.0;
f64_t math_powy( f64_t x ) { return ( ((-HUGE_EXP < x) && (x < +HUGE_EXP)) ? math_pow( x, m_y ) : 0./0. ); }
f64_t      powy( f64_t x ) { return ( ((-HUGE_EXP < x) && (x < +HUGE_EXP)) ?      pow( x, m_y ) : 0./0. ); }

f64_t      sinx( f64_t x ) { return ( ((-MATH_SIN_RANGE < x) && (x < +MATH_SIN_RANGE)) ? sin(x) : 0./0. ); }
f64_t      cosx( f64_t x ) { return ( ((-MATH_SIN_RANGE < x) && (x < +MATH_SIN_RANGE)) ? cos(x) : 0./0. ); }

static int load_array_from_file(const char *filename, double *array, int size)
{
    FILE *file = fopen(filename, "r");

    if (!file)
    {
        perror("Error opening file");
        return -1;
    }

    int i = 0;
    while (i < size && fscanf(file, "%lf,", &array[i]) == 1)
    {
        i++;
    }

    fclose(file);

    if (i != size)
    {
        fprintf(stderr, "Warning: Expected %d elements, but read %d.\n", size, i);
        return -1;
    }

    return 0;
}

// Function to read the Pearson correlation coefficient from a file

double read_pearson_correlation(const char *filename)
{
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Error opening file");
        return -1.0; // Return an error value
    }

    double value;
    if (fscanf(file, "%lf", &value) != 1) {
        fprintf(stderr, "Error reading the Pearson correlation coefficient.\n");
        fclose(file);
        return -1.0; // Return an error value
    }

    fclose(file);
    return value;
}

void write_xy(const char *filename, const double * x, const double * y, unsigned n)
{
    FILE *file = fopen(filename, "w");
    if (!file) {
        perror("Error opening file");
        return; // Return an error value
    }

    for(unsigned i = 0; i < n; i++)
    {
        fprintf(file,"%+.23e, %+.23e\n", x[i], y[i] );
    }

    fclose(file);
}

void test_math_cwsetexp( f64_t (*tested_cwsetexp)(f64_t, i16_t), bool_t active )
{
    if( active )
    {
        #define expect_cwsetexp(_f_,_n_) math_cwnormalize((_f_), FINITE).f.d * pow(2.0, (_n_))

        /* test math_cwsetexp( f, n ), with 0.5 <= f < 1, and -1021 <= n <= +1023 */
        for( int n = -1021; n <= 1023; n++ )
        for( double f =-10.; f < 10.; f += 0.1 )
        TEST_ASSERT_EQUAL_DOUBLE( expect_cwsetexp(f, n), tested_cwsetexp(f, n) );

        /* test some specific points */
        TEST_ASSERT_EQUAL_DOUBLE( expect_cwsetexp(+0.5, 2), tested_cwsetexp(+0.5, 2) );
        TEST_ASSERT_EQUAL_DOUBLE( expect_cwsetexp(-0.5, 2), tested_cwsetexp(-0.5, 2) );

        TEST_ASSERT_EQUAL_DOUBLE( expect_cwsetexp(+0.0, 2), tested_cwsetexp(+0.0, 2) );
        TEST_ASSERT_EQUAL_DOUBLE( expect_cwsetexp(-0.0, 2), tested_cwsetexp(-0.0, 2) );

        #undef expect_cwsetexp
    }
}
