#include <assert.h>
#include <stdio.h>
#include "../src/private_math.h"

#define INTRND_ASSERT_MSG "INTRND failed\n"

short _Dint_v0(double *px, short xexp)
{	/* test and drop (scaled) fraction bits */
	unsigned short *ps = (unsigned short *)px;

	unsigned short frac = ( ps[ W0] & DFRAC ) || ps[ W1] || ps[ W2] || ps[ W3];

	short exponent = ( ps[ W0 ] & DMASK ) >> DOFF;

	if ( exponent == 0 && !frac )
	{
		return (ZERO);	/* zero */
	}
	else if ( exponent != DMAX )
	{
		;	/* finite */
	}
	else if (!frac)
	{
		return (INF);
	}
	else
	{	/* NaN */
		return (NIL);
	}

	exponent = (DBIAS+48+DOFF+1) - exponent - xexp;
	if ( exponent <= 0 )
	{
		return (0);	/* no frac bits to drop */
	}
	else if ( (48+DOFF) < exponent )
	{	/* all frac bits */
		ps[W0] = 0, ps[W1] = 0;
		ps[W2] = 0, ps[W3] = 0;
		return (FINITE);
	}
	else
	{	/* strip out frac bits */
		static const unsigned short mask[] = {
			0x0000, 0x0001, 0x0003, 0x0007,
			0x000f, 0x001f, 0x003f, 0x007f,
			0x00ff, 0x01ff, 0x03ff, 0x07ff,
			0x0fff, 0x1fff, 0x3fff, 0x7fff};
		static const u64_t sub[] = {W3, W2, W1, W0};

        if(*px<0.0)
        {
            *px-=0.5;
        }
        else
        {
            *px+=0.5;
        }

		frac = mask[exponent & 0xf];
		exponent >>= 4;
		frac &= ps[sub[exponent]];
		ps[sub[exponent]] ^= frac;
		switch (exponent)
		{	/* cascade through! */
			case 3:
				frac |= ps[W1], ps[W1] = 0;
			case 2:
				frac |= ps[W2], ps[W2] = 0;
			case 1:
				frac |= ps[W3], ps[W3] = 0;
		}

		return (frac ? FINITE : 0);
	}
}

double math_intrnd( double x )
{
	if(x < 0.0)
    {
        x-=0.5;
    }
    else
    {
        x+=0.5;
    }

    unsigned short *ps = (unsigned short *)&x;

	unsigned short frac = ( ps[ W0] & DFRAC ) || ps[ W1] || ps[ W2] || ps[ W3];

	short exponent = ( ps[ W0 ] & DMASK ) >> DOFF;

	if ( exponent == 0 && !frac )
	{
		return (x);	/* zero */
	}
	else if ( exponent != DMAX )
	{
		;	/* finite */
	}
	else if (!frac)
	{
		return (x);
	}
	else
	{	/* NaN */
		return (x);
	}

	exponent = (DBIAS+48+DOFF+1) - exponent - 0;
	if ( exponent <= 0 )
	{
		return (0);	/* no frac bits to drop */
	}
	else if ( (48+DOFF) < exponent )
	{	/* all frac bits */
		ps[W0] = 0, ps[W1] = 0;
		ps[W2] = 0, ps[W3] = 0;
		return (x);
	}
	else
	{	/* strip out frac bits */
		static const unsigned short mask[] = {
			0x0000, 0x0001, 0x0003, 0x0007,
			0x000f, 0x001f, 0x003f, 0x007f,
			0x00ff, 0x01ff, 0x03ff, 0x07ff,
			0x0fff, 0x1fff, 0x3fff, 0x7fff};
		static const u64_t sub[] = {W3, W2, W1, W0};

		frac = mask[exponent & 0xf];
		exponent >>= 4;
		frac &= ps[sub[exponent]];
		ps[sub[exponent]] ^= frac;
		switch (exponent)
		{	/* cascade through! */
			case 3:
				frac |= ps[W1], ps[W1] = 0;
			case 2:
				frac |= ps[W2], ps[W2] = 0;
			case 1:
				frac |= ps[W3], ps[W3] = 0;
		}

		return (frac ? x : 0.0);
	}
}

double _Dint( dw_t x, short xexp)
{
	_Bool frac = ( x.w[ W0] & DFRAC ) != 0 || x.w[ W1 ] != 0 || x.w[ W2 ] != 0 || x.w[ W3 ] != 0;

	short exponent = ( x.w[ W0 ] & DMASK ) >> DOFF;

    u16_t dtype = math_type( x.d );

    if( dtype == ZERO || dtype == INF || dtype == NIL )
    {
        return x.d;
    }

    if(x.d < 0.0)
    {
        x.d-=0.5;
    }
    else
    {
        x.d+=0.5;
    }

	exponent = (DBIAS+48+DOFF+1) - exponent - xexp;

	if ( exponent <= 0 )
	{
		return (0.0);	/* no frac bits to drop */
	}

	if ( (48+DOFF) < exponent )
	{	/* all frac bits */
		return (0.0);
	}

    static const unsigned short mask[] = {
        0x0000, 0x0001, 0x0003, 0x0007,
        0x000f, 0x001f, 0x003f, 0x007f,
        0x00ff, 0x01ff, 0x03ff, 0x07ff,
        0x0fff, 0x1fff, 0x3fff, 0x7fff};

    static const u64_t sub[] = {W3, W2, W1, W0};

    frac = mask[exponent & 0xf];

    exponent >>= 4;

    frac &= x.w[sub[exponent]];
    x.w[sub[exponent]] ^= frac;

    switch (exponent)
    {	/* cascade through! */
        case 3:
            frac |= x.w[W1], x.w[W1] = 0;
        case 2:
            frac |= x.w[W2], x.w[W2] = 0;
        case 1:
            frac |= x.w[W3], x.w[W3] = 0;
    }

    return (frac ? x.d : 0.0);
}

void test_intrnd( void )
{
    assert( +0.0 == _Dint((dw_t){ .d=-0.1   }, 0) && INTRND_ASSERT_MSG );
    assert( +0.0 == _Dint((dw_t){ .d=+0.1   }, 0) && INTRND_ASSERT_MSG );
    assert( +1.5 == _Dint((dw_t){ .d=+1.1   }, 0) && INTRND_ASSERT_MSG );
    //assert( +2.0 == _Dint((dw_t){ .d=+2.001 }, 0) && INTRND_ASSERT_MSG );
    assert( -2.0 == _Dint((dw_t){ .d=-1.57  }, 0) && INTRND_ASSERT_MSG );
    assert( +2.0 == _Dint((dw_t){ .d=+1.57  }, 0) && INTRND_ASSERT_MSG );
}

int main( void )
{
    test_intrnd();

    double          x = 10.;
    double          C = 1.5707963267948966192313216916398;

    double          k = _Dint( (dw_t){.d=x/C}, 0);

    double          g = x - (k*C);

    double          e = 0.550618933926680444456560297951;

    printf("%.17e, %.17e, %.17f\n", g, e, k);

    return ( 0 );
}