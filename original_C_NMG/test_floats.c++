#include <cstdio>
#include <float.h>
#include <math.h>
#include <stdio.h>

void TestFunction()
{
    float_t num(1.0f);
    num.i -= 1;
    printf("Float value, representation, sign, exponent, mantissa\n");
    for (;;)
    {
        // Breakpoint here.
        printf("%1.8e, 0x%08X, %d, %d, 0x%06X\n",
               num.f, num.i,
               num.parts.sign, num.parts.exponent, num.parts.mantissa);
    }
}