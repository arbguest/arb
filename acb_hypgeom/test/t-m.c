/*=============================================================================

    This file is part of ARB.

    ARB is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    ARB is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with ARB; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2015 Fredrik Johansson

******************************************************************************/

#include "acb_hypgeom.h"

void
acb_randtest_maybe_half_int(acb_t x, flint_rand_t state, long prec, long size)
{
    if (n_randint(state, 8) == 0)
    {
        fmpz_t t;
        fmpz_init(t);
        fmpz_randtest(t, state, 1 + n_randint(state, prec));
        arb_set_fmpz(acb_realref(x), t);
        arb_zero(acb_imagref(x));
        acb_mul_2exp_si(x, x, -1);
        fmpz_clear(t);
    }
    else
    {
        acb_randtest(x, state, prec, size);
    }
}

int main()
{
    long iter;
    flint_rand_t state;

    printf("m....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 2000; iter++)
    {
        acb_t a0, a1, a2, b, z, w0, w1, w2, t, u;
        long prec0, prec1, prec2;
        int regularized, ebits;

        acb_init(a0);
        acb_init(a1);
        acb_init(a2);
        acb_init(b);
        acb_init(z);
        acb_init(w0);
        acb_init(w1);
        acb_init(w2);
        acb_init(t);
        acb_init(u);

        prec0 = 2 + n_randint(state, 700);
        prec1 = 2 + n_randint(state, 700);
        prec2 = 2 + n_randint(state, 700);

        if (n_randint(state, 5) == 0)
            ebits = 100;
        else
            ebits = 10;

        acb_randtest_maybe_half_int(a0, state, 1 + n_randint(state, 1000), 1 + n_randint(state, ebits));
        acb_randtest_maybe_half_int(b, state, 1 + n_randint(state, 1000), 1 + n_randint(state, ebits));
        acb_randtest(z, state, 1 + n_randint(state, 1000), 1 + n_randint(state, ebits));
        acb_randtest(w0, state, 1 + n_randint(state, 1000), 1 + n_randint(state, ebits));
        acb_randtest(w1, state, 1 + n_randint(state, 1000), 1 + n_randint(state, ebits));
        acb_randtest(w2, state, 1 + n_randint(state, 1000), 1 + n_randint(state, ebits));
        regularized = n_randint(state, 2);

        acb_add_ui(a1, a0, 1, prec0);
        acb_add_ui(a2, a0, 2, prec0);

        switch (n_randint(state, 3))
        {
            case 0:
                acb_hypgeom_m_asymp(w0, a0, b, z, regularized, prec0);
                break;
            case 1:
                acb_hypgeom_m_1f1(w0, a0, b, z, regularized, prec0);
                break;
            default:
                acb_hypgeom_m(w0, a0, b, z, regularized, prec0);
        }

        switch (n_randint(state, 3))
        {
            case 0:
                acb_hypgeom_m_asymp(w1, a0, b, z, regularized, prec1);
                break;
            case 1:
                acb_hypgeom_m_1f1(w1, a0, b, z, regularized, prec1);
                break;
            default:
                acb_hypgeom_m(w1, a0, b, z, regularized, prec1);
        }

        if (!acb_overlaps(w0, w1))
        {
            printf("FAIL: consistency\n\n");
            printf("regularized = %d\n\n", regularized);
            printf("a = "); acb_printd(a0, 30); printf("\n\n");
            printf("b = "); acb_printd(b, 30); printf("\n\n");
            printf("z = "); acb_printd(z, 30); printf("\n\n");
            printf("w0 = "); acb_printd(w0, 30); printf("\n\n");
            printf("w1 = "); acb_printd(w1, 30); printf("\n\n");
            abort();
        }

        switch (n_randint(state, 3))
        {
            case 0:
                acb_hypgeom_m_asymp(w1, a1, b, z, regularized, prec1);
                break;
            case 1:
                acb_hypgeom_m_1f1(w1, a1, b, z, regularized, prec1);
                break;
            default:
                acb_hypgeom_m(w1, a1, b, z, regularized, prec1);
        }

        switch (n_randint(state, 3))
        {
            case 0:
                acb_hypgeom_m_asymp(w2, a2, b, z, regularized, prec2);
                break;
            case 1:
                acb_hypgeom_m_1f1(w2, a2, b, z, regularized, prec2);
                break;
            default:
                acb_hypgeom_m(w2, a2, b, z, regularized, prec2);
        }

        /* 1F1(a,b,z) * (a-b+1) */
        acb_sub(u, a1, b, prec0);
        acb_mul(t, w0, u, prec0);

        /* 1F1(a+1,b,z) * (2a - b + z + 2) */
        acb_mul_2exp_si(u, a0, 1);
        acb_sub(u, u, b, prec0);
        acb_add(u, u, z, prec0);
        acb_add_ui(u, u, 2, prec0);
        acb_submul(t, w1, u, prec0);

        /* 1F1(a+2,b,z) * -(a+1) */
        acb_addmul(t, w2, a1, prec0);

        if (!acb_contains_zero(t))
        {
            printf("FAIL: contiguous relation\n\n");
            printf("a = "); acb_printd(a0, 30); printf("\n\n");
            printf("b = "); acb_printd(b, 30); printf("\n\n");
            printf("z = ");  acb_printd(z, 30); printf("\n\n");
            printf("w0 = "); acb_printd(w0, 30); printf("\n\n");
            printf("w1 = "); acb_printd(w1, 30); printf("\n\n");
            printf("w2 = "); acb_printd(w2, 30); printf("\n\n");
            printf("t = "); acb_printd(t, 30); printf("\n\n");
            abort();
        }

        acb_add_ui(a1, b, 1, prec0);
        acb_add_ui(a2, b, 2, prec0);

        switch (n_randint(state, 3))
        {
            case 0:
                acb_hypgeom_m_asymp(w1, a0, a1, z, regularized, prec1);
                break;
            case 1:
                acb_hypgeom_m_1f1(w1, a0, a1, z, regularized, prec1);
                break;
            default:
                acb_hypgeom_m(w1, a0, a1, z, regularized, prec1);
        }

        switch (n_randint(state, 3))
        {
            case 0:
                acb_hypgeom_m_asymp(w2, a0, a2, z, regularized, prec2);
                break;
            case 1:
                acb_hypgeom_m_1f1(w2, a0, a2, z, regularized, prec2);
                break;
            default:
                acb_hypgeom_m(w2, a0, a2, z, regularized, prec2);
        }

        if (regularized)
        {
            acb_set(t, w0);

            acb_add(u, b, z, prec0);
            acb_submul(t, w1, u, prec0);

            acb_sub(u, a0, a1, prec0);
            acb_mul(u, u, z, prec0);
            acb_submul(t, w2, u, prec0);
        }
        else
        {
            acb_mul(t, w0, b, prec0);
            acb_mul(t, t, a1, prec0);

            acb_add(u, b, z, prec0);
            acb_mul(u, u, a1, prec0);
            acb_submul(t, w1, u, prec0);

            acb_sub(u, a0, a1, prec0);
            acb_mul(u, u, z, prec0);
            acb_submul(t, w2, u, prec0);
        }

        if (!acb_contains_zero(t))
        {
            printf("FAIL: contiguous relation 2\n\n");
            printf("a = "); acb_printd(a0, 30); printf("\n\n");
            printf("b = "); acb_printd(b, 30); printf("\n\n");
            printf("z = ");  acb_printd(z, 30); printf("\n\n");
            printf("w0 = "); acb_printd(w0, 30); printf("\n\n");
            printf("w1 = "); acb_printd(w1, 30); printf("\n\n");
            printf("w2 = "); acb_printd(w2, 30); printf("\n\n");
            printf("t = "); acb_printd(t, 30); printf("\n\n");
            abort();
        }

        switch (n_randint(state, 3))
        {
            case 0:
                acb_hypgeom_m(a0, a0, b, z, regularized, prec0);
                acb_set(t, a0);
                break;
            case 1:
                acb_hypgeom_m(b, a0, b, z, regularized, prec0);
                acb_set(t, b);
                break;
            default:
                acb_hypgeom_m(z, a0, b, z, regularized, prec2);
                acb_set(t, z);
        }

        if (!acb_overlaps(t, w0))
        {
            printf("FAIL: aliasing\n\n");
            printf("a = "); acb_printd(a0, 30); printf("\n\n");
            printf("b = "); acb_printd(b, 30); printf("\n\n");
            printf("z = ");  acb_printd(z, 30); printf("\n\n");
            printf("w0 = "); acb_printd(w0, 30); printf("\n\n");
            printf("t = "); acb_printd(t, 30); printf("\n\n");
            abort();
        }

        acb_clear(a0);
        acb_clear(a1);
        acb_clear(a2);
        acb_clear(b);
        acb_clear(z);
        acb_clear(w0);
        acb_clear(w1);
        acb_clear(w2);
        acb_clear(t);
        acb_clear(u);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
