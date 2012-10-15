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

    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "fmprb_poly.h"

void
_fmprb_poly_mullow_classical(fmprb_struct * res,
    const fmprb_struct * poly1, long len1,
    const fmprb_struct * poly2, long len2, long n, long prec)
{
    if ((len1 == 1 && len2 == 1) || n == 1)
    {
        fmprb_mul(res, poly1, poly2, prec);
    }
    else /* Ordinary case */
    {
        long i;

        /* Set res[i] = poly1[i]*poly2[0] */
        _fmprb_vec_scalar_mul(res, poly1, FLINT_MIN(len1, n), poly2, prec);

        /* Set res[i+len1-1] = in1[len1-1]*in2[i] */
        if (n > len1)
            _fmprb_vec_scalar_mul(res + len1, poly2 + 1, n - len1,
                                      poly1 + len1 - 1, prec);

        /* out[i+j] += in1[i]*in2[j] */
        for (i = 0; i < FLINT_MIN(len1, n) - 1; i++)
            _fmprb_vec_scalar_addmul(res + i + 1, poly2 + 1,
                                         FLINT_MIN(len2, n - i) - 1,
                                         poly1 + i, prec);
    }
}

void
fmprb_poly_mullow_classical(fmprb_poly_t res, const fmprb_poly_t poly1,
                                            const fmprb_poly_t poly2,
                                                long n, long prec)
{
    long len_out;

    if (poly1->length == 0 || poly2->length == 0 || n == 0)
    {
        fmprb_poly_zero(res);
        return;
    }

    len_out = poly1->length + poly2->length - 1;
    if (n > len_out)
        n = len_out;

    if (res == poly1 || res == poly2)
    {
        fmprb_poly_t t;
        fmprb_poly_init2(t, n);
        _fmprb_poly_mullow_classical(t->coeffs, poly1->coeffs, poly1->length,
                                    poly2->coeffs, poly2->length, n, prec);
        fmprb_poly_swap(res, t);
        fmprb_poly_clear(t);
    }
    else
    {
        fmprb_poly_fit_length(res, n);
        _fmprb_poly_mullow_classical(res->coeffs, poly1->coeffs, poly1->length,
                                    poly2->coeffs, poly2->length, n, prec);
    }

    _fmprb_poly_set_length(res, n);
    _fmprb_poly_normalise(res);
}