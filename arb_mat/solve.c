/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"

void
_arb_mat_get_mid(arb_mat_t B, const arb_mat_t A)
{
    slong i, j;

    for (i = 0; i < arb_mat_nrows(A); i++)
        for (j = 0; j < arb_mat_ncols(A); j++)
            arb_get_mid_arb(arb_mat_entry(B, i, j), arb_mat_entry(A, i, j));
}

int
_arb_mat_solve_hansen_method_1(arb_mat_t X, const arb_mat_t A, const arb_mat_t B, slong prec)
{
    int result;
    slong n, m, *perm;
    arb_mat_t LU;

    n = arb_mat_nrows(A);
    m = arb_mat_ncols(X);

    if (n == 0 || m == 0)
        return 1;

    perm = _perm_init(n);
    arb_mat_init(LU, n, n);

    result = arb_mat_lu(perm, LU, A, prec);

    if (result)
        arb_mat_solve_lu_precomp(X, perm, LU, B, prec);

    arb_mat_clear(LU);
    _perm_clear(perm);

    return result;
}

int
_arb_mat_solve_hansen_method_4(arb_mat_t X, const arb_mat_t A, const arb_mat_t B, slong prec)
{
    int result;
    slong m, n;
    arb_mat_t I, C;

    n = arb_mat_nrows(A);
    m = arb_mat_ncols(X);

    if (n == 0 || m == 0)
        return 1;

    arb_mat_init(I, n, n);
    arb_mat_init(C, n, n);

    arb_mat_one(I);
    result = _arb_mat_solve_hansen_method_1(C, A, I, prec);

    if (result)
    {
        arb_mat_t CA, CB;

        _arb_mat_get_mid(C, C);

        arb_mat_init(CA, n, n);
        arb_mat_init(CB, n, m);
        arb_mat_mul(CA, C, A, prec);
        arb_mat_mul(CB, C, B, prec);

        result = _arb_mat_solve_hansen_method_1(X, CA, CB, prec);

        arb_mat_clear(CA);
        arb_mat_clear(CB);
    }

    arb_mat_clear(I);
    arb_mat_clear(C);

    return result;
}

int
arb_mat_solve(arb_mat_t X, const arb_mat_t A, const arb_mat_t B, slong prec)
{
    return _arb_mat_solve_hansen_method_4(X, A, B, prec);
}
