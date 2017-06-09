/*************************************************************************
ALGLIB 3.11.0 (source code generated 2017-05-11)
Copyright (c) Sergey Bochkanov (ALGLIB project).

>>> SOURCE LICENSE >>>
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation (www.fsf.org); either version 2 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License is available at
http://www.fsf.org/licensing/licenses
>>> END OF LICENSE >>>
*************************************************************************/
#include "stdafx.h"
#include "interpolation.h"

// disable some irrelevant warnings
#if (AE_COMPILER==AE_MSVC)
#pragma warning(disable:4100)
#pragma warning(disable:4127)
#pragma warning(disable:4702)
#pragma warning(disable:4996)
#endif
using namespace std;

/////////////////////////////////////////////////////////////////////////
//
// THIS SECTION CONTAINS IMPLEMENTATION OF C++ INTERFACE
//
/////////////////////////////////////////////////////////////////////////
namespace alglib
{
_spline1dinterpolant_owner::_spline1dinterpolant_owner()
{
    p_struct = (alglib_impl::spline1dinterpolant*)alglib_impl::ae_malloc(sizeof(alglib_impl::spline1dinterpolant), NULL);
    if( p_struct==NULL )
        throw ap_error("ALGLIB: malloc error");
    alglib_impl::_spline1dinterpolant_init(p_struct, NULL);
}

_spline1dinterpolant_owner::_spline1dinterpolant_owner(const _spline1dinterpolant_owner &rhs)
{
    p_struct = (alglib_impl::spline1dinterpolant*)alglib_impl::ae_malloc(sizeof(alglib_impl::spline1dinterpolant), NULL);
    if( p_struct==NULL )
        throw ap_error("ALGLIB: malloc error");
    alglib_impl::_spline1dinterpolant_init_copy(p_struct, const_cast<alglib_impl::spline1dinterpolant*>(rhs.p_struct), NULL);
}

_spline1dinterpolant_owner& _spline1dinterpolant_owner::operator=(const _spline1dinterpolant_owner &rhs)
{
    if( this==&rhs )
        return *this;
    alglib_impl::_spline1dinterpolant_clear(p_struct);
    alglib_impl::_spline1dinterpolant_init_copy(p_struct, const_cast<alglib_impl::spline1dinterpolant*>(rhs.p_struct), NULL);
    return *this;
}

_spline1dinterpolant_owner::~_spline1dinterpolant_owner()
{
    alglib_impl::_spline1dinterpolant_clear(p_struct);
    ae_free(p_struct);
}

alglib_impl::spline1dinterpolant* _spline1dinterpolant_owner::c_ptr()
{
    return p_struct;
}

alglib_impl::spline1dinterpolant* _spline1dinterpolant_owner::c_ptr() const
{
    return const_cast<alglib_impl::spline1dinterpolant*>(p_struct);
}
spline1dinterpolant::spline1dinterpolant() : _spline1dinterpolant_owner()
{
}

spline1dinterpolant::spline1dinterpolant(const spline1dinterpolant &rhs):_spline1dinterpolant_owner(rhs)
{
}

spline1dinterpolant& spline1dinterpolant::operator=(const spline1dinterpolant &rhs)
{
    if( this==&rhs )
        return *this;
    _spline1dinterpolant_owner::operator=(rhs);
    return *this;
}

spline1dinterpolant::~spline1dinterpolant()
{
}

void spline1dconvcubic(const real_1d_array &x, const real_1d_array &y, const real_1d_array &x2, real_1d_array &y2)
{
    alglib_impl::ae_state _alglib_env_state;
    ae_int_t n;
    ae_int_t boundltype;
    double boundl;
    ae_int_t boundrtype;
    double boundr;
    ae_int_t n2;
    if( (x.length()!=y.length()))
        throw ap_error("Error while calling 'spline1dconvcubic': looks like one of arguments has wrong size");
    n = x.length();
    boundltype = 0;
    boundl = 0;
    boundrtype = 0;
    boundr = 0;
    n2 = x2.length();
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::spline1dconvcubic(const_cast<alglib_impl::ae_vector*>(x.c_ptr()), const_cast<alglib_impl::ae_vector*>(y.c_ptr()), n, boundltype, boundl, boundrtype, boundr, const_cast<alglib_impl::ae_vector*>(x2.c_ptr()), n2, const_cast<alglib_impl::ae_vector*>(y2.c_ptr()), &_alglib_env_state);

        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
}


}

namespace alglib_impl
{
void _spline1dinterpolant_init(void* _p, ae_state *_state)
{
    spline1dinterpolant *p = (spline1dinterpolant*)_p;
    ae_touch_ptr((void*)p);
    ae_vector_init(&p->x, 0, DT_REAL, _state);
    ae_vector_init(&p->c, 0, DT_REAL, _state);
}


void _spline1dinterpolant_init_copy(void* _dst, void* _src, ae_state *_state)
{
    spline1dinterpolant *dst = (spline1dinterpolant*)_dst;
    spline1dinterpolant *src = (spline1dinterpolant*)_src;
    dst->periodic = src->periodic;
    dst->n = src->n;
    dst->k = src->k;
    dst->continuity = src->continuity;
    ae_vector_init_copy(&dst->x, &src->x, _state);
    ae_vector_init_copy(&dst->c, &src->c, _state);
}


void _spline1dinterpolant_clear(void* _p)
{
    spline1dinterpolant *p = (spline1dinterpolant*)_p;
    ae_touch_ptr((void*)p);
    ae_vector_clear(&p->x);
    ae_vector_clear(&p->c);
}


void _spline1dinterpolant_destroy(void* _p)
{
    spline1dinterpolant *p = (spline1dinterpolant*)_p;
    ae_touch_ptr((void*)p);
    ae_vector_destroy(&p->x);
    ae_vector_destroy(&p->c);
}
static void spline1d_spline1dgriddiffcubicinternal(/* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_int_t n,
     ae_int_t boundltype,
     double boundl,
     ae_int_t boundrtype,
     double boundr,
     /* Real    */ ae_vector* d,
     /* Real    */ ae_vector* a1,
     /* Real    */ ae_vector* a2,
     /* Real    */ ae_vector* a3,
     /* Real    */ ae_vector* b,
     /* Real    */ ae_vector* dt,
     ae_state *_state);
static void spline1d_heapsortppoints(/* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     /* Integer */ ae_vector* p,
     ae_int_t n,
     ae_state *_state);
static void spline1d_solvetridiagonal(/* Real    */ ae_vector* a,
     /* Real    */ ae_vector* b,
     /* Real    */ ae_vector* c,
     /* Real    */ ae_vector* d,
     ae_int_t n,
     /* Real    */ ae_vector* x,
     ae_state *_state);
static void spline1d_solvecyclictridiagonal(/* Real    */ ae_vector* a,
     /* Real    */ ae_vector* b,
     /* Real    */ ae_vector* c,
     /* Real    */ ae_vector* d,
     ae_int_t n,
     /* Real    */ ae_vector* x,
     ae_state *_state);
void spline1dconvcubic(/* Real    */ ae_vector* x,
                       /* Real    */ ae_vector* y,
                       ae_int_t n,
                       ae_int_t boundltype,
                       double boundl,
                       ae_int_t boundrtype,
                       double boundr,
                       /* Real    */ ae_vector* x2,
                       ae_int_t n2,
                       /* Real    */ ae_vector* y2,
                       ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector _x;
    ae_vector _y;
    ae_vector _x2;
    ae_vector a1;
    ae_vector a2;
    ae_vector a3;
    ae_vector b;
    ae_vector d;
    ae_vector dt;
    ae_vector d1;
    ae_vector d2;
    ae_vector p;
    ae_vector p2;
    ae_int_t i;
    ae_int_t ylen;
    double t;
    double t2;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init_copy(&_x, x, _state);
    x = &_x;
    ae_vector_init_copy(&_y, y, _state);
    y = &_y;
    ae_vector_init_copy(&_x2, x2, _state);
    x2 = &_x2;
    ae_vector_clear(y2);
    ae_vector_init(&a1, 0, DT_REAL, _state);
    ae_vector_init(&a2, 0, DT_REAL, _state);
    ae_vector_init(&a3, 0, DT_REAL, _state);
    ae_vector_init(&b, 0, DT_REAL, _state);
    ae_vector_init(&d, 0, DT_REAL, _state);
    ae_vector_init(&dt, 0, DT_REAL, _state);
    ae_vector_init(&d1, 0, DT_REAL, _state);
    ae_vector_init(&d2, 0, DT_REAL, _state);
    ae_vector_init(&p, 0, DT_INT, _state);
    ae_vector_init(&p2, 0, DT_INT, _state);


    /*
     * check correctness of boundary conditions
     */
    ae_assert(((boundltype==-1||boundltype==0)||boundltype==1)||boundltype==2, "Spline1DConvCubic: incorrect BoundLType!", _state);
    ae_assert(((boundrtype==-1||boundrtype==0)||boundrtype==1)||boundrtype==2, "Spline1DConvCubic: incorrect BoundRType!", _state);
    ae_assert((boundrtype==-1&&boundltype==-1)||(boundrtype!=-1&&boundltype!=-1), "Spline1DConvCubic: incorrect BoundLType/BoundRType!", _state);
    if( boundltype==1||boundltype==2 )
    {
        ae_assert(ae_isfinite(boundl, _state), "Spline1DConvCubic: BoundL is infinite or NAN!", _state);
    }
    if( boundrtype==1||boundrtype==2 )
    {
        ae_assert(ae_isfinite(boundr, _state), "Spline1DConvCubic: BoundR is infinite or NAN!", _state);
    }

    /*
     * check lengths of arguments
     */
    ae_assert(n>=2, "Spline1DConvCubic: N<2!", _state);
    ae_assert(x->cnt>=n, "Spline1DConvCubic: Length(X)<N!", _state);
    ae_assert(y->cnt>=n, "Spline1DConvCubic: Length(Y)<N!", _state);
    ae_assert(n2>=2, "Spline1DConvCubic: N2<2!", _state);
    ae_assert(x2->cnt>=n2, "Spline1DConvCubic: Length(X2)<N2!", _state);

    /*
     * check and sort X/Y
     */
    ylen = n;
    if( boundltype==-1 )
    {
        ylen = n-1;
    }
    ae_assert(isfinitevector(x, n, _state), "Spline1DConvCubic: X contains infinite or NAN values!", _state);
    ae_assert(isfinitevector(y, ylen, _state), "Spline1DConvCubic: Y contains infinite or NAN values!", _state);
    ae_assert(isfinitevector(x2, n2, _state), "Spline1DConvCubic: X2 contains infinite or NAN values!", _state);
    spline1d_heapsortppoints(x, y, &p, n, _state);
    ae_assert(aredistinct(x, n, _state), "Spline1DConvCubic: at least two consequent points are too close!", _state);

    /*
     * set up DT (we will need it below)
     */
    ae_vector_set_length(&dt, ae_maxint(n, n2, _state), _state);

    /*
     * sort X2:
     * * use fake array DT because HeapSortPPoints() needs both integer AND real arrays
     * * if we have periodic problem, wrap points
     * * sort them, store permutation at P2
     */
    if( boundrtype==-1&&boundltype==-1 )
    {
        for(i=0; i<=n2-1; i++)
        {
            t = x2->ptr.p_double[i];
            apperiodicmap(&t, x->ptr.p_double[0], x->ptr.p_double[n-1], &t2, _state);
            x2->ptr.p_double[i] = t;
        }
    }
    spline1d_heapsortppoints(x2, &dt, &p2, n2, _state);

    /*
     * Now we've checked and preordered everything, so we:
     * * call internal GridDiff() function to get Hermite form of spline
     * * convert using internal Conv() function
     * * convert Y2 back to original order
     */
    spline1d_spline1dgriddiffcubicinternal(x, y, n, boundltype, boundl, boundrtype, boundr, &d, &a1, &a2, &a3, &b, &dt, _state);
    spline1dconvdiffinternal(x, y, &d, n, x2, n2, y2, ae_true, &d1, ae_false, &d2, ae_false, _state);
    ae_assert(dt.cnt>=n2, "Spline1DConvCubic: internal error!", _state);
    for(i=0; i<=n2-1; i++)
    {
        dt.ptr.p_double[p2.ptr.p_int[i]] = y2->ptr.p_double[i];
    }
    ae_v_move(&y2->ptr.p_double[0], 1, &dt.ptr.p_double[0], 1, ae_v_len(0,n2-1));
    ae_frame_leave(_state);
}

static void spline1d_heapsortppoints(/* Real    */ ae_vector* x,
                                     /* Real    */ ae_vector* y,
                                     /* Integer */ ae_vector* p,
                                     ae_int_t n,
                                     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector rbuf;
    ae_vector ibuf;
    ae_int_t i;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init(&rbuf, 0, DT_REAL, _state);
    ae_vector_init(&ibuf, 0, DT_INT, _state);

    if( p->cnt<n )
    {
        ae_vector_set_length(p, n, _state);
    }
    ae_vector_set_length(&rbuf, n, _state);
    for(i=0; i<=n-1; i++)
    {
        p->ptr.p_int[i] = i;
    }
    tagsortfasti(x, p, &rbuf, &ibuf, n, _state);
    for(i=0; i<=n-1; i++)
    {
        rbuf.ptr.p_double[i] = y->ptr.p_double[p->ptr.p_int[i]];
    }
    ae_v_move(&y->ptr.p_double[0], 1, &rbuf.ptr.p_double[0], 1, ae_v_len(0,n-1));
    ae_frame_leave(_state);
}
static void spline1d_spline1dgriddiffcubicinternal(/* Real    */ ae_vector* x,
                                                   /* Real    */ ae_vector* y,
                                                   ae_int_t n,
                                                   ae_int_t boundltype,
                                                   double boundl,
                                                   ae_int_t boundrtype,
                                                   double boundr,
                                                   /* Real    */ ae_vector* d,
                                                   /* Real    */ ae_vector* a1,
                                                   /* Real    */ ae_vector* a2,
                                                   /* Real    */ ae_vector* a3,
                                                   /* Real    */ ae_vector* b,
                                                   /* Real    */ ae_vector* dt,
                                                   ae_state *_state)
{
    ae_int_t i;



    /*
     * allocate arrays
     */
    if( d->cnt<n )
    {
        ae_vector_set_length(d, n, _state);
    }
    if( a1->cnt<n )
    {
        ae_vector_set_length(a1, n, _state);
    }
    if( a2->cnt<n )
    {
        ae_vector_set_length(a2, n, _state);
    }
    if( a3->cnt<n )
    {
        ae_vector_set_length(a3, n, _state);
    }
    if( b->cnt<n )
    {
        ae_vector_set_length(b, n, _state);
    }
    if( dt->cnt<n )
    {
        ae_vector_set_length(dt, n, _state);
    }

    /*
     * Special cases:
     * * N=2, parabolic terminated boundary condition on both ends
     * * N=2, periodic boundary condition
     */
    if( (n==2&&boundltype==0)&&boundrtype==0 )
    {
        d->ptr.p_double[0] = (y->ptr.p_double[1]-y->ptr.p_double[0])/(x->ptr.p_double[1]-x->ptr.p_double[0]);
        d->ptr.p_double[1] = d->ptr.p_double[0];
        return;
    }
    if( (n==2&&boundltype==-1)&&boundrtype==-1 )
    {
        d->ptr.p_double[0] = (double)(0);
        d->ptr.p_double[1] = (double)(0);
        return;
    }

    /*
     * Periodic and non-periodic boundary conditions are
     * two separate classes
     */
    if( boundrtype==-1&&boundltype==-1 )
    {

        /*
         * Periodic boundary conditions
         */
        y->ptr.p_double[n-1] = y->ptr.p_double[0];

        /*
         * Boundary conditions at N-1 points
         * (one point less because last point is the same as first point).
         */
        a1->ptr.p_double[0] = x->ptr.p_double[1]-x->ptr.p_double[0];
        a2->ptr.p_double[0] = 2*(x->ptr.p_double[1]-x->ptr.p_double[0]+x->ptr.p_double[n-1]-x->ptr.p_double[n-2]);
        a3->ptr.p_double[0] = x->ptr.p_double[n-1]-x->ptr.p_double[n-2];
        b->ptr.p_double[0] = 3*(y->ptr.p_double[n-1]-y->ptr.p_double[n-2])/(x->ptr.p_double[n-1]-x->ptr.p_double[n-2])*(x->ptr.p_double[1]-x->ptr.p_double[0])+3*(y->ptr.p_double[1]-y->ptr.p_double[0])/(x->ptr.p_double[1]-x->ptr.p_double[0])*(x->ptr.p_double[n-1]-x->ptr.p_double[n-2]);
        for(i=1; i<=n-2; i++)
        {

            /*
             * Altough last point is [N-2], we use X[N-1] and Y[N-1]
             * (because of periodicity)
             */
            a1->ptr.p_double[i] = x->ptr.p_double[i+1]-x->ptr.p_double[i];
            a2->ptr.p_double[i] = 2*(x->ptr.p_double[i+1]-x->ptr.p_double[i-1]);
            a3->ptr.p_double[i] = x->ptr.p_double[i]-x->ptr.p_double[i-1];
            b->ptr.p_double[i] = 3*(y->ptr.p_double[i]-y->ptr.p_double[i-1])/(x->ptr.p_double[i]-x->ptr.p_double[i-1])*(x->ptr.p_double[i+1]-x->ptr.p_double[i])+3*(y->ptr.p_double[i+1]-y->ptr.p_double[i])/(x->ptr.p_double[i+1]-x->ptr.p_double[i])*(x->ptr.p_double[i]-x->ptr.p_double[i-1]);
        }

        /*
         * Solve, add last point (with index N-1)
         */
        spline1d_solvecyclictridiagonal(a1, a2, a3, b, n-1, dt, _state);
        ae_v_move(&d->ptr.p_double[0], 1, &dt->ptr.p_double[0], 1, ae_v_len(0,n-2));
        d->ptr.p_double[n-1] = d->ptr.p_double[0];
    }
    else
    {

        /*
         * Non-periodic boundary condition.
         * Left boundary conditions.
         */
        if( boundltype==0 )
        {
            a1->ptr.p_double[0] = (double)(0);
            a2->ptr.p_double[0] = (double)(1);
            a3->ptr.p_double[0] = (double)(1);
            b->ptr.p_double[0] = 2*(y->ptr.p_double[1]-y->ptr.p_double[0])/(x->ptr.p_double[1]-x->ptr.p_double[0]);
        }
        if( boundltype==1 )
        {
            a1->ptr.p_double[0] = (double)(0);
            a2->ptr.p_double[0] = (double)(1);
            a3->ptr.p_double[0] = (double)(0);
            b->ptr.p_double[0] = boundl;
        }
        if( boundltype==2 )
        {
            a1->ptr.p_double[0] = (double)(0);
            a2->ptr.p_double[0] = (double)(2);
            a3->ptr.p_double[0] = (double)(1);
            b->ptr.p_double[0] = 3*(y->ptr.p_double[1]-y->ptr.p_double[0])/(x->ptr.p_double[1]-x->ptr.p_double[0])-0.5*boundl*(x->ptr.p_double[1]-x->ptr.p_double[0]);
        }

        /*
         * Central conditions
         */
        for(i=1; i<=n-2; i++)
        {
            a1->ptr.p_double[i] = x->ptr.p_double[i+1]-x->ptr.p_double[i];
            a2->ptr.p_double[i] = 2*(x->ptr.p_double[i+1]-x->ptr.p_double[i-1]);
            a3->ptr.p_double[i] = x->ptr.p_double[i]-x->ptr.p_double[i-1];
            b->ptr.p_double[i] = 3*(y->ptr.p_double[i]-y->ptr.p_double[i-1])/(x->ptr.p_double[i]-x->ptr.p_double[i-1])*(x->ptr.p_double[i+1]-x->ptr.p_double[i])+3*(y->ptr.p_double[i+1]-y->ptr.p_double[i])/(x->ptr.p_double[i+1]-x->ptr.p_double[i])*(x->ptr.p_double[i]-x->ptr.p_double[i-1]);
        }

        /*
         * Right boundary conditions
         */
        if( boundrtype==0 )
        {
            a1->ptr.p_double[n-1] = (double)(1);
            a2->ptr.p_double[n-1] = (double)(1);
            a3->ptr.p_double[n-1] = (double)(0);
            b->ptr.p_double[n-1] = 2*(y->ptr.p_double[n-1]-y->ptr.p_double[n-2])/(x->ptr.p_double[n-1]-x->ptr.p_double[n-2]);
        }
        if( boundrtype==1 )
        {
            a1->ptr.p_double[n-1] = (double)(0);
            a2->ptr.p_double[n-1] = (double)(1);
            a3->ptr.p_double[n-1] = (double)(0);
            b->ptr.p_double[n-1] = boundr;
        }
        if( boundrtype==2 )
        {
            a1->ptr.p_double[n-1] = (double)(1);
            a2->ptr.p_double[n-1] = (double)(2);
            a3->ptr.p_double[n-1] = (double)(0);
            b->ptr.p_double[n-1] = 3*(y->ptr.p_double[n-1]-y->ptr.p_double[n-2])/(x->ptr.p_double[n-1]-x->ptr.p_double[n-2])+0.5*boundr*(x->ptr.p_double[n-1]-x->ptr.p_double[n-2]);
        }

        /*
         * Solve
         */
        spline1d_solvetridiagonal(a1, a2, a3, b, n, d, _state);
    }
}
static void spline1d_solvecyclictridiagonal(/* Real    */ ae_vector* a,
     /* Real    */ ae_vector* b,
     /* Real    */ ae_vector* c,
     /* Real    */ ae_vector* d,
     ae_int_t n,
     /* Real    */ ae_vector* x,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector _b;
    ae_int_t k;
    double alpha;
    double beta;
    double gamma;
    ae_vector y;
    ae_vector z;
    ae_vector u;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init_copy(&_b, b, _state);
    b = &_b;
    ae_vector_init(&y, 0, DT_REAL, _state);
    ae_vector_init(&z, 0, DT_REAL, _state);
    ae_vector_init(&u, 0, DT_REAL, _state);

    if( x->cnt<n )
    {
        ae_vector_set_length(x, n, _state);
    }
    beta = a->ptr.p_double[0];
    alpha = c->ptr.p_double[n-1];
    gamma = -b->ptr.p_double[0];
    b->ptr.p_double[0] = 2*b->ptr.p_double[0];
    b->ptr.p_double[n-1] = b->ptr.p_double[n-1]-alpha*beta/gamma;
    ae_vector_set_length(&u, n, _state);
    for(k=0; k<=n-1; k++)
    {
        u.ptr.p_double[k] = (double)(0);
    }
    u.ptr.p_double[0] = gamma;
    u.ptr.p_double[n-1] = alpha;
    spline1d_solvetridiagonal(a, b, c, d, n, &y, _state);
    spline1d_solvetridiagonal(a, b, c, &u, n, &z, _state);
    for(k=0; k<=n-1; k++)
    {
        x->ptr.p_double[k] = y.ptr.p_double[k]-(y.ptr.p_double[0]+beta/gamma*y.ptr.p_double[n-1])/(1+z.ptr.p_double[0]+beta/gamma*z.ptr.p_double[n-1])*z.ptr.p_double[k];
    }
    ae_frame_leave(_state);
}
static void spline1d_solvetridiagonal(/* Real    */ ae_vector* a,
                                      /* Real    */ ae_vector* b,
                                      /* Real    */ ae_vector* c,
                                      /* Real    */ ae_vector* d,
                                      ae_int_t n,
                                      /* Real    */ ae_vector* x,
                                      ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector _b;
    ae_vector _d;
    ae_int_t k;
    double t;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init_copy(&_b, b, _state);
    b = &_b;
    ae_vector_init_copy(&_d, d, _state);
    d = &_d;

    if( x->cnt<n )
    {
        ae_vector_set_length(x, n, _state);
    }
    for(k=1; k<=n-1; k++)
    {
        t = a->ptr.p_double[k]/b->ptr.p_double[k-1];
        b->ptr.p_double[k] = b->ptr.p_double[k]-t*c->ptr.p_double[k-1];
        d->ptr.p_double[k] = d->ptr.p_double[k]-t*d->ptr.p_double[k-1];
    }
    x->ptr.p_double[n-1] = d->ptr.p_double[n-1]/b->ptr.p_double[n-1];
    for(k=n-2; k>=0; k--)
    {
        x->ptr.p_double[k] = (d->ptr.p_double[k]-c->ptr.p_double[k]*x->ptr.p_double[k+1])/b->ptr.p_double[k];
    }
    ae_frame_leave(_state);
}
void spline1dconvdiffinternal(/* Real    */ ae_vector* xold,
                              /* Real    */ ae_vector* yold,
                              /* Real    */ ae_vector* dold,
                              ae_int_t n,
                              /* Real    */ ae_vector* x2,
                              ae_int_t n2,
                              /* Real    */ ae_vector* y,
                              ae_bool needy,
                              /* Real    */ ae_vector* d1,
                              ae_bool needd1,
                              /* Real    */ ae_vector* d2,
                              ae_bool needd2,
                              ae_state *_state)
{
    ae_int_t intervalindex;
    ae_int_t pointindex;
    ae_bool havetoadvance;
    double c0;
    double c1;
    double c2;
    double c3;
    double a;
    double b;
    double w;
    double w2;
    double w3;
    double fa;
    double fb;
    double da;
    double db;
    double t;



    /*
     * Prepare space
     */
    if( needy&&y->cnt<n2 )
    {
        ae_vector_set_length(y, n2, _state);
    }
    if( needd1&&d1->cnt<n2 )
    {
        ae_vector_set_length(d1, n2, _state);
    }
    if( needd2&&d2->cnt<n2 )
    {
        ae_vector_set_length(d2, n2, _state);
    }

    /*
     * These assignments aren't actually needed
     * (variables are initialized in the loop below),
     * but without them compiler will complain about uninitialized locals
     */
    c0 = (double)(0);
    c1 = (double)(0);
    c2 = (double)(0);
    c3 = (double)(0);
    a = (double)(0);
    b = (double)(0);

    /*
     * Cycle
     */
    intervalindex = -1;
    pointindex = 0;
    for(;;)
    {

        /*
         * are we ready to exit?
         */
        if( pointindex>=n2 )
        {
            break;
        }
        t = x2->ptr.p_double[pointindex];

        /*
         * do we need to advance interval?
         */
        havetoadvance = ae_false;
        if( intervalindex==-1 )
        {
            havetoadvance = ae_true;
        }
        else
        {
            if( intervalindex<n-2 )
            {
                havetoadvance = ae_fp_greater_eq(t,b);
            }
        }
        if( havetoadvance )
        {
            intervalindex = intervalindex+1;
            a = xold->ptr.p_double[intervalindex];
            b = xold->ptr.p_double[intervalindex+1];
            w = b-a;
            w2 = w*w;
            w3 = w*w2;
            fa = yold->ptr.p_double[intervalindex];
            fb = yold->ptr.p_double[intervalindex+1];
            da = dold->ptr.p_double[intervalindex];
            db = dold->ptr.p_double[intervalindex+1];
            c0 = fa;
            c1 = da;
            c2 = (3*(fb-fa)-2*da*w-db*w)/w2;
            c3 = (2*(fa-fb)+da*w+db*w)/w3;
            continue;
        }

        /*
         * Calculate spline and its derivatives using power basis
         */
        t = t-a;
        if( needy )
        {
            y->ptr.p_double[pointindex] = c0+t*(c1+t*(c2+t*c3));
        }
        if( needd1 )
        {
            d1->ptr.p_double[pointindex] = c1+2*t*c2+3*t*t*c3;
        }
        if( needd2 )
        {
            d2->ptr.p_double[pointindex] = 2*c2+6*t*c3;
        }
        pointindex = pointindex+1;
    }
}
}
