/* ---------------------------------------------------------------------------

  carray/carray_mathfunc_gsl.c

  This file is part of Ruby/CArray extension library.
  You can redistribute it and/or modify it under the terms of
  the GNU General Public License (GPL).

  Copyright (C) 2005-2008  Hiroki Motoyoshi

---------------------------------------------------------------------------- */

#include "ruby.h"
#include "carray.h"
#include <math.h>

#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

static gsl_rng *camath_gsl_rng;

/* ----------------------------------------------------------------------- */

#define rb_camath_gsl_func_d_d(grp, name)              \
  static void                                      \
  mathfunc_##name (void *p0, void *p1)             \
  {                                                \
    *(double *)p0 = gsl_##grp##_##name(*(double*)p1);   \
  }                                                \
  static VALUE                                     \
  rb_camath_##name (VALUE mod, VALUE rx)           \
  {                                                \
    return ca_call_cfunc_1_1(CA_DOUBLE, CA_DOUBLE, mathfunc_##name, rx); \
  }                                                \
  static VALUE                                     \
  rb_ca_##name (VALUE rx)           \
  {                                                \
    return ca_call_cfunc_1_1(CA_DOUBLE, CA_DOUBLE, mathfunc_##name, rx); \
  }                                                \

#define rb_camath_gsl_func_d_i(grp, name)              \
  static void                                      \
  mathfunc_##name (void *p0, void *p1)             \
  {                                                \
    *(double *)p0 = gsl_##grp##_##name(*(int*)p1);   \
  }                                                \
  static VALUE                                     \
  rb_camath_##name (VALUE mod, VALUE rx)           \
  {                                                \
    return ca_call_cfunc_1_1(CA_DOUBLE, CA_INT, mathfunc_##name, rx); \
  }                                                \
  static VALUE                                     \
  rb_ca_##name (VALUE rx)           \
  {                                                \
    return ca_call_cfunc_1_1(CA_DOUBLE, CA_INT, mathfunc_##name, rx); \
  }                                                \

#define rb_camath_gsl_func_d_dd(grp, name)           \
  static void                                      \
  mathfunc_##name (void *p0, void *p1, void *p2)            \
  {                                                \
    *(double *)p0 = gsl_##grp##_##name(*(double*)p1, *(double*)p2);   \
  }                                                \
  static VALUE                                     \
  rb_camath_##name (VALUE mod, VALUE rx1, VALUE rx2)      \
  {                                                \
    return ca_call_cfunc_1_2(CA_DOUBLE, CA_DOUBLE, CA_DOUBLE, mathfunc_##name, rx1, rx2); \
  }                                                \
  static VALUE                                     \
  rb_ca_##name (VALUE rx1, VALUE rx2)      \
  {                                                \
    return ca_call_cfunc_1_2(CA_DOUBLE, CA_DOUBLE, CA_DOUBLE, mathfunc_##name, rx1, rx2); \
  }                                                \

#define rb_camath_gsl_func_d_id(grp, name)           \
  static void                                      \
  mathfunc_##name (void *p0, void *p1, void *p2)            \
  {                                                \
    *(double *)p0 = gsl_##grp##_##name(*(int*)p1, *(double*)p2);   \
  }                                                \
  static VALUE                                     \
  rb_camath_##name (VALUE mod, VALUE rx1, VALUE rx2)      \
  {                                                \
    return ca_call_cfunc_1_2(CA_DOUBLE, CA_INT, CA_DOUBLE, mathfunc_##name, rx1, rx2); \
  }                                                \
  static VALUE                                     \
  rb_ca_##name (VALUE rx1, VALUE rx2)      \
  {                                                \
    return ca_call_cfunc_1_2(CA_DOUBLE, CA_INT, CA_DOUBLE, mathfunc_##name, rx1, rx2); \
  }                                                \

#define rb_camath_gsl_func_d_di(grp, name)           \
  static void                                      \
  mathfunc_##name (void *p0, void *p1, void *p2)            \
  {                                                \
    *(double *)p0 = gsl_##grp##_##name(*(double*)p1, *(int*)p2);   \
  }                                                \
  static VALUE                                     \
  rb_camath_##name (VALUE mod, VALUE rx1, VALUE rx2)      \
  {                                                \
    return ca_call_cfunc_1_2(CA_DOUBLE, CA_DOUBLE, CA_INT, mathfunc_##name, rx1, rx2); \
  }                                                \
  static VALUE                                     \
  rb_ca_##name (VALUE rx1, VALUE rx2)      \
  {                                                \
    return ca_call_cfunc_1_2(CA_DOUBLE, CA_DOUBLE, CA_INT, mathfunc_##name, rx1, rx2); \
  }                                                \

#define rb_camath_gsl_func_d_ii(grp, name)           \
  static void                                      \
  mathfunc_##name (void *p0, void *p1, void *p2)            \
  {                                                \
    *(double *)p0 = gsl_##grp##_##name(*(int*)p1, *(int*)p2);   \
  }                                                \
  static VALUE                                     \
  rb_camath_##name (VALUE mod, VALUE rx1, VALUE rx2)      \
  {                                                \
    return ca_call_cfunc_1_2(CA_DOUBLE, CA_INT, CA_INT, mathfunc_##name, rx1, rx2); \
  }                                                \
  static VALUE                                     \
  rb_ca_##name (VALUE rx1, VALUE rx2)      \
  {                                                \
    return ca_call_cfunc_1_2(CA_DOUBLE, CA_INT, CA_INT, mathfunc_##name, rx1, rx2); \
  }                                                \

#define rb_camath_gsl_func_d_ddd(grp, name)          \
  static void                                      \
  mathfunc_##name (void *p0, void *p1, void *p2, void *p3)  \
  {                                                \
    *(double *)p0 = gsl_##grp##_##name(*(double*)p1, *(double*)p2, *(double*)p3);   \
  }                                                \
  static VALUE                                     \
  rb_camath_##name (VALUE mod, VALUE rx1, VALUE rx2, VALUE rx3)   \
  {                                                \
    return ca_call_cfunc_1_3(CA_DOUBLE, CA_DOUBLE, CA_DOUBLE, CA_DOUBLE, mathfunc_##name, rx1, rx2, rx3); \
  }                                                \
  static VALUE                                     \
  rb_ca_##name (VALUE rx1, VALUE rx2, VALUE rx3)   \
  {                                                \
    return ca_call_cfunc_1_3(CA_DOUBLE, CA_DOUBLE, CA_DOUBLE, CA_DOUBLE, mathfunc_##name, rx1, rx2, rx3); \
  }                                                \

#define rb_camath_gsl_func_d_ddi(grp, name)          \
  static void                                      \
  mathfunc_##name (void *p0, void *p1, void *p2, void *p3)  \
  {                                                \
    *(double *)p0 = gsl_##grp##_##name(*(double*)p1, *(double*)p2, *(int*)p3);   \
  }                                                \
  static VALUE                                     \
  rb_camath_##name (VALUE mod, VALUE rx1, VALUE rx2, VALUE rx3)   \
  {                                                \
    return ca_call_cfunc_1_3(CA_DOUBLE, CA_DOUBLE, CA_DOUBLE, CA_INT, mathfunc_##name, rx1, rx2, rx3); \
  }                                                \


#define rb_camath_gsl_func_d_idd(grp, name)          \
  static void                                      \
  mathfunc_##name (void *p0, void *p1, void *p2, void *p3)  \
  {                                                \
    *(double *)p0 = gsl_##grp##_##name(*(int*)p1, *(double*)p2, *(double*)p3);   \
  }                                                \
  static VALUE                                     \
  rb_camath_##name (VALUE mod, VALUE rx1, VALUE rx2, VALUE rx3)   \
  {                                                \
    return ca_call_cfunc_1_3(CA_DOUBLE, CA_INT, CA_DOUBLE, CA_DOUBLE, mathfunc_##name, rx1, rx2, rx3); \
  }                                                \
  static VALUE                                     \
  rb_ca_##name (VALUE rx1, VALUE rx2, VALUE rx3)   \
  {                                                \
    return ca_call_cfunc_1_3(CA_DOUBLE, CA_INT, CA_DOUBLE, CA_DOUBLE, mathfunc_##name, rx1, rx2, rx3); \
  }                                                \

#define rb_camath_gsl_func_d_idi(grp, name)          \
  static void                                      \
  mathfunc_##name (void *p0, void *p1, void *p2, void *p3)  \
  {                                                \
    *(double *)p0 = gsl_##grp##_##name(*(int*)p1, *(double*)p2, *(int*)p3);   \
  }                                                \
  static VALUE                                     \
  rb_camath_##name (VALUE mod, VALUE rx1, VALUE rx2, VALUE rx3)   \
  {                                                \
    return ca_call_cfunc_1_3(CA_DOUBLE, CA_INT, CA_DOUBLE, CA_INT, mathfunc_##name, rx1, rx2, rx3); \
  }                                                \

#define rb_camath_gsl_func_d_iid(grp, name)          \
  static void                                      \
  mathfunc_##name (void *p0, void *p1, void *p2, void *p3)  \
  {                                                \
    *(double *)p0 = gsl_##grp##_##name(*(int*)p1, *(int*)p2, *(double*)p3);   \
  }                                                \
  static VALUE                                     \
  rb_camath_##name (VALUE mod, VALUE rx1, VALUE rx2, VALUE rx3)   \
  {                                                \
    return ca_call_cfunc_1_3(CA_DOUBLE, CA_INT, CA_INT, CA_DOUBLE, mathfunc_##name, rx1, rx2, rx3); \
  }                                                \
  static VALUE                                     \
  rb_ca_##name (VALUE rx1, VALUE rx2, VALUE rx3)   \
  {                                                \
    return ca_call_cfunc_1_3(CA_DOUBLE, CA_INT, CA_INT, CA_DOUBLE, mathfunc_##name, rx1, rx2, rx3); \
  }                                                \

#define rb_camath_gsl_func_d_dddd(grp, name)          \
  static void                                      \
  mathfunc_##name (void *p0, void *p1, void *p2, void *p3, void *p4)  \
  {                                                \
    *(double *)p0 = gsl_##grp##_##name(*(double*)p1, *(double*)p2, *(double*)p3, *(double*)p4);   \
  }                                                \
  static VALUE                                     \
  rb_camath_##name (VALUE mod, VALUE rx1, VALUE rx2, VALUE rx3, VALUE rx4)   \
  {                                                \
    return ca_call_cfunc_1_4(CA_DOUBLE, CA_DOUBLE, CA_DOUBLE, CA_DOUBLE, CA_DOUBLE, mathfunc_##name, rx1, rx2, rx3, rx4); \
  }                                                \
  static VALUE                                     \
  rb_ca_##name (VALUE rx1, VALUE rx2, VALUE rx3, VALUE rx4)   \
  {                                                \
    return ca_call_cfunc_1_4(CA_DOUBLE, CA_DOUBLE, CA_DOUBLE, CA_DOUBLE, CA_DOUBLE, mathfunc_##name, rx1, rx2, rx3, rx4); \
  }                                                \

#define rb_camath_gsl_func_d_dddi(grp, name)          \
  static void                                      \
  mathfunc_##name (void *p0, void *p1, void *p2, void *p3, void *p4)  \
  {                                                \
    *(double *)p0 = gsl_##grp##_##name(*(double*)p1, *(double*)p2, *(double*)p3, *(int*)p4);   \
  }                                                \
  static VALUE                                     \
  rb_camath_##name (VALUE mod, VALUE rx1, VALUE rx2, VALUE rx3, VALUE rx4)   \
  {                                                \
    return ca_call_cfunc_1_4(CA_DOUBLE, CA_DOUBLE, CA_DOUBLE, CA_DOUBLE, CA_INT, mathfunc_##name, rx1, rx2, rx3, rx4); \
  }                                                \


#define rb_camath_gsl_func_d_iidd(grp, name)          \
  static void                                      \
  mathfunc_##name (void *p0, void *p1, void *p2, void *p3, void *p4)  \
  {                                                \
    *(double *)p0 = gsl_##grp##_##name(*(int*)p1, *(int*)p2, *(double*)p3, *(double*)p4);   \
  }                                                \
  static VALUE                                     \
  rb_camath_##name (VALUE mod, VALUE rx1, VALUE rx2, VALUE rx3, VALUE rx4)   \
  {                                                \
    return ca_call_cfunc_1_4(CA_DOUBLE, CA_INT, CA_INT, CA_DOUBLE, CA_DOUBLE, mathfunc_##name, rx1, rx2, rx3, rx4); \
  }                                                \
  static VALUE                                     \
  rb_ca_##name (VALUE rx1, VALUE rx2, VALUE rx3, VALUE rx4)   \
  {                                                \
    return ca_call_cfunc_1_4(CA_DOUBLE, CA_INT, CA_INT, CA_DOUBLE, CA_DOUBLE, mathfunc_##name, rx1, rx2, rx3, rx4); \
  }                                                \

#define rb_camath_gsl_func_d_iiii(grp, name)          \
  static void                                      \
  mathfunc_##name (void *p0, void *p1, void *p2, void *p3, void *p4)  \
  {                                                \
    *(double *)p0 = gsl_##grp##_##name(*(int*)p1, *(int*)p2, *(int*)p3, *(int*)p4);   \
  }                                                \
  static VALUE                                     \
  rb_camath_##name (VALUE mod, VALUE rx1, VALUE rx2, VALUE rx3, VALUE rx4)   \
  {                                                \
    return ca_call_cfunc_1_4(CA_DOUBLE, CA_INT, CA_INT, CA_INT, CA_INT, mathfunc_##name, rx1, rx2, rx3, rx4); \
  }                                                \

#define rb_camath_gsl_func_d_ddddi(grp, name)          \
  static void                                      \
  mathfunc_##name (void *p0, void *p1, void *p2, void *p3, void *p4, void *p5)  \
  {                                                \
    *(double *)p0 = gsl_##grp##_##name(*(double*)p1, *(double*)p2, *(double*)p3, *(double*)p4, *(int*)p5);   \
  }                                                \
  static VALUE                                     \
  rb_camath_##name (VALUE mod, VALUE rx1, VALUE rx2, VALUE rx3, VALUE rx4, VALUE rx5)   \
  {                                                \
    return ca_call_cfunc_1_5(CA_DOUBLE, CA_DOUBLE, CA_DOUBLE, CA_DOUBLE, CA_DOUBLE, CA_INT, mathfunc_##name, rx1, rx2, rx3, rx4, rx5); \
  }                                                \


#define rb_camath_gsl_func_d_ddddd(grp, name)          \
  static void                                      \
  mathfunc_##name (void *p0, void *p1, void *p2, void *p3, void *p4, void *p5)  \
  {                                                \
    *(double *)p0 = gsl_##grp##_##name(*(double*)p1, *(double*)p2, *(double*)p3, *(double*)p4, *(double*)p5);   \
  }                                                \
  static VALUE                                     \
  rb_camath_##name (VALUE mod, VALUE rx1, VALUE rx2, VALUE rx3, VALUE rx4, VALUE rx5)   \
  {                                                \
    return ca_call_cfunc_1_5(CA_DOUBLE, CA_DOUBLE, CA_DOUBLE, CA_DOUBLE, CA_DOUBLE, CA_DOUBLE, mathfunc_##name, rx1, rx2, rx3, rx4, rx5); \
  }                                                \

#define rb_camath_gsl_func_d_iiiiii(grp, name)          \
  static void                                      \
  mathfunc_##name (void *p0, void *p1, void *p2, void *p3, void *p4, void *p5, void *p6)  \
  {                                                \
    *(double *)p0 = gsl_##grp##_##name(*(int*)p1, *(int*)p2, *(int*)p3, *(int*)p4, *(int*)p5, *(int*)p6);   \
  }                                                \
  static VALUE                                     \
  rb_camath_##name (VALUE mod, VALUE rx1, VALUE rx2, VALUE rx3, VALUE rx4, VALUE rx5, VALUE rx6)   \
  {                                                \
    return ca_call_cfunc_1_6(CA_DOUBLE, CA_INT, CA_INT, CA_INT, CA_INT, CA_INT, CA_INT, mathfunc_##name, rx1, rx2, rx3, rx4, rx5, rx6); \
  }                                                \

/* ----------------------------------------------------------------------- */

#define rb_camath_gsl_func_dd_dd(grp, name, ext)           \
  static void                                      \
  mathfunc_##name (void *p0, void *p1, void *p2, void *p3)            \
  {                                          \
    gsl_sf_result r0, r1;                                              \
    gsl_##grp##_##name##ext(*(double*)p2, *(double*)p3, &r0, &r1);   \
    *(double*) p0 = r0.val; *(double*) p1 = r1.val; \
  }                                                \
  static VALUE                                     \
  rb_camath_##name (VALUE mod, VALUE rx1, VALUE rx2)      \
  {                                                \
    return ca_call_cfunc_2_2(CA_DOUBLE, CA_DOUBLE, CA_DOUBLE, CA_DOUBLE, mathfunc_##name, rx1, rx2); \
  }                                                \

#define rb_camath_gsl_func_ddd_dd(grp, name, ext)           \
  static void                                      \
  mathfunc_##name (void *p0, void *p1, void *p2, void *p3, void *p4)            \
  {                                          \
    gsl_sf_result r0, r1, r2;                                              \
    gsl_##grp##_##name##ext(*(double*)p2, *(double*)p3, (double*)&r0, (double*)&r1, (double*)&r2);   \
    *(double*) p0 = r0.val; *(double*) p1 = r1.val; *(double*) p2 = r2.val; \
  }                                                \
  static VALUE                                     \
  rb_camath_##name (VALUE mod, VALUE rx1, VALUE rx2)      \
  {                                                \
    return ca_call_cfunc_3_2(CA_DOUBLE, CA_DOUBLE, CA_DOUBLE, CA_DOUBLE, CA_DOUBLE, mathfunc_##name, rx1, rx2); \
  }                                                \

/* ----------------------------------------------------------------------- */

#define rb_camath_gsl_rng_d_d(name)              \
  static void                                      \
  mathfunc_random_##name (void *p0, void *p1)             \
  {                                                \
    *(double *)p0 = gsl_ran_##name(camath_gsl_rng);    \
  }                                                \
  static VALUE                                     \
  rb_camath_random_##name (VALUE mod, VALUE rx)           \
  {                                                \
    return ca_call_cfunc_1_1(CA_DOUBLE, CA_DOUBLE, mathfunc_random_##name, rx);     \
  }                                                \

#define rb_camath_gsl_rng_d_dd(name)              \
  static void                                      \
  mathfunc_random_##name (void *p0, void *p1, void *p2)            \
  {                                                \
    *(double *)p0 = gsl_ran_##name(camath_gsl_rng, *(double*)p2);    \
  }                                                \
  static VALUE                                     \
  rb_camath_random_##name (VALUE mod, VALUE rx1, VALUE rx2)    \
  {                                                \
    return ca_call_cfunc_1_2(CA_DOUBLE, CA_DOUBLE, CA_DOUBLE, mathfunc_random_##name, rx1, rx2); \
  }                                                \

#define rb_camath_gsl_rng_d_ddd(name)              \
  static void                                      \
  mathfunc_random_##name (void *p0, void *p1, void *p2, void *p3)  \
  {                                                \
    *(double *)p0 = gsl_ran_##name(camath_gsl_rng, *(double*)p2, *(double*)p3);  \
  }                                                \
  static VALUE                                     \
  rb_camath_random_##name (VALUE mod, VALUE rx1, VALUE rx2, VALUE rx3)     \
  {                                                \
    return ca_call_cfunc_1_3(CA_DOUBLE, CA_DOUBLE, CA_DOUBLE, CA_DOUBLE, mathfunc_random_##name, rx1, rx2, rx3); \
  }                                                \

#define rb_camath_gsl_rng_d_ddi(name)              \
  static void                                      \
  mathfunc_random_##name (void *p0, void *p1, void *p2, void *p3)  \
  {                                                \
    *(double *)p0 = gsl_ran_##name(camath_gsl_rng, *(double*)p2, *(int*)p3);  \
  }                                                \
  static VALUE                                     \
  rb_camath_random_##name (VALUE mod, VALUE rx1, VALUE rx2, VALUE rx3)     \
  {                                                \
    return ca_call_cfunc_1_3(CA_DOUBLE, CA_DOUBLE, CA_DOUBLE, CA_INT, mathfunc_random_##name, rx1, rx2, rx3); \
  }                                                \

#define rb_camath_gsl_rng_d_dddd(name)              \
  static void                                      \
  mathfunc_random_##name (void *p0, void *p1, void *p2, void *p3, void *p4)  \
  {                                                \
    *(double *)p0 = gsl_ran_##name(camath_gsl_rng, *(double*)p2, *(double*)p3, *(double *)p4);  \
  }                                                \
  static VALUE                                     \
  rb_camath_random_##name (VALUE mod, VALUE rx1, VALUE rx2, VALUE rx3, VALUE rx4)     \
  {                                                \
    return ca_call_cfunc_1_4(CA_DOUBLE, CA_DOUBLE, CA_DOUBLE, CA_DOUBLE, CA_DOUBLE, mathfunc_random_##name, rx1, rx2, rx3, rx4); \
  }                                                \

#define rb_camath_gsl_rng_d_diii(name)              \
  static void                                      \
  mathfunc_random_##name (void *p0, void *p1, void *p2, void *p3, void *p4)  \
  {                                                \
    *(double *)p0 = gsl_ran_##name(camath_gsl_rng, *(int*)p2, *(int*)p3, *(int *)p4);  \
  }                                                \
  static VALUE                                     \
  rb_camath_random_##name (VALUE mod, VALUE rx1, VALUE rx2, VALUE rx3, VALUE rx4)     \
  {                                                \
    return ca_call_cfunc_1_4(CA_DOUBLE, CA_DOUBLE, CA_INT, CA_INT, CA_INT, mathfunc_random_##name, rx1, rx2, rx3, rx4); \
  }                                                \

#define rb_camath_gsl_rng_dd_ddd(name)              \
  static void                                      \
  mathfunc_random_##name (void *p0, void *p1, void *p2, void *p3, void *p4)  \
  {                                                \
    gsl_ran_##name(camath_gsl_rng, *(double*)p2, *(double*)p3, *(double *)p4, (double*)p0, (double*)p1);  \
  }                                                \
  static VALUE                                     \
  rb_camath_random_##name (VALUE mod, VALUE rx1, VALUE rx2, VALUE rx3)     \
  {                                                \
    return ca_call_cfunc_2_3(CA_DOUBLE, CA_DOUBLE, CA_DOUBLE, CA_DOUBLE, CA_DOUBLE, mathfunc_random_##name, rx1, rx2, rx3); \
  }                                                \

/* ----------------------------------------------------------------------- */

#define rb_define_camath_gsl_func(name, n)                             \
  rb_define_module_function(rb_mCAMath, #name, rb_camath_##name, n)

#define rb_define_camath_gsl_func2(name, n)                             \
  rb_define_module_function(rb_mCAMath, #name, rb_camath_##name, n); \
  rb_define_method(rb_cCArray, #name, rb_ca_##name, n-1)

#define rb_define_camath_gsl_func2_(name, n)                             \
  rb_define_module_function(rb_mCAMath, #name "_", rb_camath_##name, n); \
  rb_define_method(rb_cCArray, #name "_", rb_ca_##name, n-1)

#define rb_define_camath_gsl_func2u(name, n)                             \
  rb_define_module_function(rb_mCAMath, #name, rb_camath_##name, n); \
  rb_define_method(rb_cCArray, #name, rb_ca_##name, n-1); \
  rb_undef(rb_cCArray, rb_intern(#name)); \

/* ----------------------------------------------------------------------- */

static VALUE
rb_camath_random_dir_2d (VALUE mod, VALUE vn)
{
  CArray *cx, *cy;
  double *px, *py;
  int32_t elements = NUM2INT(vn), i;
  cx = carray_new(CA_DOUBLE, 1, &elements, 0, NULL);
  cy = carray_new(CA_DOUBLE, 1, &elements, 0, NULL);
  px = (double*) cx->ptr;
  py = (double*) cy->ptr;
  for (i=0; i<elements; i++) {
    gsl_ran_dir_2d(camath_gsl_rng, px, py) ;
    px++; py++;
  }
  if ( elements == 1 ) {
    volatile VALUE out = rb_ary_new3(2, rb_float_new(*(double*)cx->ptr),
                          rb_float_new(*(double*)cy->ptr));
    ca_free(cx);
    ca_free(cy);
    return out;
  }
  else {
    return rb_ary_new3(2, ca_wrap_struct(cx), ca_wrap_struct(cy));
  }
}

static VALUE
rb_camath_random_dir_3d (VALUE mod, VALUE vn)
{
  CArray *cx, *cy, *cz;
  double *px, *py, *pz;
  int32_t elements = NUM2INT(vn), i;
  cx = carray_new(CA_DOUBLE, 1, &elements, 0, NULL);
  cy = carray_new(CA_DOUBLE, 1, &elements, 0, NULL);
  cz = carray_new(CA_DOUBLE, 1, &elements, 0, NULL);
  px = (double*) cx->ptr;
  py = (double*) cy->ptr;
  pz = (double*) cz->ptr;
  for (i=0; i<elements; i++) {
    gsl_ran_dir_3d(camath_gsl_rng, px, py, pz) ;
    px++; py++; pz++;
  }
  if ( elements == 1 ) {
    volatile VALUE out = rb_ary_new3(3, rb_float_new(*(double*)cx->ptr),
                                  rb_float_new(*(double*)cy->ptr),
                                  rb_float_new(*(double*)cz->ptr));
    ca_free(cx);
    ca_free(cy);
    ca_free(cz);
    return out;
  }
  else {
    return rb_ary_new3(3, ca_wrap_struct(cx), 
                          ca_wrap_struct(cy), 
                          ca_wrap_struct(cz));
  }
}

static VALUE
rb_camath_random_dir_nd (VALUE mod, VALUE vn, VALUE vd)
{
  volatile VALUE out;
  CArray **ca;
  double **p;
  double *pv;
  int32_t d = NUM2INT(vd), elements = NUM2INT(vn), k, i;
  if ( d < 0 ) {
    rb_raise(rb_eArgError, "dimension should be positive number");
  }
  if ( elements < 0 ) {
    rb_raise(rb_eArgError, "elements should be positive number");
  }
  ca = ALLOC_N(CArray*, d);
  p  = ALLOC_N(double*, d);
  pv = ALLOC_N(double, d);
  for (k=0; k<d; k++) {
    ca[k] = carray_new(CA_DOUBLE, 1, &elements, 0, NULL);
    p[k] = (double*) ca[k]->ptr;
  }
  for (i=0; i<elements; i++) {
    gsl_ran_dir_nd(camath_gsl_rng, d, pv);
    for (k=0; k<d; k++) {
      *(p[k]) = pv[k];
      p[k]++;
    }
  }
  out = rb_ary_new();
  if ( elements == 1 ) {
    for (k=0; k<d; k++) {
      rb_ary_store(out, k, rb_float_new(*(double*)ca[k]->ptr));
      ca_free(ca[k]);
    }
   }
  else {
    for (k=0; k<d; k++) {
      rb_ary_store(out, k, ca_wrap_struct(ca[k]));
    }
  }
  free(ca);
  free(p);
  free(pv);
  return out;
}

static VALUE
rb_camath_random_multinomial (VALUE mod, 
                              VALUE vn, VALUE vm, volatile VALUE vpk)
{
  volatile VALUE out;
  CArray *cpk;
  CArray **ca;
  unsigned int **pn;
  unsigned int *pv;
  int32_t elements, i, k, d, m;
  m = NUM2INT(vm);
  if ( m <0 ) {
    rb_raise(rb_eArgError, "# of tries should be positive number");
  }
  elements = NUM2INT(vn);
  if ( elements < 0 ) {
    rb_raise(rb_eArgError, "elements should be positive number");
  }

  cpk = ca_wrap_readonly(vpk, CA_DOUBLE);
  d = cpk->elements;

  ca = ALLOC_N(CArray*, d);
  pn = ALLOC_N(unsigned int*, d);
  pv = ALLOC_N(unsigned int, d);

  for (k=0; k<d; k++) {
    ca[k] = carray_new(CA_INT32, 1, &elements, 0, NULL);
    pn[k] = (unsigned int*) ca[k]->ptr;
  }
  
  ca_attach(cpk);
  for (i=0; i<elements; i++) {
    gsl_ran_multinomial(camath_gsl_rng, d, m, (double*)cpk->ptr, pv);
    for (k=0; k<d; k++) {
      *(pn[k]) = pv[k];
      pn[k]++;
    }
  }
  ca_detach(cpk);
  
  out = rb_ary_new();
  if ( elements == 1 ) {
    for (k=0; k<d; k++) {
      rb_ary_store(out, k, INT2NUM(*(unsigned int*)ca[k]->ptr));
      ca_free(ca[k]);
    }
   }
  else {
    for (k=0; k<d; k++) {
      rb_ary_store(out, k, ca_wrap_struct(ca[k]));
    }
  }
  free(ca);
  free(pn);
  free(pv);
  return out;
}

static VALUE
rb_camath_multinomial_pdf (VALUE mod, 
                           VALUE vnk, volatile VALUE vpk)
{
  volatile VALUE out, tmp;
  CArray *cpk, *cy;
  CArray **ca;
  unsigned int **pn;
  unsigned int *pv;
  double *py;
  int32_t elements, i, k, d;

  Check_Type(vnk, T_ARRAY);

  cpk = ca_wrap_readonly(vpk, CA_DOUBLE);
  d = cpk->elements;
  
  if ( RARRAY_LEN(vnk) != d ) {
    rb_raise(rb_eArgError, "# of nk should equal # of pk");
  }

  ca = ALLOC_N(CArray*, d);
  pn = ALLOC_N(unsigned int*, d);
  pv = ALLOC_N(unsigned int, d);

  for (k=0; k<d; k++) {
    tmp = rb_ary_entry(vnk, k);
    ca[k] = ca_wrap_readonly(tmp, CA_INT32);
  }

  elements = ca[0]->elements;
  for (k=1; k<d; k++) {
    if ( elements != ca[k]->elements ) {
      rb_raise(rb_eArgError, "element size mismatch in nk");
    }
  }
    
  for (k=0; k<d; k++) {
    ca_attach(ca[k]);
    pn[k] = (unsigned int*) ca[k]->ptr;
  }

  cy = carray_new(CA_DOUBLE, 1, &elements, 0, NULL);
  py = (double*)cy->ptr;

  ca_attach(cpk);
  for (i=0; i<elements; i++) {
    for (k=0; k<d; k++) {
      pv[k] = *(pn[k]);
      pn[k]++;
    }
    *py = gsl_ran_multinomial_pdf(d, (double*)cpk->ptr, pv);
    py++;
  }
  ca_detach(cpk);

  if ( elements == 1 ) {
    out = rb_float_new(*(double*)cy->ptr);
    ca_free(cy);
  }
  else {
    out = ca_wrap_struct(cy);
  }

  for (k=0; k<d; k++) {
    ca_detach(ca[k]);
  }

  free(ca);
  free(pn);
  free(pv);
  return out;
}

static VALUE
rb_camath_multinomial_lnpdf (VALUE mod, 
                           VALUE vnk, volatile VALUE vpk)
{
  volatile VALUE out, tmp;
  CArray *cpk, *cy;
  CArray **ca;
  unsigned int **pn;
  unsigned int *pv;
  double *py;
  int32_t elements, i, k, d;

  Check_Type(vnk, T_ARRAY);

  cpk = ca_wrap_readonly(vpk, CA_DOUBLE);
  d = cpk->elements;
  
  if ( RARRAY_LEN(vnk) != d ) {
    rb_raise(rb_eArgError, "# of nk should equal # of pk");
  }

  ca = ALLOC_N(CArray*, d);
  pn = ALLOC_N(unsigned int*, d);
  pv = ALLOC_N(unsigned int, d);

  for (k=0; k<d; k++) {
    tmp = rb_ary_entry(vnk, k);
    ca[k] = ca_wrap_readonly(tmp, CA_INT32);
  }

  elements = ca[0]->elements;
  for (k=1; k<d; k++) {
    if ( elements != ca[k]->elements ) {
      rb_raise(rb_eArgError, "element size mismatch in nk");
    }
  }
    
  for (k=0; k<d; k++) {
    ca_attach(ca[k]);
    pn[k] = (unsigned int*) ca[k]->ptr;
  }

  cy = carray_new(CA_DOUBLE, 1, &elements, 0, NULL);
  py = (double*)cy->ptr;

  ca_attach(cpk);
  for (i=0; i<elements; i++) {
    for (k=0; k<d; k++) {
      pv[k] = *(pn[k]);
      pn[k]++;
    }
    *py = gsl_ran_multinomial_lnpdf(d, (double*)cpk->ptr, pv);
    py++;
  }
  ca_detach(cpk);

  if ( elements == 1 ) {
    out = rb_float_new(*(double*)cy->ptr);
    ca_free(cy);
  }
  else {
    out = ca_wrap_struct(cy);
  }

  for (k=0; k<d; k++) {
    ca_detach(ca[k]);
  }

  free(ca);
  free(pn);
  free(pv);
  return out;
}

/* ----------------------------------------------------------------------- */

rb_camath_gsl_func_d_di(sf, airy_Ai);
rb_camath_gsl_func_d_di(sf, airy_Bi);
rb_camath_gsl_func_d_di(sf, airy_Ai_scaled);
rb_camath_gsl_func_d_di(sf, airy_Bi_scaled);
rb_camath_gsl_func_d_di(sf, airy_Ai_deriv);
rb_camath_gsl_func_d_di(sf, airy_Bi_deriv);
rb_camath_gsl_func_d_di(sf, airy_Ai_deriv_scaled);
rb_camath_gsl_func_d_di(sf, airy_Bi_deriv_scaled);
rb_camath_gsl_func_d_i(sf, airy_zero_Ai);
rb_camath_gsl_func_d_i(sf, airy_zero_Bi);

rb_camath_gsl_func_d_d(sf, bessel_J0);
rb_camath_gsl_func_d_d(sf, bessel_J1);
rb_camath_gsl_func_d_id(sf, bessel_Jn);

rb_camath_gsl_func_d_d(sf, bessel_Y0);
rb_camath_gsl_func_d_d(sf, bessel_Y1);
rb_camath_gsl_func_d_id(sf, bessel_Yn);

rb_camath_gsl_func_d_d(sf, bessel_I0);
rb_camath_gsl_func_d_d(sf, bessel_I1);
rb_camath_gsl_func_d_id(sf, bessel_In);
rb_camath_gsl_func_d_d(sf, bessel_I0_scaled);
rb_camath_gsl_func_d_d(sf, bessel_I1_scaled);
rb_camath_gsl_func_d_id(sf, bessel_In_scaled);

rb_camath_gsl_func_d_d(sf, bessel_K0);
rb_camath_gsl_func_d_d(sf, bessel_K1);
rb_camath_gsl_func_d_id(sf, bessel_Kn);
rb_camath_gsl_func_d_d(sf, bessel_K0_scaled);
rb_camath_gsl_func_d_d(sf, bessel_K1_scaled);
rb_camath_gsl_func_d_id(sf, bessel_Kn_scaled);

rb_camath_gsl_func_d_d(sf, bessel_j0);
rb_camath_gsl_func_d_d(sf, bessel_j1);
rb_camath_gsl_func_d_d(sf, bessel_j2);
rb_camath_gsl_func_d_id(sf, bessel_jl);

rb_camath_gsl_func_d_d(sf, bessel_y0);
rb_camath_gsl_func_d_d(sf, bessel_y1);
rb_camath_gsl_func_d_d(sf, bessel_y2);
rb_camath_gsl_func_d_id(sf, bessel_yl);

rb_camath_gsl_func_d_d(sf, bessel_i0_scaled);
rb_camath_gsl_func_d_d(sf, bessel_i1_scaled);
rb_camath_gsl_func_d_d(sf, bessel_i2_scaled);
rb_camath_gsl_func_d_id(sf, bessel_il_scaled);

rb_camath_gsl_func_d_d(sf, bessel_k0_scaled);
rb_camath_gsl_func_d_d(sf, bessel_k1_scaled);
rb_camath_gsl_func_d_d(sf, bessel_k2_scaled);
rb_camath_gsl_func_d_id(sf, bessel_kl_scaled);

rb_camath_gsl_func_d_dd(sf, bessel_Jnu);
rb_camath_gsl_func_d_dd(sf, bessel_Ynu);
rb_camath_gsl_func_d_dd(sf, bessel_Inu);
rb_camath_gsl_func_d_dd(sf, bessel_Inu_scaled);
rb_camath_gsl_func_d_dd(sf, bessel_Knu);
rb_camath_gsl_func_d_dd(sf, bessel_Knu_scaled);

rb_camath_gsl_func_d_i(sf, bessel_zero_J0);
rb_camath_gsl_func_d_i(sf, bessel_zero_J1);
rb_camath_gsl_func_d_di(sf, bessel_zero_Jnu);

rb_camath_gsl_func_d_d(sf, clausen);

rb_camath_gsl_func_d_dd(sf, hydrogenicR_1);
rb_camath_gsl_func_d_iidd(sf, hydrogenicR);
/* not implemented :  Coulomb Wave Functions */
/* not implemented :  Coulomb Wave Function Normalization Constant */

rb_camath_gsl_func_d_iiiiii(sf, coupling_3j);
rb_camath_gsl_func_d_iiiiii(sf, coupling_6j);

/* rb_camath_gsl_func_d_iiiiiiiii(sf, coupling_9j); */

rb_camath_gsl_func_d_d(sf, dawson);

rb_camath_gsl_func_d_d(sf, debye_1);
rb_camath_gsl_func_d_d(sf, debye_2);
rb_camath_gsl_func_d_d(sf, debye_3);
rb_camath_gsl_func_d_d(sf, debye_4);
rb_camath_gsl_func_d_d(sf, debye_5);
rb_camath_gsl_func_d_d(sf, debye_6);

rb_camath_gsl_func_d_d(sf, dilog);
/* rb_camath_gsl_func_d_d(sf, complex_dilog); */

rb_camath_gsl_func_d_di(sf, ellint_Kcomp);
rb_camath_gsl_func_d_di(sf, ellint_Ecomp);
rb_camath_gsl_func_d_ddi(sf, ellint_Pcomp);
rb_camath_gsl_func_d_ddi(sf, ellint_F);
rb_camath_gsl_func_d_ddi(sf, ellint_E);
rb_camath_gsl_func_d_dddi(sf, ellint_P);
rb_camath_gsl_func_d_dddi(sf, ellint_D);
rb_camath_gsl_func_d_ddi(sf, ellint_RC);
rb_camath_gsl_func_d_dddi(sf, ellint_RD);
rb_camath_gsl_func_d_dddi(sf, ellint_RF);
rb_camath_gsl_func_d_ddddi(sf, ellint_RJ);

rb_camath_gsl_func_ddd_dd(sf, elljac, _e);

rb_camath_gsl_func_d_d(sf, erf);
rb_camath_gsl_func_d_d(sf, erf_Z);
rb_camath_gsl_func_d_d(sf, erf_Q);
rb_camath_gsl_func_d_d(sf, erfc);
rb_camath_gsl_func_d_d(sf, log_erfc);
rb_camath_gsl_func_d_d(sf, hazard);

/* rb_camath_gsl_func_d_dd(sf, exp); */
/* rb_camath_gsl_func_d_dd(sf, exp_mult); */
/* rb_camath_gsl_func_d_dd(sf, expm1); */
/* rb_camath_gsl_func_d_dd(sf, exprel); */
/* rb_camath_gsl_func_d_dd(sf, exprel_2); */
/* rb_camath_gsl_func_d_dd(sf, exprel_n); */

rb_camath_gsl_func_d_d(sf, expint_E1);
rb_camath_gsl_func_d_d(sf, expint_E2);
rb_camath_gsl_func_d_id(sf, expint_En);
rb_camath_gsl_func_d_d(sf, expint_Ei);
rb_camath_gsl_func_d_d(sf, Shi);
rb_camath_gsl_func_d_d(sf, Chi);
rb_camath_gsl_func_d_d(sf, expint_3);
rb_camath_gsl_func_d_d(sf, Si);
rb_camath_gsl_func_d_d(sf, Ci);
rb_camath_gsl_func_d_d(sf, atanint);

rb_camath_gsl_func_d_d(sf, fermi_dirac_m1);
rb_camath_gsl_func_d_d(sf, fermi_dirac_0);
rb_camath_gsl_func_d_d(sf, fermi_dirac_1);
rb_camath_gsl_func_d_d(sf, fermi_dirac_2);
rb_camath_gsl_func_d_id(sf, fermi_dirac_int);
rb_camath_gsl_func_d_d(sf, fermi_dirac_mhalf);
rb_camath_gsl_func_d_d(sf, fermi_dirac_half);
rb_camath_gsl_func_d_d(sf, fermi_dirac_3half);
rb_camath_gsl_func_d_dd(sf, fermi_dirac_inc_0);

rb_camath_gsl_func_d_d(sf, gamma);
rb_camath_gsl_func_d_d(sf, lngamma);
rb_camath_gsl_func_d_d(sf, gammastar);
rb_camath_gsl_func_d_d(sf, gammainv);

rb_camath_gsl_func_d_dd(sf, gamma_inc);
rb_camath_gsl_func_d_dd(sf, gamma_inc_Q);
rb_camath_gsl_func_d_dd(sf, gamma_inc_P);

rb_camath_gsl_func_d_i(sf, fact);
rb_camath_gsl_func_d_i(sf, doublefact);
rb_camath_gsl_func_d_i(sf, lnfact);
rb_camath_gsl_func_d_i(sf, lndoublefact);
rb_camath_gsl_func_d_ii(sf, choose);
rb_camath_gsl_func_d_ii(sf, lnchoose);
rb_camath_gsl_func_d_id(sf, taylorcoeff);

rb_camath_gsl_func_d_dd(sf, beta);
rb_camath_gsl_func_d_dd(sf, lnbeta);
rb_camath_gsl_func_d_ddd(sf, beta_inc);

rb_camath_gsl_func_d_dd(sf, gegenpoly_1);
rb_camath_gsl_func_d_dd(sf, gegenpoly_2);
rb_camath_gsl_func_d_dd(sf, gegenpoly_3);
rb_camath_gsl_func_d_idd(sf, gegenpoly_n);

rb_camath_gsl_func_d_dd(sf, hyperg_0F1);
rb_camath_gsl_func_d_iid(sf, hyperg_1F1_int);
rb_camath_gsl_func_d_ddd(sf, hyperg_1F1);
rb_camath_gsl_func_d_iid(sf, hyperg_U_int);
rb_camath_gsl_func_d_ddd(sf, hyperg_U);
rb_camath_gsl_func_d_dddd(sf, hyperg_2F1);
rb_camath_gsl_func_d_dddd(sf, hyperg_2F1_conj);
rb_camath_gsl_func_d_dddd(sf, hyperg_2F1_renorm);
rb_camath_gsl_func_d_dddd(sf, hyperg_2F1_conj_renorm);
rb_camath_gsl_func_d_ddd(sf, hyperg_2F0);

rb_camath_gsl_func_d_dd(sf, laguerre_1);
rb_camath_gsl_func_d_dd(sf, laguerre_2);
rb_camath_gsl_func_d_dd(sf, laguerre_3);
rb_camath_gsl_func_d_idd(sf, laguerre_n);

rb_camath_gsl_func_d_d(sf, lambert_W0);
rb_camath_gsl_func_d_d(sf, lambert_Wm1);

rb_camath_gsl_func_d_d(sf, legendre_P1);
rb_camath_gsl_func_d_d(sf, legendre_P2);
rb_camath_gsl_func_d_d(sf, legendre_P3);
rb_camath_gsl_func_d_id(sf, legendre_Pl);
rb_camath_gsl_func_d_d(sf, legendre_Q0);
rb_camath_gsl_func_d_d(sf, legendre_Q1);
rb_camath_gsl_func_d_id(sf, legendre_Ql);
rb_camath_gsl_func_d_iid(sf, legendre_Plm);
rb_camath_gsl_func_d_iid(sf, legendre_sphPlm);

rb_camath_gsl_func_d_dd(sf, conicalP_half);
rb_camath_gsl_func_d_dd(sf, conicalP_mhalf);
rb_camath_gsl_func_d_dd(sf, conicalP_0);
rb_camath_gsl_func_d_dd(sf, conicalP_1);
rb_camath_gsl_func_d_idd(sf, conicalP_sph_reg);
rb_camath_gsl_func_d_idd(sf, conicalP_cyl_reg);

rb_camath_gsl_func_d_dd(sf, legendre_H3d_0);
rb_camath_gsl_func_d_dd(sf, legendre_H3d_1);
rb_camath_gsl_func_d_idd(sf, legendre_H3d);

/* rb_camath_gsl_func_d_d(sf, log); */
/* rb_camath_gsl_func_d_d(sf, log_abs); */
/* rb_camath_gsl_func_dd_dd(sf, complex_log_e); */
/* rb_camath_gsl_func_dd_dd(sf, complex_log_1plusx); */
/* rb_camath_gsl_func_dd_dd(sf, complex_log_1plusx_mx); */

/* not implelmeted : Mathieu Functions */

/* rb_camath_gsl_func_d_di(sf, pow_int); */

rb_camath_gsl_func_d_d(sf, psi);
rb_camath_gsl_func_d_d(sf, psi_1piy);
rb_camath_gsl_func_d_d(sf, psi_1);
rb_camath_gsl_func_d_id(sf, psi_n);

rb_camath_gsl_func_d_d(sf, synchrotron_1);
rb_camath_gsl_func_d_d(sf, synchrotron_2);

rb_camath_gsl_func_d_d(sf, transport_2);
rb_camath_gsl_func_d_d(sf, transport_3);
rb_camath_gsl_func_d_d(sf, transport_4);
rb_camath_gsl_func_d_d(sf, transport_5);

/* rb_camath_gsl_func_d_d(sf, sin); */
/* rb_camath_gsl_func_d_d(sf, cos); */
/* rb_camath_gsl_func_d_dd(sf, hypot); */
rb_camath_gsl_func_d_d(sf, sinc);
/* rb_camath_gsl_func_dd_dd(sf, complex_sin); */
/* rb_camath_gsl_func_dd_dd(sf, complex_cos); */
/* rb_camath_gsl_func_dd_dd(sf, complex_logsin); */
/* rb_camath_gsl_func_d_d(sf, lnsinh); */
/* rb_camath_gsl_func_d_d(sf, lncosh); */
rb_camath_gsl_func_dd_dd(sf, polar_to_rect,);
rb_camath_gsl_func_dd_dd(sf, rect_to_polar,);
rb_camath_gsl_func_d_d(sf, angle_restrict_symm);
rb_camath_gsl_func_d_d(sf, angle_restrict_pos);

rb_camath_gsl_func_d_d(sf, zeta);
rb_camath_gsl_func_d_d(sf, zetam1);
rb_camath_gsl_func_d_dd(sf, hzeta);
rb_camath_gsl_func_d_d(sf, eta);

/* random distributions */

rb_camath_gsl_rng_d_dd(gaussian);
rb_camath_gsl_func_d_dd(ran, gaussian_pdf);
rb_camath_gsl_func_d_dd(cdf, gaussian_P);
rb_camath_gsl_func_d_dd(cdf, gaussian_Q);
rb_camath_gsl_func_d_dd(cdf, gaussian_Pinv);
rb_camath_gsl_func_d_dd(cdf, gaussian_Qinv);

rb_camath_gsl_rng_d_d(ugaussian);
rb_camath_gsl_func_d_d(ran, ugaussian_pdf);
rb_camath_gsl_func_d_d(cdf, ugaussian_P);
rb_camath_gsl_func_d_d(cdf, ugaussian_Q);
rb_camath_gsl_func_d_d(cdf, ugaussian_Pinv);
rb_camath_gsl_func_d_d(cdf, ugaussian_Qinv);

rb_camath_gsl_rng_d_ddd(gaussian_tail);
rb_camath_gsl_func_d_ddd(ran, gaussian_tail_pdf);

rb_camath_gsl_rng_d_dd(ugaussian_tail);
rb_camath_gsl_func_d_dd(ran, ugaussian_tail_pdf);

rb_camath_gsl_rng_dd_ddd(bivariate_gaussian);
rb_camath_gsl_func_d_ddddd(ran, bivariate_gaussian_pdf);

rb_camath_gsl_rng_d_dd(exponential);
rb_camath_gsl_func_d_dd(ran, exponential_pdf);
rb_camath_gsl_func_d_dd(cdf, exponential_P);
rb_camath_gsl_func_d_dd(cdf, exponential_Q);
rb_camath_gsl_func_d_dd(cdf, exponential_Pinv);
rb_camath_gsl_func_d_dd(cdf, exponential_Qinv);

rb_camath_gsl_rng_d_dd(laplace);
rb_camath_gsl_func_d_dd(ran, laplace_pdf);
rb_camath_gsl_func_d_dd(cdf, laplace_P);
rb_camath_gsl_func_d_dd(cdf, laplace_Q);
rb_camath_gsl_func_d_dd(cdf, laplace_Pinv);
rb_camath_gsl_func_d_dd(cdf, laplace_Qinv);

rb_camath_gsl_rng_d_ddd(exppow);
rb_camath_gsl_func_d_ddd(ran, exppow_pdf);
rb_camath_gsl_func_d_ddd(cdf, exppow_P);
rb_camath_gsl_func_d_ddd(cdf, exppow_Q);

rb_camath_gsl_rng_d_dd(cauchy);
rb_camath_gsl_func_d_dd(ran, cauchy_pdf);
rb_camath_gsl_func_d_dd(cdf, cauchy_P);
rb_camath_gsl_func_d_dd(cdf, cauchy_Q);
rb_camath_gsl_func_d_dd(cdf, cauchy_Pinv);
rb_camath_gsl_func_d_dd(cdf, cauchy_Qinv);

rb_camath_gsl_rng_d_dd(rayleigh);
rb_camath_gsl_func_d_dd(ran, rayleigh_pdf);
rb_camath_gsl_func_d_dd(cdf, rayleigh_P);
rb_camath_gsl_func_d_dd(cdf, rayleigh_Q);
rb_camath_gsl_func_d_dd(cdf, rayleigh_Pinv);
rb_camath_gsl_func_d_dd(cdf, rayleigh_Qinv);

rb_camath_gsl_rng_d_ddd(rayleigh_tail);
rb_camath_gsl_func_d_ddd(ran, rayleigh_tail_pdf);

rb_camath_gsl_rng_d_d(landau);
rb_camath_gsl_func_d_d(ran, landau_pdf);

rb_camath_gsl_rng_d_ddd(levy);
rb_camath_gsl_rng_d_dddd(levy_skew);

rb_camath_gsl_rng_d_ddd(gamma);
rb_camath_gsl_func_d_ddd(ran, gamma_pdf);
rb_camath_gsl_func_d_ddd(cdf, gamma_P);
rb_camath_gsl_func_d_ddd(cdf, gamma_Q);
rb_camath_gsl_func_d_ddd(cdf, gamma_Pinv);
rb_camath_gsl_func_d_ddd(cdf, gamma_Qinv);

rb_camath_gsl_rng_d_ddd(erlang);
rb_camath_gsl_func_d_ddd(ran, erlang_pdf);

rb_camath_gsl_rng_d_ddd(flat);
rb_camath_gsl_func_d_ddd(ran, flat_pdf);
rb_camath_gsl_func_d_ddd(cdf, flat_P);
rb_camath_gsl_func_d_ddd(cdf, flat_Q);
rb_camath_gsl_func_d_ddd(cdf, flat_Pinv);
rb_camath_gsl_func_d_ddd(cdf, flat_Qinv);

rb_camath_gsl_rng_d_ddd(lognormal);
rb_camath_gsl_func_d_ddd(ran, lognormal_pdf);
rb_camath_gsl_func_d_ddd(cdf, lognormal_P);
rb_camath_gsl_func_d_ddd(cdf, lognormal_Q);
rb_camath_gsl_func_d_ddd(cdf, lognormal_Pinv);
rb_camath_gsl_func_d_ddd(cdf, lognormal_Qinv);

rb_camath_gsl_rng_d_dd(chisq);
rb_camath_gsl_func_d_dd(ran, chisq_pdf);
rb_camath_gsl_func_d_dd(cdf, chisq_P);
rb_camath_gsl_func_d_dd(cdf, chisq_Q);
rb_camath_gsl_func_d_dd(cdf, chisq_Pinv);
rb_camath_gsl_func_d_dd(cdf, chisq_Qinv);

rb_camath_gsl_rng_d_ddd(fdist);
rb_camath_gsl_func_d_ddd(ran, fdist_pdf);
rb_camath_gsl_func_d_ddd(cdf, fdist_P);
rb_camath_gsl_func_d_ddd(cdf, fdist_Q);
rb_camath_gsl_func_d_ddd(cdf, fdist_Pinv);
rb_camath_gsl_func_d_ddd(cdf, fdist_Qinv);

rb_camath_gsl_rng_d_dd(tdist);
rb_camath_gsl_func_d_dd(ran, tdist_pdf);
rb_camath_gsl_func_d_dd(cdf, tdist_P);
rb_camath_gsl_func_d_dd(cdf, tdist_Q);
rb_camath_gsl_func_d_dd(cdf, tdist_Pinv);
rb_camath_gsl_func_d_dd(cdf, tdist_Qinv);

rb_camath_gsl_rng_d_ddd(beta);
rb_camath_gsl_func_d_ddd(ran, beta_pdf);
rb_camath_gsl_func_d_ddd(cdf, beta_P);
rb_camath_gsl_func_d_ddd(cdf, beta_Q);
rb_camath_gsl_func_d_ddd(cdf, beta_Pinv);
rb_camath_gsl_func_d_ddd(cdf, beta_Qinv);

rb_camath_gsl_rng_d_dd(logistic);
rb_camath_gsl_func_d_dd(ran, logistic_pdf);
rb_camath_gsl_func_d_dd(cdf, logistic_P);
rb_camath_gsl_func_d_dd(cdf, logistic_Q);
rb_camath_gsl_func_d_dd(cdf, logistic_Pinv);
rb_camath_gsl_func_d_dd(cdf, logistic_Qinv);

rb_camath_gsl_rng_d_ddd(pareto);
rb_camath_gsl_func_d_ddd(ran, pareto_pdf);
rb_camath_gsl_func_d_ddd(cdf, pareto_P);
rb_camath_gsl_func_d_ddd(cdf, pareto_Q);
rb_camath_gsl_func_d_ddd(cdf, pareto_Pinv);
rb_camath_gsl_func_d_ddd(cdf, pareto_Qinv);

rb_camath_gsl_rng_d_ddd(weibull);
rb_camath_gsl_func_d_ddd(ran, weibull_pdf);
rb_camath_gsl_func_d_ddd(cdf, weibull_P);
rb_camath_gsl_func_d_ddd(cdf, weibull_Q);
rb_camath_gsl_func_d_ddd(cdf, weibull_Pinv);
rb_camath_gsl_func_d_ddd(cdf, weibull_Qinv);

rb_camath_gsl_rng_d_ddd(gumbel1);
rb_camath_gsl_func_d_ddd(ran, gumbel1_pdf);
rb_camath_gsl_func_d_ddd(cdf, gumbel1_P);
rb_camath_gsl_func_d_ddd(cdf, gumbel1_Q);
rb_camath_gsl_func_d_ddd(cdf, gumbel1_Pinv);
rb_camath_gsl_func_d_ddd(cdf, gumbel1_Qinv);

rb_camath_gsl_rng_d_ddd(gumbel2);
rb_camath_gsl_func_d_ddd(ran, gumbel2_pdf);
rb_camath_gsl_func_d_ddd(cdf, gumbel2_P);
rb_camath_gsl_func_d_ddd(cdf, gumbel2_Q);
rb_camath_gsl_func_d_ddd(cdf, gumbel2_Pinv);
rb_camath_gsl_func_d_ddd(cdf, gumbel2_Qinv);

rb_camath_gsl_rng_d_dd(poisson);
rb_camath_gsl_func_d_id(ran, poisson_pdf);
rb_camath_gsl_func_d_id(cdf, poisson_P);
rb_camath_gsl_func_d_id(cdf, poisson_Q);

rb_camath_gsl_rng_d_dd(bernoulli);
rb_camath_gsl_func_d_id(ran, bernoulli_pdf);

rb_camath_gsl_rng_d_ddi(binomial);
rb_camath_gsl_func_d_idi(ran, binomial_pdf);
rb_camath_gsl_func_d_idi(cdf, binomial_P);
rb_camath_gsl_func_d_idi(cdf, binomial_Q);

rb_camath_gsl_rng_d_ddd(negative_binomial);
rb_camath_gsl_func_d_idd(ran, negative_binomial_pdf);
rb_camath_gsl_func_d_idd(cdf, negative_binomial_P);
rb_camath_gsl_func_d_idd(cdf, negative_binomial_Q);

rb_camath_gsl_rng_d_ddi(pascal);
rb_camath_gsl_func_d_idi(ran, pascal_pdf);
rb_camath_gsl_func_d_idi(cdf, pascal_P);
rb_camath_gsl_func_d_idi(cdf, pascal_Q);

rb_camath_gsl_rng_d_dd(geometric);
rb_camath_gsl_func_d_id(ran, geometric_pdf);
rb_camath_gsl_func_d_id(cdf, geometric_P);
rb_camath_gsl_func_d_id(cdf, geometric_Q);

rb_camath_gsl_rng_d_diii(hypergeometric);
rb_camath_gsl_func_d_iiii(ran, hypergeometric_pdf);
rb_camath_gsl_func_d_iiii(cdf, hypergeometric_P);
rb_camath_gsl_func_d_iiii(cdf, hypergeometric_Q);

rb_camath_gsl_rng_d_dd(logarithmic);
rb_camath_gsl_func_d_id(ran, logarithmic_pdf);

void
Init_carray_mathfunc_gsl ()
{
  uint32_t seed;

  camath_gsl_rng = gsl_rng_alloc(gsl_rng_mt19937);

  seed = NUM2ULONG(rb_funcall(rb_mKernel,
                              rb_intern("rand"), 1, ULONG2NUM(0xffffffff)));

  gsl_rng_set(camath_gsl_rng, seed);

  /* special functions */

  rb_define_camath_gsl_func2_(airy_Ai, 2);
  rb_define_camath_gsl_func2_(airy_Bi, 2);
  rb_define_camath_gsl_func2_(airy_Ai_scaled, 2);
  rb_define_camath_gsl_func2_(airy_Bi_scaled, 2);
  rb_define_camath_gsl_func2_(airy_Ai_deriv, 2);
  rb_define_camath_gsl_func2_(airy_Bi_deriv, 2);
  rb_define_camath_gsl_func2_(airy_Ai_deriv_scaled, 2);
  rb_define_camath_gsl_func2_(airy_Bi_deriv_scaled, 2);
  rb_define_camath_gsl_func2u(airy_zero_Ai, 1); /* not func2 */
  rb_define_camath_gsl_func2u(airy_zero_Bi, 1); /* not func2 */

  rb_define_camath_gsl_func2(bessel_J0, 1);
  rb_define_camath_gsl_func2(bessel_J1, 1);
  rb_define_camath_gsl_func2(bessel_Jn, 2);

  rb_define_camath_gsl_func2(bessel_Y0, 1);
  rb_define_camath_gsl_func2(bessel_Y1, 1);
  rb_define_camath_gsl_func2(bessel_Yn, 2);

  rb_define_camath_gsl_func2(bessel_I0, 1);
  rb_define_camath_gsl_func2(bessel_I1, 1);
  rb_define_camath_gsl_func2(bessel_In, 2);
  rb_define_camath_gsl_func2(bessel_I0_scaled, 1);
  rb_define_camath_gsl_func2(bessel_I1_scaled, 1);
  rb_define_camath_gsl_func2(bessel_In_scaled, 2);

  rb_define_camath_gsl_func2(bessel_K0, 1);
  rb_define_camath_gsl_func2(bessel_K1, 1);
  rb_define_camath_gsl_func2(bessel_Kn, 2);
  rb_define_camath_gsl_func2(bessel_K0_scaled, 1);
  rb_define_camath_gsl_func2(bessel_K1_scaled, 1);
  rb_define_camath_gsl_func2(bessel_Kn_scaled, 2);

  rb_define_camath_gsl_func2(bessel_j0, 1);
  rb_define_camath_gsl_func2(bessel_j1, 1);
  rb_define_camath_gsl_func2(bessel_j2, 1);
  rb_define_camath_gsl_func2(bessel_jl, 2);

  rb_define_camath_gsl_func2(bessel_y0, 1);
  rb_define_camath_gsl_func2(bessel_y1, 1);
  rb_define_camath_gsl_func2(bessel_y2, 1);
  rb_define_camath_gsl_func2(bessel_yl, 2);

  rb_define_camath_gsl_func2(bessel_i0_scaled, 1);
  rb_define_camath_gsl_func2(bessel_i1_scaled, 1);
  rb_define_camath_gsl_func2(bessel_i2_scaled, 1);
  rb_define_camath_gsl_func2(bessel_il_scaled, 2);

  rb_define_camath_gsl_func2(bessel_k0_scaled, 1);
  rb_define_camath_gsl_func2(bessel_k1_scaled, 1);
  rb_define_camath_gsl_func2(bessel_k2_scaled, 1);
  rb_define_camath_gsl_func2(bessel_kl_scaled, 2);

  rb_define_camath_gsl_func2(bessel_Jnu, 2);
  rb_define_camath_gsl_func2(bessel_Ynu, 2);
  rb_define_camath_gsl_func2(bessel_Inu, 2);
  rb_define_camath_gsl_func2(bessel_Inu_scaled, 2);
  rb_define_camath_gsl_func2(bessel_Knu, 2);
  rb_define_camath_gsl_func2(bessel_Knu_scaled, 2);

  rb_define_camath_gsl_func2u(bessel_zero_J0, 1); /* not func2 */
  rb_define_camath_gsl_func2u(bessel_zero_J1, 1); /* not func2 */
  rb_define_camath_gsl_func2u(bessel_zero_Jnu, 2); /* not func2 */

  rb_define_camath_gsl_func2(gamma, 1);
  rb_define_camath_gsl_func2(lngamma, 1);
  rb_define_camath_gsl_func2(gammastar, 1);
  rb_define_camath_gsl_func2(gammainv, 1);

  rb_define_camath_gsl_func2(gamma_inc, 2);
  rb_define_camath_gsl_func2(gamma_inc_P, 2);
  rb_define_camath_gsl_func2(gamma_inc_Q, 2);

  rb_define_camath_gsl_func2(beta, 2);
  rb_define_camath_gsl_func2(lnbeta, 2);
  rb_define_camath_gsl_func2(beta_inc, 3);

  rb_define_camath_gsl_func2(clausen, 1);

  rb_define_camath_gsl_func2(hydrogenicR_1, 2);
  rb_define_camath_gsl_func2(hydrogenicR, 4);

  rb_define_camath_gsl_func(coupling_3j, 6); /* not func2 */
  rb_define_camath_gsl_func(coupling_6j, 6); /* not func2 */

  rb_define_camath_gsl_func2(dawson, 1);

  rb_define_camath_gsl_func2(debye_1, 1);
  rb_define_camath_gsl_func2(debye_2, 1);
  rb_define_camath_gsl_func2(debye_3, 1);
  rb_define_camath_gsl_func2(debye_4, 1);
  rb_define_camath_gsl_func2(debye_5, 1);
  rb_define_camath_gsl_func2(debye_6, 1);

  rb_define_camath_gsl_func2(dilog, 1);

  rb_define_camath_gsl_func2u(ellint_Kcomp, 2); /* not func2 */
  rb_define_camath_gsl_func2u(ellint_Ecomp, 2); /* not func2 */
  rb_define_camath_gsl_func(ellint_Pcomp, 3); /* not func2 */
  rb_define_camath_gsl_func(ellint_F, 3); /* not func2 */
  rb_define_camath_gsl_func(ellint_E, 3); /* not func2 */
  rb_define_camath_gsl_func(ellint_P, 4); /* not func2 */
  rb_define_camath_gsl_func(ellint_D, 4); /* not func2 */
  rb_define_camath_gsl_func(ellint_RC, 3); /* not func2 */
  rb_define_camath_gsl_func(ellint_RD, 4); /* not func2 */
  rb_define_camath_gsl_func(ellint_RF, 4); /* not func2 */
  rb_define_camath_gsl_func(ellint_RJ, 5); /* not func2 */

  rb_define_camath_gsl_func(elljac, 2); /* not func2 */

  rb_define_camath_gsl_func2(fact, 1);
  rb_define_camath_gsl_func2(doublefact, 1);
  rb_define_camath_gsl_func2(lnfact, 1);
  rb_define_camath_gsl_func2(lndoublefact, 1);
  rb_define_camath_gsl_func2(choose, 2);
  rb_define_camath_gsl_func2(lnchoose, 2);
  rb_define_camath_gsl_func2(taylorcoeff, 2);

  rb_define_camath_gsl_func2(erf, 1);
  rb_define_camath_gsl_func2(erf_Z, 1);
  rb_define_camath_gsl_func2(erf_Q, 1);
  rb_define_camath_gsl_func2(erfc, 1);
  rb_define_camath_gsl_func2(log_erfc, 1);
  rb_define_camath_gsl_func2(hazard, 1);

  rb_define_camath_gsl_func2(expint_E1, 1);
  rb_define_camath_gsl_func2(expint_E2, 1);
  rb_define_camath_gsl_func2(expint_En, 2);
  rb_define_camath_gsl_func2(expint_Ei, 1);
  rb_define_camath_gsl_func2(Shi, 1);
  rb_define_camath_gsl_func2(Chi, 1);
  rb_define_camath_gsl_func2(expint_3, 1);
  rb_define_camath_gsl_func2(Si, 1);
  rb_define_camath_gsl_func2(Ci, 1);
  rb_define_camath_gsl_func2(atanint, 1);

  rb_define_camath_gsl_func2(fermi_dirac_m1, 1);
  rb_define_camath_gsl_func2(fermi_dirac_0, 1);
  rb_define_camath_gsl_func2(fermi_dirac_1, 1);
  rb_define_camath_gsl_func2(fermi_dirac_2, 1);
  rb_define_camath_gsl_func2(fermi_dirac_int, 2);
  rb_define_camath_gsl_func2(fermi_dirac_mhalf, 1);
  rb_define_camath_gsl_func2(fermi_dirac_half, 1);
  rb_define_camath_gsl_func2(fermi_dirac_3half, 1);
  rb_define_camath_gsl_func2(fermi_dirac_inc_0, 2);

  rb_define_camath_gsl_func2(gegenpoly_1, 2);
  rb_define_camath_gsl_func2(gegenpoly_2, 2);
  rb_define_camath_gsl_func2(gegenpoly_3, 2);
  rb_define_camath_gsl_func2(gegenpoly_n, 3);

  rb_define_camath_gsl_func2(hyperg_0F1, 2);
  rb_define_camath_gsl_func2(hyperg_1F1_int, 3);
  rb_define_camath_gsl_func2(hyperg_1F1, 3);
  rb_define_camath_gsl_func2(hyperg_U_int, 3);
  rb_define_camath_gsl_func2(hyperg_U, 3);
  rb_define_camath_gsl_func2(hyperg_2F1, 4);
  rb_define_camath_gsl_func2(hyperg_2F1_conj, 4);
  rb_define_camath_gsl_func2(hyperg_2F1_renorm, 4);
  rb_define_camath_gsl_func2(hyperg_2F1_conj_renorm, 4);
  rb_define_camath_gsl_func2(hyperg_2F0, 3);

  rb_define_camath_gsl_func2(laguerre_1, 2);
  rb_define_camath_gsl_func2(laguerre_2, 2);
  rb_define_camath_gsl_func2(laguerre_3, 2);
  rb_define_camath_gsl_func2(laguerre_n, 3);

  rb_define_camath_gsl_func2(lambert_W0, 1);
  rb_define_camath_gsl_func2(lambert_Wm1, 1);

  rb_define_camath_gsl_func2(legendre_P1, 1);
  rb_define_camath_gsl_func2(legendre_P2, 1);
  rb_define_camath_gsl_func2(legendre_P3, 1);
  rb_define_camath_gsl_func2(legendre_Pl, 2);

  rb_define_camath_gsl_func2(legendre_Q0, 1);
  rb_define_camath_gsl_func2(legendre_Q1, 1);
  rb_define_camath_gsl_func2(legendre_Ql, 2);

  rb_define_camath_gsl_func2(legendre_Plm, 3);
  rb_define_camath_gsl_func2(legendre_sphPlm, 3);

  rb_define_camath_gsl_func2(conicalP_half, 2);
  rb_define_camath_gsl_func2(conicalP_mhalf, 2);
  rb_define_camath_gsl_func2(conicalP_0, 2);
  rb_define_camath_gsl_func2(conicalP_1, 2);
  rb_define_camath_gsl_func2(conicalP_sph_reg, 3);
  rb_define_camath_gsl_func2(conicalP_cyl_reg, 3);

  rb_define_camath_gsl_func2(legendre_H3d_0, 2);
  rb_define_camath_gsl_func2(legendre_H3d_1, 2);
  rb_define_camath_gsl_func2(legendre_H3d, 3);

  rb_define_camath_gsl_func2(psi, 1);
  rb_define_camath_gsl_func2(psi_1piy, 1);
  rb_define_camath_gsl_func2(psi_1, 1);
  rb_define_camath_gsl_func2(psi_n, 2);

  rb_define_camath_gsl_func2(synchrotron_1, 1);
  rb_define_camath_gsl_func2(synchrotron_2, 1);

  rb_define_camath_gsl_func2(transport_2, 1);
  rb_define_camath_gsl_func2(transport_3, 1);
  rb_define_camath_gsl_func2(transport_4, 1);
  rb_define_camath_gsl_func2(transport_5, 1);

  rb_define_camath_gsl_func2(sinc, 1);
  rb_define_camath_gsl_func(polar_to_rect, 2); /* not func2 */
  rb_define_camath_gsl_func(rect_to_polar, 2); /* not func2 */
  rb_define_camath_gsl_func2(angle_restrict_symm, 1);
  rb_define_camath_gsl_func2(angle_restrict_pos, 1);

  rb_define_camath_gsl_func2(zeta, 1);
  rb_define_camath_gsl_func2(zetam1, 1);
  rb_define_camath_gsl_func2(hzeta, 2);
  rb_define_camath_gsl_func2(eta, 1);

  /* random distributions */

  rb_define_camath_gsl_func(random_gaussian, 2);
  rb_define_camath_gsl_func2(gaussian_pdf, 2);
  rb_define_camath_gsl_func2(gaussian_P, 2);
  rb_define_camath_gsl_func2(gaussian_Q, 2);
  rb_define_camath_gsl_func2(gaussian_Pinv, 2);
  rb_define_camath_gsl_func2(gaussian_Qinv, 2);

  rb_define_camath_gsl_func(random_ugaussian, 1);
  rb_define_camath_gsl_func2(ugaussian_pdf, 1);
  rb_define_camath_gsl_func2(ugaussian_P, 1);
  rb_define_camath_gsl_func2(ugaussian_Q, 1);
  rb_define_camath_gsl_func2(ugaussian_Pinv, 1);
  rb_define_camath_gsl_func2(ugaussian_Qinv, 1);

  rb_define_camath_gsl_func(random_gaussian_tail, 3);
  rb_define_camath_gsl_func2(gaussian_tail_pdf, 3);

  rb_define_camath_gsl_func(random_ugaussian_tail, 2);
  rb_define_camath_gsl_func2(ugaussian_tail_pdf, 2);

  rb_define_camath_gsl_func(random_bivariate_gaussian, 3);
  rb_define_camath_gsl_func(bivariate_gaussian_pdf, 5);   /* not func2 */

  rb_define_camath_gsl_func(random_exponential, 2);
  rb_define_camath_gsl_func2(exponential_pdf, 2);
  rb_define_camath_gsl_func2(exponential_P, 2);
  rb_define_camath_gsl_func2(exponential_Q, 2);
  rb_define_camath_gsl_func2(exponential_Pinv, 2);
  rb_define_camath_gsl_func2(exponential_Qinv, 2);

  rb_define_camath_gsl_func(random_laplace, 2);
  rb_define_camath_gsl_func2(laplace_pdf, 2);
  rb_define_camath_gsl_func2(laplace_P, 2);
  rb_define_camath_gsl_func2(laplace_Q, 2);
  rb_define_camath_gsl_func2(laplace_Pinv, 2);
  rb_define_camath_gsl_func2(laplace_Qinv, 2);

  rb_define_camath_gsl_func(random_exppow, 3);
  rb_define_camath_gsl_func2(exppow_pdf, 3);
  rb_define_camath_gsl_func2(exppow_P, 3);
  rb_define_camath_gsl_func2(exppow_Q, 3);

  rb_define_camath_gsl_func(random_cauchy, 2);
  rb_define_camath_gsl_func2(cauchy_pdf, 2);
  rb_define_camath_gsl_func2(cauchy_P, 2);
  rb_define_camath_gsl_func2(cauchy_Q, 2);
  rb_define_camath_gsl_func2(cauchy_Pinv, 2);
  rb_define_camath_gsl_func2(cauchy_Qinv, 2);

  rb_define_camath_gsl_func(random_rayleigh, 2);
  rb_define_camath_gsl_func2(rayleigh_pdf, 2);
  rb_define_camath_gsl_func2(rayleigh_P, 2);
  rb_define_camath_gsl_func2(rayleigh_Q, 2);
  rb_define_camath_gsl_func2(rayleigh_Pinv, 2);
  rb_define_camath_gsl_func2(rayleigh_Qinv, 2);

  rb_define_camath_gsl_func(random_rayleigh_tail, 3);
  rb_define_camath_gsl_func2(rayleigh_tail_pdf, 3);

  rb_define_camath_gsl_func(random_landau, 1);
  rb_define_camath_gsl_func2(landau_pdf, 1);
  
  rb_define_camath_gsl_func(random_levy, 3);
  rb_define_camath_gsl_func(random_levy_skew, 4);

  rb_define_camath_gsl_func(random_gamma, 3);
  rb_define_camath_gsl_func2(gamma_pdf, 3);
  rb_define_camath_gsl_func2(gamma_P, 3);
  rb_define_camath_gsl_func2(gamma_Q, 3);
  rb_define_camath_gsl_func2(gamma_Pinv, 3);
  rb_define_camath_gsl_func2(gamma_Qinv, 3);

  rb_define_camath_gsl_func(random_erlang, 3);
  rb_define_camath_gsl_func2(erlang_pdf, 3);

  rb_define_camath_gsl_func(random_flat, 3);
  rb_define_camath_gsl_func2(flat_pdf, 3);
  rb_define_camath_gsl_func2(flat_P, 3);
  rb_define_camath_gsl_func2(flat_Q, 3);
  rb_define_camath_gsl_func2(flat_Pinv, 3);
  rb_define_camath_gsl_func2(flat_Qinv, 3);

  rb_define_camath_gsl_func(random_lognormal, 3);
  rb_define_camath_gsl_func2(lognormal_pdf, 3);
  rb_define_camath_gsl_func2(lognormal_P, 3);
  rb_define_camath_gsl_func2(lognormal_Q, 3);
  rb_define_camath_gsl_func2(lognormal_Pinv, 3);
  rb_define_camath_gsl_func2(lognormal_Qinv, 3);

  rb_define_camath_gsl_func(random_chisq, 2);
  rb_define_camath_gsl_func2(chisq_pdf, 2);
  rb_define_camath_gsl_func2(chisq_P, 2);
  rb_define_camath_gsl_func2(chisq_Q, 2);
  rb_define_camath_gsl_func2(chisq_Pinv, 2);
  rb_define_camath_gsl_func2(chisq_Qinv, 2);

  rb_define_camath_gsl_func(random_fdist, 3);
  rb_define_camath_gsl_func2(fdist_pdf, 3);
  rb_define_camath_gsl_func2(fdist_P, 3);
  rb_define_camath_gsl_func2(fdist_Q, 3);
  rb_define_camath_gsl_func2(fdist_Pinv, 3);
  rb_define_camath_gsl_func2(fdist_Qinv, 3);

  rb_define_camath_gsl_func(random_tdist, 2);
  rb_define_camath_gsl_func2(tdist_pdf, 2);
  rb_define_camath_gsl_func2(tdist_P, 2);
  rb_define_camath_gsl_func2(tdist_Q, 2);
  rb_define_camath_gsl_func2(tdist_Pinv, 2);
  rb_define_camath_gsl_func2(tdist_Qinv, 2);
  
  rb_define_camath_gsl_func(random_beta, 3);
  rb_define_camath_gsl_func2(beta_pdf, 3);
  rb_define_camath_gsl_func2(beta_P, 3);
  rb_define_camath_gsl_func2(beta_Q, 3);
  rb_define_camath_gsl_func2(beta_Pinv, 3);
  rb_define_camath_gsl_func2(beta_Qinv, 3);

  rb_define_camath_gsl_func(random_logistic, 2);
  rb_define_camath_gsl_func2(logistic_pdf, 2);
  rb_define_camath_gsl_func2(logistic_P, 2);
  rb_define_camath_gsl_func2(logistic_Q, 2);
  rb_define_camath_gsl_func2(logistic_Pinv, 2);
  rb_define_camath_gsl_func2(logistic_Qinv, 2);
  
  rb_define_camath_gsl_func(random_pareto, 3);
  rb_define_camath_gsl_func2(pareto_pdf, 3);
  rb_define_camath_gsl_func2(pareto_P, 3);
  rb_define_camath_gsl_func2(pareto_Q, 3);
  rb_define_camath_gsl_func2(pareto_Pinv, 3);
  rb_define_camath_gsl_func2(pareto_Qinv, 3);

  rb_define_camath_gsl_func(random_dir_2d, 1); /* ran */
  rb_define_camath_gsl_func(random_dir_3d, 1); /* ran */
  rb_define_camath_gsl_func(random_dir_nd, 2); /* ran */

  rb_define_camath_gsl_func(random_weibull, 3);
  rb_define_camath_gsl_func2(weibull_pdf, 3);
  rb_define_camath_gsl_func2(weibull_P, 3);
  rb_define_camath_gsl_func2(weibull_Q, 3);
  rb_define_camath_gsl_func2(weibull_Pinv, 3);
  rb_define_camath_gsl_func2(weibull_Qinv, 3);

  rb_define_camath_gsl_func(random_gumbel1, 3);
  rb_define_camath_gsl_func2(gumbel1_pdf, 3);
  rb_define_camath_gsl_func2(gumbel1_P, 3);
  rb_define_camath_gsl_func2(gumbel1_Q, 3);
  rb_define_camath_gsl_func2(gumbel1_Pinv, 3);
  rb_define_camath_gsl_func2(gumbel1_Qinv, 3);

  rb_define_camath_gsl_func(random_gumbel2, 3);
  rb_define_camath_gsl_func2(gumbel2_pdf, 3);
  rb_define_camath_gsl_func2(gumbel2_P, 3);
  rb_define_camath_gsl_func2(gumbel2_Q, 3);
  rb_define_camath_gsl_func2(gumbel2_Pinv, 3);
  rb_define_camath_gsl_func2(gumbel2_Qinv, 3);

/*
  rb_define_camath_gsl_func(random_dirichlet, -1); */

/*
  rb_define_camath_gsl_func(random_discrete_preproc, -1); */

  rb_define_camath_gsl_func(random_poisson, 2);
  rb_define_camath_gsl_func2(poisson_pdf, 2);
  rb_define_camath_gsl_func2(poisson_P, 2); 
  rb_define_camath_gsl_func2(poisson_Q, 2); 

  rb_define_camath_gsl_func(random_bernoulli, 2);
  rb_define_camath_gsl_func2u(bernoulli_pdf, 2); /* not func2 */
  
  rb_define_camath_gsl_func(random_binomial, 3);
  rb_define_camath_gsl_func(binomial_pdf, 3); /* not func2 */
  rb_define_camath_gsl_func(binomial_P, 3); /* not func2 */
  rb_define_camath_gsl_func(binomial_Q, 3); /* not func2 */

  rb_define_camath_gsl_func(random_multinomial, 3); 
  rb_define_camath_gsl_func(multinomial_pdf, 2);/* not func2 */
  rb_define_camath_gsl_func(multinomial_lnpdf, 2); /* not func2 */

  rb_define_camath_gsl_func(random_negative_binomial, 3);
  rb_define_camath_gsl_func2u(negative_binomial_pdf, 3); /* not func2 */
  rb_define_camath_gsl_func2u(negative_binomial_P, 3); /* not func2 */
  rb_define_camath_gsl_func2u(negative_binomial_Q, 3); /* not func2 */

  rb_define_camath_gsl_func(random_pascal, 3);
  rb_define_camath_gsl_func(pascal_pdf, 3); /* not func2 */
  rb_define_camath_gsl_func(pascal_P, 3); /* not func2 */
  rb_define_camath_gsl_func(pascal_Q, 3); /* not func2 */

  rb_define_camath_gsl_func(random_geometric, 2);
  rb_define_camath_gsl_func2u(geometric_pdf, 2); /* not func2 */
  rb_define_camath_gsl_func2u(geometric_P, 2); /* not func2 */
  rb_define_camath_gsl_func2u(geometric_Q, 2); /* not func2 */

  rb_define_camath_gsl_func(random_hypergeometric, 4);
  rb_define_camath_gsl_func(hypergeometric_pdf, 4); /* not func2 */
  rb_define_camath_gsl_func(hypergeometric_P, 4); /* not func2 */
  rb_define_camath_gsl_func(hypergeometric_Q, 4); /* not func2 */

  rb_define_camath_gsl_func(random_logarithmic, 2);
  rb_define_camath_gsl_func2u(logarithmic_pdf, 2); /* not func2 */

}
