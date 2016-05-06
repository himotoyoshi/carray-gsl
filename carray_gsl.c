/* ---------------------------------------------------------------------------

  carray/carray_gsl.c

  This file is part of Ruby/CArray extension library.
  You can redistribute it and/or modify it under the terms of
  the GNU General Public License (GPL).

  Copyright (C) 2005-2008  Hiroki Motoyoshi

---------------------------------------------------------------------------- */

#include "carray.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

static VALUE rb_mGSL;
static VALUE rb_cVector, rb_cVectorInt, rb_cVectorComplex;
static VALUE rb_cVectorCol, rb_cVectorIntCol, rb_cVectorComplexCol;
static VALUE rb_cMatrix, rb_cMatrixInt, rb_cMatrixComplex;
static VALUE rb_cCAMatrix, rb_cCAVector;

static VALUE
rb_gsl_vector_ca (VALUE self)
{
  volatile VALUE obj, block;
  gsl_vector *v;
  int data_type;
  int32_t dim[2], dim0;

  Data_Get_Struct(self, gsl_vector, v);

  if ( rb_obj_is_kind_of(self, rb_cVectorComplex) ) {
    data_type = CA_CMPLX128;
  }
  else if ( rb_obj_is_kind_of(self, rb_cVectorInt) ) {
    data_type = CA_INT32;
  }
  else if ( rb_obj_is_kind_of(self, rb_cVector) ) {
    data_type = CA_FLOAT64;
  }
  else {
    rb_raise(rb_eRuntimeError, "can't create CArray from GSL::Vector");
  }

  if ( v->stride > 1 ) {
    int32_t start0, step0, count0;
    if ( rb_obj_is_kind_of(self, rb_cVectorCol) ||
         rb_obj_is_kind_of(self, rb_cVectorIntCol) ||
         rb_obj_is_kind_of(self, rb_cVectorComplexCol) ) {
      dim0 = dim[0] = v->size / ca_sizeof[data_type];
      dim[1] = 1;
      block = rb_carray_wrap_ptr(data_type, 2, dim, 0, NULL, 
                                 (char *) v->data, self);
    }
    else {
      dim0 = dim[0] = v->size / ca_sizeof[data_type];
      block = rb_carray_wrap_ptr(data_type, 1, dim, 0, NULL, 
                                 (char *) v->data, self);
    }
    start0 = 0;
    step0  = v->stride;
    count0 = v->size;
    obj = rb_ca_block_new(block, 1, &dim0, &start0, &step0, &count0, 0);
  }
  else {
    if ( rb_obj_is_kind_of(self, rb_cVectorCol) ||
         rb_obj_is_kind_of(self, rb_cVectorIntCol) ||
         rb_obj_is_kind_of(self, rb_cVectorComplexCol) ) {
      dim[0] = v->size;
      dim[1] = 1;
      obj = rb_carray_wrap_ptr(data_type, 2, dim, 0, NULL, 
                                 (char *) v->data, self);
    }
    else {
      dim[0] = v->size;
      obj = rb_carray_wrap_ptr(data_type, 1, dim, 0, NULL, 
                                 (char *) v->data, self);
    }
  }

  return obj;
}

static VALUE
rb_ca_gv_i (VALUE self)
{
  volatile VALUE klass, obj;
  CArray *ca;

  Data_Get_Struct(self, CArray, ca);

  if ( ! ca_is_attached(ca) ) {
    rb_raise(rb_eRuntimeError, 
             "can't create GSL vector refers detached array");
  }

  switch ( ca->data_type ) {
  case CA_INT32: {
    gsl_vector_int *v;
    if ( ( ca->rank == 2 ) && ( ca->dim[1] == 1 ) ) {
      klass = rb_cVectorIntCol;      
    }
    else {
      klass = rb_cVectorInt;
    }
    v = ALLOC(gsl_vector_int);
    v->data   = (int *) ca_ptr_at_addr(ca, 0);
    v->size   = ca->elements;
    v->stride = 1;
    v->block  = NULL;
    v->owner  = 0;
    obj = Data_Wrap_Struct(klass, 0, gsl_vector_int_free, v);
    break;
  }
  case CA_FLOAT64: {
    gsl_vector *v;
    if ( ( ca->rank == 2 ) && ( ca->dim[1] == 1 ) ) {
      klass = rb_cVectorCol;      
    }
    else {
      klass = rb_cVector;
    }
    v = ALLOC(gsl_vector);
    v->data   = (double*) ca_ptr_at_addr(ca, 0);
    v->size   = ca->elements;
    v->stride = 1;
    v->block  = NULL;
    v->owner  = 0;
    obj = Data_Wrap_Struct(klass, 0, gsl_vector_free, v);
    break;
  }
  case CA_CMPLX128: {
    gsl_vector_complex *v;
    if ( ( ca->rank == 2 ) && ( ca->dim[1] == 1 ) ) {
      klass = rb_cVectorComplexCol;      
    }
    else {
      klass = rb_cVectorComplex;
    }
    v = ALLOC(gsl_vector_complex);
    v->data   = (double *) ca_ptr_at_addr(ca, 0);
    v->size   = 2*ca->elements;
    v->stride = 1;
    v->block  = NULL;
    v->owner  = 0;
    obj = Data_Wrap_Struct(klass, 0, gsl_vector_complex_free, v);
    break;
  }
  default:
    rb_raise(rb_eRuntimeError, "invalid data_type");
  }

  rb_ivar_set(obj, rb_intern("referred_object"), self);

  return obj;
}

static VALUE
rb_ca_gv (VALUE self)
{
  return rb_ca_gv_i(self);
}

static VALUE
rb_gsl_matrix_ca (VALUE self)
{
  volatile VALUE obj, block;
  gsl_matrix *v;
  int data_type;
  int32_t dim[2];

  Data_Get_Struct(self, gsl_matrix, v);

  if ( rb_obj_is_kind_of(self, rb_cMatrixComplex) ) {
    data_type = CA_CMPLX128;
  }
  else if ( rb_obj_is_kind_of(self, rb_cMatrix) ) {
    data_type = CA_FLOAT64;
  }
  else if ( rb_obj_is_kind_of(self, rb_cMatrixInt) ) {
    data_type = CA_INT32;
  }
  else {
    rb_raise(rb_eRuntimeError, "can't create CArray from GSL::Matrix");    
  }

  if ( v->tda != v->size2 ) {
    int32_t start[2], step[2], count[2];
    dim[0] = v->size1;
    dim[1] = v->tda;
    block = rb_carray_wrap_ptr(data_type, 2, dim, 0, NULL, 
                                           (char *) v->data, self);
    start[0] = 0;
    start[1] = 0;
    step[0]  = 1;
    step[1]  = 1;
    count[0] = v->size1;
    count[1] = v->size2;
    obj = rb_ca_block_new(block, 2, dim, start, step, count, 0);
  }
  else {
    dim[0] = v->size1;
    dim[1] = v->size2;
    obj = rb_carray_wrap_ptr(data_type, 2, dim, 0, NULL, 
                             (char *) v->data, self);
  }

  return obj;
}

static VALUE
rb_ca_gm_i (VALUE self)
{
  volatile VALUE klass, obj;
  CArray *ca;

  Data_Get_Struct(self, CArray, ca);

  if ( ! ca_is_attached(ca) ) {
    rb_raise(rb_eRuntimeError, 
             "can't create GSL matrix for detached array");
  }

  if ( ca->rank < 2 && ca->dim[0] == 0 ) {
    rb_raise(rb_eRuntimeError, "invalid rank or dim to convert GSL::Matrix");
  }

  switch ( ca->data_type ) {
  case CA_INT32: {
    gsl_matrix_int *v;
    klass = rb_cMatrixInt;
    v = ALLOC(gsl_matrix_int);
    v->data   = (int *) ca_ptr_at_addr(ca, 0);
    v->size1  = ca->dim[0];
    v->size2  = ca->elements/ca->dim[0];
    v->tda    = v->size2;
    v->block  = NULL;
    v->owner  = 0;
    obj = Data_Wrap_Struct(klass, 0, gsl_matrix_int_free, v);
    break;
  }
  case CA_FLOAT64: {
    gsl_matrix *v;    
    klass = rb_cMatrix;
    v = ALLOC(gsl_matrix);
    v->data   = (double *) ca_ptr_at_addr(ca, 0);
    v->size1  = ca->dim[0];
    v->size2  = ca->elements/ca->dim[0];
    v->tda    = v->size2;
    v->block  = NULL;
    v->owner  = 0;
    obj = Data_Wrap_Struct(klass, 0, gsl_matrix_free, v);
    break;
  }
  case CA_CMPLX128: {
    gsl_matrix_complex *v;    
    klass = rb_cMatrixComplex;
    v = ALLOC(gsl_matrix_complex);
    v->data   = (double *) ca_ptr_at_addr(ca, 0);
    v->size1  = ca->dim[0];
    v->size2  = 2*ca->elements/ca->dim[0];
    v->tda    = v->size2;
    v->block  = NULL;
    v->owner  = 0;
    obj = Data_Wrap_Struct(klass, 0, gsl_matrix_complex_free, v);
    break;
  }
  default:
    rb_raise(rb_eRuntimeError, "invalid data_type");
  }

  rb_ivar_set(obj, rb_intern("referred_object"), self);

  return obj;
}

static VALUE
rb_ca_gm (VALUE self)
{
  return rb_ca_gm_i(self);
}

static VALUE
rb_ca_to_gv (VALUE self)
{
  volatile VALUE klass, obj;
  CArray *ca;

  Data_Get_Struct(self, CArray, ca);

  ca_attach(ca);

  switch ( ca->data_type ) {
  case CA_INT32: {
    gsl_vector_int *v;
    if ( ( ca->rank == 2 ) && ( ca->dim[1] == 1 ) ) {
      klass = rb_cVectorIntCol;      
    }
    else {
      klass = rb_cVectorInt;
    }
    v = gsl_vector_int_alloc(ca->elements);
    memcpy(v->data, ca->ptr, v->size*sizeof(int));
    obj = Data_Wrap_Struct(klass, 0, gsl_vector_int_free, v);
    break;
  }
  case CA_FLOAT64: {
    gsl_vector *v;
    if ( ( ca->rank == 2 ) && ( ca->dim[1] == 1 ) ) {
      klass = rb_cVectorCol;      
    }
    else {
      klass = rb_cVector;
    }
    v = gsl_vector_alloc(ca->elements);
    memcpy(v->data, ca->ptr, v->size*sizeof(double));
    obj = Data_Wrap_Struct(klass, 0, gsl_vector_free, v);
    break;
  }
  case CA_CMPLX128: {
    gsl_vector_complex *v;
    if ( ( ca->rank == 2 ) && ( ca->dim[1] == 1 ) ) {
      klass = rb_cVectorComplexCol;      
    }
    else {
      klass = rb_cVectorComplex;
    }
    v = gsl_vector_complex_alloc(ca->elements);
    memcpy(v->data, ca->ptr, v->size*sizeof(gsl_complex));
    obj = Data_Wrap_Struct(klass, 0, gsl_vector_complex_free, v);
    break;
  }
  default:
    rb_raise(rb_eRuntimeError, "invalid data_type");
  }

  ca_detach(ca);

  return obj;
}

static VALUE
rb_ca_to_gm (VALUE self)
{
  volatile VALUE klass, obj;
  CArray *ca;

  Data_Get_Struct(self, CArray, ca);

  if ( ca->rank < 2 && ca->dim[0] == 0 ) {
    rb_raise(rb_eRuntimeError, "invalid rank or dim to convert GSL::Matrix");
  }

  ca_attach(ca);

  switch ( ca->data_type ) {
  case CA_INT32: {
    gsl_matrix_int *m;
    klass = rb_cMatrixInt;
    m = gsl_matrix_int_alloc(ca->dim[0], ca->elements/ca->dim[0]);
    memcpy(m->data, ca->ptr, m->size1*m->size2*sizeof(int));
    obj = Data_Wrap_Struct(klass, 0, gsl_matrix_int_free, m);
    break;
  }
  case CA_FLOAT64: {
    gsl_matrix *m;
    klass = rb_cMatrix;
    m = gsl_matrix_alloc(ca->dim[0], ca->elements/ca->dim[0]);
    memcpy(m->data, ca->ptr, m->size1*m->size2*sizeof(double));
    obj = Data_Wrap_Struct(klass, 0, gsl_matrix_free, m);
    break;
  }
  case CA_CMPLX128: {
    gsl_matrix_complex *m;
    klass = rb_cMatrixComplex;
    m = gsl_matrix_complex_alloc(ca->dim[0], ca->elements/ca->dim[0]);
    memcpy(m->data, ca->ptr, m->size1*m->size2*sizeof(gsl_complex));
    obj = Data_Wrap_Struct(klass, 0, gsl_matrix_complex_free, m);
    break;
  }
  default:
    rb_raise(rb_eRuntimeError, "invalid data_type");
  }

  ca_detach(ca);

  return obj;
}

static VALUE
rb_ca_m (VALUE self)
{
  volatile VALUE obj;
  CArray *ca;
  Data_Get_Struct(self, CArray, ca);
  if ( rb_obj_is_kind_of(self, rb_cCAMatrix) ) {
    return self;
  }
  if ( ca->data_type == CA_FIXLEN ) {
    rb_raise(rb_eRuntimeError, "can't create fixlen type CAMatrix");
  }
  if ( ca->rank != 2 ) {
    rb_raise(rb_eRuntimeError, "rank should be 2 for CAMatrix");
  }
  /* super advanced technique
     to fake other array including virtual array as CAWrap object.
     When attach, detach ... are called, the routins for the 
     original array are executed. 
     For these, CAMatrix should be a subclass of CAWrap. */
  obj = Data_Wrap_Struct(rb_cCAMatrix, NULL, NULL, ca);
  rb_ivar_set(obj, rb_intern("referred_object"), self);
  return obj;
}

static VALUE
rb_ca_v (VALUE self)
{
  volatile VALUE obj;
  CArray *ca;
  Data_Get_Struct(self, CArray, ca);
  if ( rb_obj_is_kind_of(self, rb_cCAVector) ) {
    return self;
  }
  if ( ca->data_type == CA_FIXLEN ) {
    rb_raise(rb_eRuntimeError, "can't create fixlen type CAMatrix");
  }
  /* super advanced technique 
     to fake other array including virtual array as CAWrap object.
     When attach, detach ... are called, the routins for the 
     original array are executed.
     For these, CAVector should be a subclass of CAWrap. */
  if ( ca->rank == 1 ) {
    obj = Data_Wrap_Struct(rb_cCAVector, NULL, NULL, ca);
  }
  else if ( ca->rank == 2 && ( ca->dim[0] == 1 || ca->dim[1] == 1 ) ) {
    obj = Data_Wrap_Struct(rb_cCAVector, NULL, NULL, ca);
  }
  else {
    rb_raise(rb_eRuntimeError, "rank should be 1 for CAVector");
  }
  rb_ivar_set(obj, rb_intern("referred_object"), self);
  return obj;
}

static VALUE
rb_ca_a (VALUE self)
{
  if ( rb_obj_is_kind_of(self, rb_cCAMatrix) ||
       rb_obj_is_kind_of(self, rb_cCAVector) ) {
     return rb_ivar_get(self, rb_intern("referred_object"));         
  }
  else {
    return self;
  }
}

void Init_carray_mathfunc_gsl ();

void
Init_carray_gsl ()
{
  rb_mGSL           = rb_const_get(rb_cObject,    rb_intern("GSL"));
  rb_cVector        = rb_const_get(rb_mGSL,       rb_intern("Vector"));
  rb_cVectorInt     = rb_const_get(rb_cVector,    rb_intern("Int"));
  rb_cVectorComplex = rb_const_get(rb_cVector,    rb_intern("Complex"));
  rb_cVectorCol     = rb_const_get(rb_cVector,    rb_intern("Col"));
  rb_cVectorIntCol  = rb_const_get(rb_cVectorInt, rb_intern("Col"));
  rb_cVectorComplexCol = rb_const_get(rb_cVectorComplex, rb_intern("Col"));
  rb_cMatrix        = rb_const_get(rb_mGSL,       rb_intern("Matrix"));
  rb_cMatrixInt     = rb_const_get(rb_cMatrix,    rb_intern("Int"));
  rb_cMatrixComplex = rb_const_get(rb_cMatrix,    rb_intern("Complex"));

  rb_define_method(rb_cVector,           "ca", rb_gsl_vector_ca, 0);
  rb_define_method(rb_cVectorInt,        "ca", rb_gsl_vector_ca, 0);
  rb_define_method(rb_cVectorComplex,    "ca", rb_gsl_vector_ca, 0);
  rb_define_method(rb_cVectorCol,        "ca", rb_gsl_vector_ca, 0);
  rb_define_method(rb_cVectorIntCol,     "ca", rb_gsl_vector_ca, 0);
  rb_define_method(rb_cVectorComplexCol, "ca", rb_gsl_vector_ca, 0);

  rb_define_method(rb_cMatrix,           "ca", rb_gsl_matrix_ca, 0);
  rb_define_method(rb_cMatrixInt,        "ca", rb_gsl_matrix_ca, 0);
  rb_define_method(rb_cMatrixComplex,    "ca", rb_gsl_matrix_ca, 0);

  rb_define_method(rb_cCArray, "gv",    rb_ca_gv, 0);
  rb_define_method(rb_cCArray, "to_gv", rb_ca_to_gv, 0);

  rb_define_method(rb_cCArray, "gm",    rb_ca_gm, 0);
  rb_define_method(rb_cCArray, "to_gm", rb_ca_to_gm, 0);

  rb_cCAMatrix = rb_define_class("CAMatrix", rb_cCAWrap); /* see rb_ca_m() */
  rb_cCAVector = rb_define_class("CAVector", rb_cCAWrap); /* see rb_ca_v() */
  rb_define_method(rb_cCArray, "m", rb_ca_m, 0);  
  rb_define_method(rb_cCArray, "v", rb_ca_v, 0);  
  rb_define_method(rb_cCArray, "a", rb_ca_a, 0);  

  Init_carray_mathfunc_gsl();

  rb_define_const(rb_cCArray, "HAVE_GSL", Qtrue);
}

