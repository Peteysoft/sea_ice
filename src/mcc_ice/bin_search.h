//
// Copyright 2004 Peter Mills.  All rights reserved.
//
// This is the header file for "peteys_tmpl_lib.cc."  Only templates are
// given here.  The desired instantiation may or may not exist--check the
// main source file.

#ifndef TMPL_LIB_INCLUDED
#define TMPL_LIB_INCLUDED

#pragma interface

//uses a binary search to search an ordered list:
template <class dt>
long bin_search (dt *list, long n, dt value, long last_ind=-1);

template <class dt>
double interpolate (dt *list, long n, dt value, long last_ind=-1);

#endif
