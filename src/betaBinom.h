/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998-2011  The R Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  http://www.r-project.org/Licenses/
 */

/* Private header file for use during compilation of Mathlib */

#include "R_ext/Visibility.h"
#define DEBUG_c_qbetabinom 0 // Print debug info in C


double	attribute_hidden dbetabinom_raw(double, double, double, double, int);


/* Beta Binomial Distribution */
double	dbetabinom(double, double, double, double, int);
double	pbetabinom(double, double, double, double, int, int);
double	qbetabinom(double, double, double, double, int, int);
double	rbetabinom(double, double, double);
