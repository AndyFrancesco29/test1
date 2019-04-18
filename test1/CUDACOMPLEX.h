//only consider operations between DataType and CudaComplex
//https://stackoverflow.com/questions/9860711/cucomplex-h-and-exp
//add a method to generate a random complex numbr between as cos(N)+jsin (N), N in [0, 2pi]


/*
Copyright (C) 1994 Free Software Foundation
	written by Guido Perrone (perrone@polito.it)
This file is part of the GPClass portable library. This library is free
software; you can redistribute it and/or modify it under the terms of
the GNU Library General Public License as published by the Free
Software Foundation; either version 2 of the License, or (at your
option) any later version.  This library is distributed in the hope
that it will be useful, but WITHOUT ANY WARRANTY; without even the
implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the GNU Library General Public License for more details.
You should have received a copy of the GNU Library General Public
License along with this library; if not, write to the Free Software
Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
*/




#if !defined( __CudaComplex_H )
#define __CudaComplex_H


#if !defined(__MATH_H)
#include <math.h>
#endif
#if !defined( __IOSTREAM_H )
#include <iostream>
#endif  // __IOSTREAM_H


//#if !defined(__gERR_H)
//#include "gerr.h"
//#endif  // __gERR_H


typedef float DataType;
DataType ata(DataType x, DataType y);

class CudaComplex
{

public:
	 CudaComplex();																// default constructor
	 CudaComplex(const DataType x, const DataType y = 0.);						// constr w/ number
	 CudaComplex(const CudaComplex &other);										// copy constr
	 ~CudaComplex();																// destructor

// -------- comparison -----------
	 bool equal(const DataType x, const DataType y = 0.) const;					// equality
	 bool operator==(const CudaComplex &other) const;							// equal operator
	 const CudaComplex& operator +(void) const;
	 const CudaComplex operator -(void) const;
	 bool diff(const DataType x, const DataType y = 0.) const;					// difference
	 bool operator!=(const CudaComplex &other) const;							// diff operator

// -------- assignement -----------
	 void copy(const DataType x, const DataType y = 0.);							// copies 2 number in a CudaComplex
	 CudaComplex &operator=(const CudaComplex &other);							// assignment operator

// -------- addition -----------
	 void add(const CudaComplex &oth1, const CudaComplex &oth2);										// adds 2 CudaComplex
	 friend CudaComplex operator + (const CudaComplex &oth1, const CudaComplex &oth2);				// sum of CudaComplex   
	 friend CudaComplex operator + (const DataType oth1, const CudaComplex &oth2);					// DataType + cmplx
	 friend CudaComplex operator + (const CudaComplex &oth2, const DataType oth1);					// DataType + cmplx

	 CudaComplex &operator+=(const CudaComplex &other);												// append operator
	 CudaComplex &operator+=(const DataType other);													// append with DataType

// -------- subtraction -----------
	 void sub(const CudaComplex &oth1, const CudaComplex &oth2);										// subs 2 cmplx
	 friend CudaComplex operator - (const CudaComplex &oth1, const CudaComplex &oth2);				// sub of cmplx
	 friend CudaComplex operator - (const DataType oth1, const CudaComplex &oth2);					// DataType - cmplx
	 friend CudaComplex operator - (const CudaComplex &oth2, const DataType oth1);					// - DataType + cmplx

	 CudaComplex &operator-=(const CudaComplex &other);												// deappend operator
	 CudaComplex &operator-=(const DataType other);													// deappend with DataType

// -------- multiplication --------
	 void mul(const CudaComplex &oth1, const CudaComplex &oth2);							 // mul 2 cmplx
	 friend CudaComplex operator *(const CudaComplex &oth1, const CudaComplex &oth2);	 // mult of cmplx
	 friend CudaComplex operator * (const DataType oth1, const CudaComplex &oth2);		 // DataType * cmplx
	 friend CudaComplex operator * (const CudaComplex &oth2, const DataType oth1);		 // DataType * cmplx

	 CudaComplex &operator*=(const CudaComplex &other);									 // append operator
	 CudaComplex &operator*=(const DataType other);										 // append with DataType

// -------- division -----------
	 void div(const CudaComplex &oth1, const CudaComplex &oth2);								    	// divs 2 cmplx
	 friend CudaComplex operator / (const CudaComplex &oth1, const CudaComplex &oth2);				// div of cmplx
	 friend CudaComplex operator / (const DataType oth1, const CudaComplex &oth2);					// DataType / cmplx
	 friend CudaComplex operator / (const CudaComplex &oth2, const DataType oth1);					// DataType / cmplx

	 CudaComplex &operator/=(const CudaComplex &other);												// div append operator
	 CudaComplex &operator/=(const DataType other);													// div append with DataType

// -------- math functions -----------
	 friend DataType abs(const CudaComplex &other);							// abs
	 friend DataType arg(const CudaComplex &other);							// phase radians
	 friend DataType real(const CudaComplex &other);							// real part
	 friend DataType imag(const CudaComplex &other);							// imaginary part
	 friend CudaComplex conj(const CudaComplex &other);						// the CudaComplex conjugate
	 friend CudaComplex exp(const CudaComplex &other);						// exponential

// -------- mix functions -----------
	 void display() const { printf("(%f, %f)\n",re,im); }           // displays the number

private:
	DataType re,                               // real part
			 im;                               // imaginary part
};

#endif