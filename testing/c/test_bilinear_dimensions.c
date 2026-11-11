/*
 *	This file is part of qpOASES.
 *
 *	qpOASES -- An Implementation of the Online Active Set Strategy.
 *	Copyright (C) 2007-2015 by Hans Joachim Ferreau, Andreas Potschka,
 *	Christian Kirches et al. All rights reserved.
 *
 *	qpOASES is free software; you can redistribute it and/or
 *	modify it under the terms of the GNU Lesser General Public
 *	License as published by the Free Software Foundation; either
 *	version 2.1 of the License, or (at your option) any later version.
 *
 *	qpOASES is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *	See the GNU Lesser General Public License for more details.
 *
 *	You should have received a copy of the GNU Lesser General Public
 *	License along with qpOASES; if not, write to the Free Software
 *	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */


/**
 *	\file testing/c/test_bilinear_dimensions.c
 *	\author Prajval Arora
 *	\version 3.1embedded
 *	\date 2026
 *
 *	Regression test for https://github.com/acados/qpOASES/issues/8:
 *	DenseMatrix_bilinear uses a fixed Ax[NVCMAX*NVMAX] buffer but writes
 *	icols->length * xN entries into it. QPs with more than QPOASES_NVMAX
 *	free variables (e.g. full condensing with slack variables) must be
 *	rejected with RET_INVALID_ARGUMENTS instead of overflowing the buffer.
 */


#include <stdlib.h>
#include <stdio.h>
#include <qpOASES_e.h>


USING_NAMESPACE_QPOASES

int main( )
{
	/* dimensions from the reported acados FULL_CONDENSING_QPOASES segfault:
	 * 703 free indices, 700 free variables -> 492100 buffer writes,
	 * while the fixed Ax buffer holds NVCMAX*NVMAX = 709*287 = 203483 entries */
	const int nFR = 703;
	const int nZ  = 700;

	real_t* Hval = (real_t*)malloc( (size_t)nFR*nFR*sizeof(real_t) );
	real_t* x    = (real_t*)malloc( (size_t)nFR*nZ*sizeof(real_t) );
	real_t* y    = (real_t*)malloc( (size_t)nZ*nZ*sizeof(real_t) );
	int i;
	returnValue ret;

	for ( i = 0; i < nFR*nFR; ++i ) Hval[i] = 0.0;
	for ( i = 0; i < nFR; ++i ) Hval[i*nFR+i] = 1.0;
	for ( i = 0; i < nFR*nZ; ++i ) x[i] = 1.0;

	DenseMatrix H;
	H.val = Hval;
	H.nRows = nFR;
	H.nCols = nFR;
	H.leaDim = nFR;

	Indexlist* icols = Indexlist_createMemory( nFR );
	IndexlistCON( icols,nFR );
	Indexlist_init( icols,nFR );
	for ( i = 0; i < nFR; ++i )
		Indexlist_addNumber( icols,i );

	ret = DenseMatrix_bilinear( &H,icols,nZ, x,nFR, y,nZ );

	free( Hval );
	free( x );
	free( y );

	if ( ret == RET_INVALID_ARGUMENTS )
	{
		fprintf( stdFile,"test_bilinear_dimensions: PASS (oversized bilinear rejected with RET_INVALID_ARGUMENTS)\n" );
		return 0;
	}

	fprintf( stdFile,"test_bilinear_dimensions: FAIL (oversized bilinear not rejected, fixed Ax buffer was overrun)\n" );
	return 1;
}


/*
 *	end of file
 */
