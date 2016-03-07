#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../LIB/zheads.h"
#include "../LIB/zdefs.h" 
#include "../LIB/zprotos.h"  
#include "../ios.h" 
#include <sys/time.h>
#include <string.h>
#include <complex.h>

/*----------------------------------------------------------------------
| main test driver for ZILUK --> Complex version of ILUK Preconditioner |
| Report bugs / send comments to: saad@cs.umn.edu                       |
+-----------------------------------------------------------------------*/

/*-------------------- protos */
void output_header( io_t *pio );
void output_result(int lfil, io_t *pio, int iparam );
int zread_coo(complex double **VAL, int **COL, int **ROW, io_t *pio,
	      complex double **rhs, complex double **sol);
int zread_inputs( char *in_file, io_t *pio );
int zget_matrix_info(FILE *fmat, io_t *pio );
void zrandvec (complex double *v, int n);
/*-------------------- end protos */
int main(){
  int ierr = 0;
/*-------------------------------------------------------------------
 * options
 *-----------------------------------------------------------------*/
  int plotting = 0, skip_its = 0;
  char pltfile[256];

  FILE *fits = NULL;
  csptr csmat = NULL;
  SMatptr MAT;
  SPreptr PRE; 
  complex double *VAL;
  int *COL, *ROW;  
  iluptr lu = NULL;
  complex double *sol = NULL, *x = NULL, *rhs = NULL;
  int n, nnz, lfil; 
  
  FILE *flog = stdout, *fmat = NULL;
  io_t io;
  double tm1, tm2;
  
  int mat, numat, iparam, i;
  double terr;
  char line[MAX_LINE];
  MAT = (SMatptr)Malloc( sizeof(SMat), "main:MAT" );
  PRE = (SPreptr)Malloc( sizeof(SPre), "main:PRE" );
/*-------------------------------------------------------------------
 * reads matrix 
 * solves using block level of fill ILU preconditioned fgmres
 *-----------------------------------------------------------------*/
  memset( &io, 0, sizeof(io) );
/*-----------------------------------------------------------------*/
  if( zread_inputs( "inputs", &io ) != 0 ) {
    fprintf( flog, "Invalid inputs file...\n" );
    goto ERROR_HANDLE;
  }
/*-----------------------------------------------------------------*/
  
  if( NULL == ( fmat = fopen( "matfile_coo", "r" ) ) ) {
    fprintf( flog, "Can't open matfile_coo...\n" );
    goto ERROR_HANDLE;
  }
  memset( line, 0, MAX_LINE );
  fgets( line, MAX_LINE, fmat );
  if( ( numat = atoi( line ) ) <= 0 ) {
    fprintf( flog, "Invalid count of matrices...\n" );
    goto ERROR_HANDLE;
  }

/*-------------------- open file OUT/ILUK.out for all performance
                       results of this run (all matrices and params) 
                       also set io->PrecMeth */
    /* sprintf( io.outfile, "OUT/%s_ILUK.out", io.HBnameF );*/
    strcpy(io.outfile,"OUT/ILUK.out");
    strcpy(io.PrecMeth,"ILUK");
    if( NULL == ( io.fout = fopen( io.outfile, "w" ) ) ) {
      fprintf(flog,"Can't open output file %s...\n", io.outfile);
      goto ERROR_HANDLE;
   }
/*-------------------- LOOP THROUGH MATRICES-*/
  for( mat = 1; mat <= numat; mat++ ) {
    if( zget_matrix_info( fmat, &io ) != 0 ) {
      fprintf( flog, "Invalid format in matfile...\n" );
      goto ERROR_HANDLE;
    }
    fprintf( flog, "MATRIX: %s...\n", io.HBnameF );
/*-------------------- Read in matrix and allocate memory-------*/
    csmat = (csptr)Malloc( sizeof(zSparMat), "main" );
    zread_coo(&VAL, &COL, &ROW, &io, &rhs, &sol) ;
    n = io.ndim;
    nnz = io.nnz;
    fprintf(stdout,"  matrix read successfully \n");
/*-------------------- conversion from COO to CS format */
    if( ( ierr = zCOOcs( n, nnz, VAL, COL, ROW, csmat ) ) != 0 ) {
      fprintf( stderr, "ILUK: zCOOcs error\n" );
      return ierr;
    }
/*-------------------- COO arrays no longer needed -- free */    
    free(ROW); ROW = NULL;
    free(VAL); VAL = NULL;
    free(COL); COL = NULL;
/*----------------------------------------------------------------------
|  The right-hand side is generated by assuming the solution is
|  a vector of ones.  To be changed if rhs is available from data.
|---------------------------------------------------------------------*/
    x=(complex double *)Malloc( n * sizeof(complex double), "main" );
    for( i = 0; i < n; i++ ) 
      x[i] = 1.0 + 0.0*I; 
    zmatvec(csmat, x, rhs);
    output_header( &io );
    lfil = io.fill_lev;
    /*  make sure this is set to zero*/
    io.tol0 = 0.0;
/*-------------------- LOOP THROUGH PARAMETERS */
    for( iparam = 1; iparam <= io.nparam; iparam++ ) {
      fprintf( flog, "iparam = %d\n", iparam );
      lu = (iluptr)Malloc( sizeof(zILUSpar), "main" );
      fprintf( flog, "begin iluk\n" );
      tm1 = sys_timer();
      ierr = zilukC( lfil, csmat, lu, flog );
      tm2 = sys_timer();
      if( ierr == -2 ) {
	fprintf( io.fout, "zero diagonal element found...\n" );
	zcleanILU( lu );
	goto NEXT_MAT;
      } else if( ierr != 0 ) {
	fprintf( flog, "*** iluk error, ierr != 0 ***\n" );
	goto ERROR_HANDLE;
      }
      io.tm_p = tm2 - tm1;
      io.fillfact = (double)znnz_ilu( lu, flog )/(double)(io.nnz+1);
      fprintf(flog,"iluk ends, fill factor (mem used) = %f\n",io.fillfact);
      
      if( skip_its ) {
	io.its = -1;
	io.tm_i = -1;
	io.enorm = -1;
	io.rnorm = -1;
	goto NEXT_PARA;
      }
      // check condition number estimation 
      if( zcondestLU( lu, sol, x, flog ) != 0 ) {
	fprintf( flog, "Not attempting iterative solution.\n" );
	fprintf( io.fout, "Not attempting iterative solution.\n" );
	io.its = -1;
	io.tm_i = -1;
	io.enorm = -1;
	io.rnorm = -1;
	goto NEXT_PARA;
      }
      zrandvec (sol,n);  /* initial guess */      
/*-------------------- create a file for printing
  'its -- time -- res' info from fgmres */
      if (plotting ) { 
	sprintf( pltfile, "OUT/%s_ILUK_F%05d", io.HBnameF, lfil);
	if( NULL == ( fits = fopen( pltfile, "w" ) ) ) {
	  fprintf( flog, "Can't open output file %s...\n", pltfile );
	  goto ERROR_HANDLE;
	}
      } else 
	fits  =NULL;
      
      io.its = io.maxits;
      tm1 = sys_timer();
/*-------------------- set up the structs before calling fgmr */
      MAT->n = n;
      MAT->CSR = csmat;
      MAT->zmatvec = zmatvecCSR; 
      PRE->ILU = lu;
      PRE->zprecon = zpreconILU;
/*-------------------- call fgmr */
      zfgmres(MAT, PRE, rhs, sol, io.tol, io.im, &io.its, fits);
      tm2 = sys_timer();
      io.tm_i = tm2 - tm1;
       if( io.its < io.maxits ) 
 fprintf(flog,"param %03d OK: converged in %d steps...\n\n",iparam,io.its );
      else 
 fprintf( flog, "not converged in %d steps...\n\n", io.maxits );
      if( fits ) fclose( fits );
/*-------------------- calculate error and residual norms */
      terr = 0.0;
      for( i = 0; i < n; i++ )
	terr += pow(cabs(sol[i] - 1.0 ),2);
      io.enorm = sqrt(terr);
      zmatvec( csmat, sol, x);
      terr = 0.0;
      for( i = 0; i < n; i++ )
	terr += pow(cabs( rhs[i] - x[i] ),2);
      io.rnorm = sqrt(terr);
/*-------------------- Test with next params   */           
NEXT_PARA:
      output_result(lfil, &io, iparam );
      lfil += io.fill_lev_inc;
      zcleanILU( lu );
    }
/*-------------------- Test with next matrix   */              
NEXT_MAT:
    zcleanCS( csmat );
    free( sol );
    free( x );
    free( rhs );
  }
  fclose( io.fout );
  if( flog != stdout ) fclose ( flog );
  fclose( fmat );
  free(MAT) ; 
  free (PRE); 
  return 0;
  
ERROR_HANDLE:
  errexit( "main.c\n" );
  return(1) ;
}
