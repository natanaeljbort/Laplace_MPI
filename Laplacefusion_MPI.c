#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

float stencil ( float v1, float v2, float v3, float v4)
{
  return (v1 + v2 + v3 + v4) * 0.25f;
}

float max_error ( float prev_error, float old, float new )
{
  float t= fabsf( new - old );
  return t>prev_error? t: prev_error;
}

float laplace_step(float *in, float *out, int n, int me, int nproces)
{
  int i, j, intj, finj;
  float error=0.0f;
  if(me==0){intj=1;}
  else{intj=0;}
  if(me==nproces-1){finj=(n/nproces)-1;}
  else{finj=n/nproces;}
  for ( j=intj; j < finj; j++ )
    for ( i=1; i < n-1; i++ )
    {if(j==0){
       out[j*n+i]= stencil(in[j*n+i+1], in[j*n+i-1], in[(n/nproces)*n+i], in[(j+1)*n+i]);
       }
     else if(j==(n/nproces)-1){
       out[j*n+i]= stencil(in[j*n+i+1], in[j*n+i-1], in[(j-1)*n+i], in[((n/nproces)+1)*n+i]);
       }
    else{
      out[j*n+i]= stencil(in[j*n+i+1], in[j*n+i-1], in[(j-1)*n+i], in[(j+1)*n+i]);
       }
      error = max_error( error, out[j*n+i], in[j*n+i] );
    }
  return error;
}


void laplace_init(float *in, int n, int me, int nproces)
{
  int i;
  const float pi  = 2.0f * asinf(1.0f);
  memset(in, 0, (n*(n/nproces)+2)*sizeof(float));
  for (i=0; i<n/nproces; i++) {
    float V = in[i*n] = sinf(pi*(i+me*(n/nproces)) / (n-1));
    in[ i*n+n-1 ] = V*expf(-pi);
  }
}

int main(int argc, char** argv)
{
  int n = 400;
  int iter_max = 1000;
  int ri, rf;
  float *A, *temp;
  int nproces, me;
  double MPI_Wtime(), inicio, final, tiempo_total;

  inicio=MPI_Wtime();


  // Initializing the MPI environment
  MPI_Init(&argc, &argv);
  // Will use the "me" variable in order to determine the current process
  // "ndims" variable for n number of dimensions
  MPI_Comm_rank(MPI_COMM_WORLD, &me);

  // "nproc" variable for n number of processes
  MPI_Comm_size(MPI_COMM_WORLD, &nproces);


  const float tol = 1.0e-5f;
  float error= 1.0f;
  // get runtime arguments
  if (argc>1) {  n        = atoi(argv[1]); }
  if (argc>2) {  iter_max = atoi(argv[2]); }

  A    = (float*) malloc( (n*(n/nproces)+2)*sizeof(float) );
  temp = (float*) malloc( (n*(n/nproces)+2)*sizeof(float) );
 //  set boundary conditions
  laplace_init (A, n, me, nproces);
  laplace_init (temp, n, me, nproces);
  //A[(n/128)*n+n/128] = 1.0f; // set singular point

  printf("Jacobi relaxation Calculation: %d x %d mesh,"
         " maximum of %d iterations\n",
         n, n, iter_max );

  int iter = 0;
  MPI_Status status;


  while ( error > tol*tol && iter < iter_max )
  { //printf("1 %d %d\n",iter, me);
    iter++;
    if(me>0){
        MPI_Send(&A[0], n, MPI_FLOAT, me-1, 1, MPI_COMM_WORLD);
        MPI_Recv(&A[n*(n/nproces)],n, MPI_FLOAT, me-1,1,MPI_COMM_WORLD, &status);}
    if(me<nproces-1){
        MPI_Send(&A[n*((n/nproces)-1)],n,MPI_FLOAT,me+1,1,MPI_COMM_WORLD);
        MPI_Recv(&A[n*((n/nproces)+1)],n,MPI_FLOAT,me+1,1,MPI_COMM_WORLD, &status);
    }
//printf("2 %d %d\n",iter, me);
    error= laplace_step (A, temp, n, me, nproces);
    float *swap= A; A=temp; temp= swap; // swap pointers A & temp
    //printf("3 %d %d\n",iter, me);
  }
  error = sqrtf( error );
  final=MPI_Wtime();
  tiempo_total=final-inicio;
  printf("Total Iterations: %5d, ERROR: %0.6f, ", iter, error);
  printf("A[%d][%d]= %0.6f\n", n/100, n/100, A[(n/100)*n+n/100]);
  printf("The total time of the program in rank %d is %f \n", me,tiempo_total);

  free(A); free(temp);

}
