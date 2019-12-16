# Laplace Parallelization (MPI)

## Introduction

&nbsp;&nbsp;&nbsp;&nbsp;This report was written in order to demonstrate the results of using Message Passing Interface (MPI) Parallelization on a 2 dimensional laplace equation using the jacobi iteration method. MPI defines the syntax and semantics of a core of library routines such as MPI_send & MPI_recv useful to a wide range of users writing message passing programs in C, C++ and Fortan. This experiment required students to use processes, defined as subsections of a predefined matrix , to interchanging and synchronize necessary data with other processors. 
&nbsp;&nbsp;&nbsp;&nbsp;Let’s take a scenario where we have a current processor “me”. It will receive the last row of the processor that preceded “me-1. In return, processor “me” would send it’s first row (ri) back to processor “me-1”. In addition, processor “me” would send it’s last row (rf) to processor “me+1” and this messages passing process would continue until al desired messages are successfully passed. Note that, all processors except the first processor will send their “first” row to their neighbor above “me – 1” and all processors except the last processor will send their “last” row to their neighbors below “me+1”.

## Methods Adopted
&nbsp;&nbsp;&nbsp;&nbsp;Using handles defined by the MPI protocol the the following calls were made in order to successfully parallelize the Laplace Equation Numerically Solved by Jacobi Iteration.

* Processor number
```c
MPI_Comm_rank(MPI_COMM_WORLD, &me);
```
* Number of Processors
```c
MPI_Comm_size(MPI_COMM_WORLD, &nproces);
```
* Start the MPI Timer
```c
inicio=MPI_Wtime();
```
* Stop the MPI Timer
```c
final=MPI_Wtime();
```
* Calculating the change in time
```c
tiempo_total=final-inicio;
```
* Allocating Memory for Matrix A and temperature
```c
A    = (float*) malloc( (n*(n/nproces)+2)*sizeof(float) );
temp = (float*) malloc( (n*(n/nproces)+2)*sizeof(float) );
```
* Including variables “me” and “nproces” when initializing boundary conditions.
** For Matrix A
```c
laplace_init (A, n, me, nproces);
For temperature “temp”
laplace_init (temp, n, me, nproces);
```
Conditional statements that were used in order to send the default communicator called “MPI_COMM_WORLD”.
* Sending communicator with message for “me” processor to processor “me-1”. Starting at pointer &A[0]for n columns
```c
if(me>0){
   MPI_Send(&A[0], n, MPI_FLOAT, me-1, 1, MPI_COMM_WORLD);
   MPI_Recv(&A[n*(n/nproces)],n, MPI_FLOAT, me-1,1,MPI_COMM_WORLD, &status);}
```
* If processor “me” is less than “nprocess-1” then the message is passed down to “me+1". starting at pointer &A[n*((n/nproces)-1)] for n columns using datatype MPI_FLOAT.
```c
if(me<nproces-1){
  MPI_Send(&A[n*((n/nproces)-1)],n,MPI_FLOAT,me+1,1,MPI_COMM_WORLD);
  MPI_Recv(&A[n*((n/nproces)+1)],n,MPI_FLOAT,me+1,1,MPI_COMM_WORLD, &status);}
```
* Interchanging and determining global error global error
```c
error= laplace_step (A, temp, n, me, nproces);
float *swap= A; A=temp; temp= swap;
```
