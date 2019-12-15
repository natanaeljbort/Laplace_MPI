Students: 
Natanael Bort Soldevila
Gael Ruta Gatera

# Laplace Parallelization (MPI)

## __Introduction__
&nbsp;&nbsp;&nbsp;&nbsp;This report was written in order to demonstrate the results of using Message Passing Interface (MPI) Parallelization on a 2 dimensional laplace equation using the jacobi iteration method. MPI defines the syntax and semantics of a core of library routines such as MPI_send & MPI_recv useful to a wide range of users writing message passing programs in C, C++ and Fortan. This experiment required students to use processes, defined as subsections of a predefined matrix , to interchanging and synchronize necessary data with other processors. 
&nbsp;&nbsp;&nbsp;&nbsp;Let’s take a scenario where we have a current processor “me”. It will receive the last row of the processor that preceded “me-1. In return, processor “me” would send it’s first row (ri) back to processor “me-1”. In addition, processor “me” would send it’s last row (rf) to processor “me+1” and this messages passing process would continue until al desired messages are successfully passed. Note that, all processors except the first processor will send their “first” row to their neighbor above “me – 1” and all processors except the last processor will send their “last” row to their neighbors below “me+1”.

## __Methods Adopted__
&nbsp;&nbsp;&nbsp;&nbsp;Using handles defined by the MPI protocol the the following calls were made in order to successfully parallelize the Laplace Equation Numerically Solved by Jacobi Iteration.

```{c}

```
