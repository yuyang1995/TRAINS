#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <unistd.h>

#include "trains.h"

/*! \file main.c
 *  \brief start of the program
 */

/*!
 *  This function initializes the MPI communication packages.
 */

int main(int argc, char **argv)
{
  double t0 = 0, t1, dt;
  int ht, mt;
  double st;

  /* initialize MPI */
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &thistask);
  MPI_Comm_size(MPI_COMM_WORLD, &totaltask);
  MPI_Get_processor_name(proc_name, &namelen);

  if (thistask == roottask)
  {
    t0 = second();
    printf("===============TRAINS==================\n");
    printf("Starts to run...\n");
    printf("%d cores used.\n", totaltask);
  }

  if (command_line_options(argc, argv) != EXIT_SUCCESS)
  {
    MPI_Finalize();
    return 0;
  }

  if (parset.flag_help == 0)
  {
    begin_run(); /* implementation run */

    end_run(); /* end run */
  }

  MPI_Finalize(); /* clean up and finalize MPI */
  if (thistask == roottask)
  {
    t1 = second();
    dt = timediff(t0, t1);
    get_hms(dt, &ht, &mt, &st);
    printf("Time used: %dh %dm %fs.\n", ht, mt, st);
    printf("Ends successfully.\n");
    printf("===============TRAINS==================\n");
  }
  return 0;
}