#include <UnitTest++.h>
#include <TestReporterStdout.h>

#include <mpi.h>

int main(int argc, char *argv[])
{
  int status;

  MPI_Init(&argc, &argv);

  status = UnitTest::RunAllTests ();

  MPI_Finalize();

  return status;
}

