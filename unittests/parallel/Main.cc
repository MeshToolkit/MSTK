#include <UnitTest++.h>
#include <TestReporterStdout.h>

#ifdef MSTK_HAVE_MPI
#include <mpi.h>
#endif

int main(int argc, char *argv[])
{

  return UnitTest::RunAllTests ();

}

