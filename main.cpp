/*! \mainpage
 *
 * 	This project consists of several different Potential Field Source Surface (PFSS) implementations for computations
 * 	of the synoptic coronal magnetic field configuration of the sun for periods of low solar activity.
 */


#include <stdlib.h>
#include <stdio.h>
#include <locale>

#include "src/pfssSolution.h"
#include "src/carRotInfo.h"
#include "src/filenames.h"

#include "boost/filesystem.hpp"
#include "boost/regex.hpp"
#include "boost/format.hpp"

using namespace boost::filesystem;
using namespace boost;
using namespace std;

#ifdef CUDA
#include "src/cuda_interface.h"
#endif

#include "boost/filesystem.hpp"
#include "src/filenames.h"

using namespace boost::filesystem;

uint Magline::colorInvalid;
uint Magline::colorClosed;
uint Magline::colorPositive;
uint Magline::colorNegative;

pthread_mutex_t hcImageFITS::mutexFits;

void initStatics(const char *crlist_dir)
{
    crInfoList::initStaticMembers(crlist_dir);
    Magline::initStaticMembers();
    hcImageFITS::initStaticMembers();
    printf("--- Static variables initialized.\n");
}

void clearStatics()
{
    crInfoList::clearStaticMembers();
    printf("--- Static variables cleared.\n");
}



int main(void)
{
	printf("--- Main: Started --------------------------------------------------\n");
    stringstream fn;
    fn << "../data/crlist";
    initStatics(fn.str().data());
    printf("--- Main: Statics initiated ----------------------------------------\n");

	setenv("LC_ALL", "C", 1);
	try
	{	std::locale loc("");}
	catch(const std::exception &e)
	{	cout << "Failed\n";}

	#ifdef CUDA
		int numDev;
		cudaGetDeviceCount(&numDev);
		cudaSetDevice(numDev-1);
		cudaDeviceReset();
	#endif

    stringstream inDir, outdir;
	inDir  << "../data/PFSS/input";
	outdir << "../data/PFSS/output";

	const uint numRGrid         = 35;
	const uint mapThetaRes      = 200;
	const uint mapPhiRes        = 400;

	PFSSsolution solution;
	solution.setOutDir(outdir.str().data());
	solution.batchKielGrid(inDir.str().data(), 2.1*r_sol, numRGrid, mapThetaRes, mapPhiRes, false, 1.0/1.1);

	clearStatics();
    printf("--- Last line of code, exiting program.\n");
    fflush(stdout);
}
