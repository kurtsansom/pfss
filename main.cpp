#include <stdlib.h>
#include <stdio.h>
#include <locale>
#include <getopt.h>

#include "src/pfssSolution.h"
#include "src/carRotInfo.h"
#include "src/filenames.h"

#include "boost/filesystem.hpp"
#include "boost/regex.hpp"
#include "boost/format.hpp"
#include "boost/algorithm/string.hpp"

using namespace boost::filesystem;
using namespace boost::algorithm;
using namespace boost;
using namespace std;

#ifdef CUDA
#include "src/cuda_interface.h"
#endif

// global variables
uint Magline::colorInvalid;
uint Magline::colorClosed;
uint Magline::colorPositive;
uint Magline::colorNegative;

pthread_mutex_t hcImageFITS::mutexFits;

string dirData 		= "";
string dirConfig 	= "";

bool initStatics()
{
	bool retval = true;
    retval &= crInfoList::initStaticMembers();
    Magline::initStaticMembers();
    hcImageFITS::initStaticMembers();
    return retval;
}

void clearStatics()
{
    crInfoList::clearStaticMembers();
}

void printHelp()
{
	printStdOutMess(__FILE__, __LINE__, "this text shall be displayed if --help is set or no command line arguments are given");
}

bool parse(int argc, char **argv)
{
    int c;
    int digit_optind 	= 0;

    bool abort			= false;
    bool help			= false;
    bool batchMap		= false;
    bool loadAndMap 	= false;
    uint resMapTheta	= DEFAULT_RES_MAP_THETA;
    uint resMapPhi		= DEFAULT_RES_MAP_PHI;

    bool batchCompute	= false;
    bool compute		= false;
    uint resCompR		= DEFAULT_RES_COMP_R;
    hcFloat rss			= DEFAULT_RSS * r_sol;
    hcFloat height		= 0.0;
    hcFloat ellipticity	= DEFAULT_ELL;
    uint order			= DEFAULT_SHC_ORDER;
    methodID method		= METH_NUMERIC;

    string filename		= "";

	while (1)
	{
		int this_option_optind 	= optind ? optind : 1;
		int option_index 		= 0;
		static struct option long_options[] = {
			{"help",  		no_argument,	 	0,  0},
			{"map",  		required_argument, 	0,  0},
			{"height",  	required_argument, 	0,  0},
			{"compute", 	required_argument, 	0,  0},
			{"config",  	required_argument, 	0,  0},
			{"ell",			required_argument, 	0,  0},
			{"method",		required_argument, 	0,  0},
			{"order",		required_argument, 	0,  0},
			{"batchcompute",required_argument, 	0,  0},
			{"batchmap",	no_argument, 		0,  0},
			{"help",		no_argument, 		0,  0},
			{"ell",			required_argument, 	0,  0},
			{"rss",    		required_argument, 	0,  0},
			{"resCompR",	required_argument, 	0,  0},
			{"resMapTheta",	required_argument, 	0,  0},
			{"resMapPhi",	required_argument, 	0,  0},
			{0,         	0,                 	0,  0}
		};

		c = getopt_long(argc, argv, "abc:d:012", long_options, &option_index);
		if (c == -1)	break;

		string command = long_options[option_index].name;
		to_lower(command);
		switch (c)
		{
			case 0:
				if(		command == "config")		setConfig(optarg);
				else if(command == "help")			help 			= true;
				else if(command == "map"){			loadAndMap 		= true; filename = optarg;}
				else if(command == "compute"){		compute 		= true; filename = optarg;}
				else if(command == "batchcompute"){	batchCompute 	= true; filename = optarg;}
				else if(command == "batchmap"){		batchMap	 	= true;}
				else if(command == "help"){			help		 	= true;}
				else if(command == "rescompr"){		resCompR 		= atoi(optarg);}
				else if(command == "resmaptheta"){	resMapTheta 	= atoi(optarg);}
				else if(command == "resmapphi"){	resMapPhi 		= atoi(optarg);}
				else if(command == "order"){		order 			= atoi(optarg);}
				else if(command == "ell")
				{
					hcFloat val = strtof(optarg, NULL);
					if(val < 0.0)
					{
						printErrMess(__FILE__, __LINE__, "ellipticity must be >= 0.0, you supplied --ell " + string(optarg));
						abort = true;
						//return false;
					}
					else ellipticity = val;
				}
				else if(command == "method")
				{
					if(strcmp(optarg, "shc") && strcmp(optarg, "numeric"))
					{
						printErrMess(__FILE__, __LINE__, "you supplied --method " + string(optarg) + ", but only 'numeric' and 'shc' are supported");
						abort = true;
						//return false;
					}

					if(!strcmp(optarg, "shc"))			method = METH_SHC;
					else if(!strcmp(optarg, "numeric"))	method = METH_NUMERIC;
				}

				else if(command == "rss")
				{
					cout << "set fucking rss\n";fflush(stdout);
					hcFloat val = strtof(optarg, NULL);
					if(val <= 1.0)
					{
						printErrMess(__FILE__, __LINE__, "source surface radius must be > 1.0, you supplied --rss " + string(optarg));
						abort = true;
						//return false;
					}
					else rss = val * r_sol;

					cout << "set rss: " << rss << "\n";
				}
				else if(command == "height")
				{
					hcFloat val = strtof(optarg, NULL);
					if(val < 1.0 || val > rss / r_sol)
					{
						printErrMess(__FILE__, __LINE__, "you supplied --height " + string(optarg) + ", but it must hold 1.0 <= height <= r_ss");
						abort = true;
						//return false;
					}
					else height = val * r_sol;
				}
				else
				{
					printStdOutMess(__FILE__, __LINE__, "unknown option " + string(long_options[option_index].name) + " with arg " + (optarg ? string(optarg) : "none"));
					abort = true;
					//return false;
				}
				break;

			case '0':
			case '1':
			case '2':
				if (digit_optind != 0 && digit_optind != this_option_optind)
					cout << "digits occur in two different argv-elements.\n";
				digit_optind = this_option_optind;
				printf("option %c\n", c);
				break;

			case 'a':
				printf("option a\n");
				break;

			case 'b':
				printf("option b\n");
				break;

			case 'c':
				printf("option c with value '%s'\n", optarg);
				break;

			case 'd':
				printf("option d with value '%s'\n", optarg);
				break;

			case '?':
				break;

			default:
				printf("?? getopt returned character code 0%o ??\n", c);
		}
	}

	if (optind < argc)
	{
		printErrMess(__FILE__, __LINE__, "non-option ARGV-elements: ");
		while (optind < argc) printErrMess(__FILE__, __LINE__,  string(argv[optind++]) + " ");
		abort = true;
		//return false;
	}

	if(method != METH_NUMERIC && ellipticity != 1.0)
	{
		printErrMess(__FILE__, __LINE__, "--ell " + toStr(ellipticity) + "not supported with computation method " + getStringFromMethodID(method));
		abort = true;
	}

	if(dirData == "")	// no config supplied via cmd-line, search for default config pfss/config/config
	{
		string stdFilename = "../config/config";
		if(!setConfig(stdFilename))
		{
			printErrMess(__FILE__, __LINE__, "configuration file supplied (--config filename) and standard config '" + stdFilename + "' invalid");
			abort = true;
			//return false;
		}
		if(!initStatics()) return false;
		setenv("LC_ALL", "C", 1);
		try
		{	std::locale loc("");}
		catch(const std::exception &e)
		{	cout << "Failed\n";return false;}

		#ifdef CUDA
			int numDev;
			cudaGetDeviceCount(&numDev);
			cudaSetDevice(numDev-1);
			cudaDeviceReset();
		#endif
	}

	if(abort)
	{
		printErrMess(__FILE__, __LINE__, "commmand line option/configuration file error");
		return false;
	}

	cout << "rss: " << rss << "\n";
	// perform operations requested by command line
	if(help || argc==1)	printHelp();
	PFSSsolution sol;
	if(batchCompute && method==METH_NUMERIC)	sol.batchKielGrid(filename, rss, resCompR, ellipticity);
	if(batchCompute && method==METH_SHC)		sol.batchKielSHC(filename, rss, resCompR, order);
	if(batchMap)								sol.batchMap(resMapTheta, resMapPhi, height);
	if(loadAndMap)								sol.loadAndMapKielGrid(filename, resMapTheta, resMapPhi, height);
	if(compute && method==METH_NUMERIC)			sol.computeKielGrid(filename, "", rss, resCompR, ellipticity);
	if(compute && method==METH_SHC)				sol.computeKielSHC(filename, order, rss, resCompR);

	return true;
}

int main(int argc, char **argv)
{
	hcDate now;
	stringstream ss;
	ss << right << setw(2) << setfill(' ') << to_string(PFSS_VER_MAJ) << "." << to_string(PFSS_VER_MIN);
	cout << "---------------------------------------------------------------------------------------\n";fflush(stdout);
	cout << "--- PFSS Computation Suite (ver. " << ss.str() << ") ------------------------ " << now.toSpiceString() << " ---\n";fflush(stdout);
	cout << "---------------------------------------------------------------------------------------\n";fflush(stdout);
	parse(argc, argv);
	clearStatics();
	now.setFromSystemTime();
	cout << "---------------------------------------------------------------------------------------\n";fflush(stdout);
	cout << "--- PFSS Computation Suite: Exiting program. ------------------ " << now.toSpiceString() << " ---\n";fflush(stdout);
	cout << "---------------------------------------------------------------------------------------\n";fflush(stdout);
}

/*! \mainpage
 *
 *  \brief This PFSS computation suite computes the PFSS solution using photospheric magnetograms from MDI, HMI, GONG, and WSO. It can also create
 *  maps of the magnetic configuration and the expansion factor at arbitrary heights between photsphere and source surface.
 *  Commands are given via command line arguments to the PFSS suite. For an example how to use the command-line interface see below.
 *  Several solar observatories are supported and automatically recognized for supplying the photospheric magnetogram.
 *
 *  \tableofcontents
 *
 *  # Building the documentation (optional)
 *
 *  A PDF containing the documentation of this program is included in the GIT repository. The documentation is implemented via doxygen. In order to obtain the most current documentation, you first need to install doxygen
 *
 *		apt install doxygen
 *
 *  and then build the documentation:
 *
 *		cd PFSS
 *		make documentation
 *
 *  The documentation will then be produced both as HTML and PDF in the directory PFSS/doc.
 *
 *  # Optimizations (optional)
 *
 *  Several optimizations can be adjusted in the Makefile. If a CUDA-capable device is present, setting the variable CUDA to '1' will build the program to employ the CUDA device, which might decrease computation time substantially.
 *  For this optimization the CUDA runtime environment needs to be installed (please consult NVIDIA's webpage for instructions).
 *
 *  The variable NUMTHREADS limits the number of threads to be utilized for multithreaded execution of the program. For optimal performance it should be the same number as CPU cores in the system.
 *
 *  Default parameters for the computation can also be set in the Makefile. The corresponding build variables start with DEFAULT_.
 *
 *  # Building the binary
 *
 *  This program utilizes several third-party libraries which need to be installed in order for the binary to be built. These libraries are very common so there is a good chance that they can be found in your distributions repositories.
 * 	If you are on a debian system, try installing them via
 *
 * 		apt install libfreeimage-dev libcfitsio-dev libboost-dev libboost-filesystem-dev libboost-regex-dev
 *
 * 	If successful the PFSS computation suite can be built via
 *
 *		cd PFSS
 *		make
 *
 *	The binary will be placed in PFSS/bin.
 *
 *  # Running the program
 *
 *  If you used make to build the PFSS computation suite, the binary file will be stored in PFSS/bin/.
 *
 *  ## Configuration and model output files
 *
 *  Upon execution the binary reads a configuration file, which specifies the data and configuration directories.
 *  The configuration directory contains information about start and stop times of Carrington rotations. The data directory contains all the output from the PFSS computation suite. If you run the binary from the PFSS/bin/ directory without
 *  specifying a configuration file, the default file PFSS/config/config will be used. The default data directoy is then PFSS/data. Please consult this default config file to set your own data directory at a location with enough disk space
 *  if the default location is not suitable. Absolute paths in your config file allows the binary to be executed from arbitrary shell locations. Manipulation of the configuration file is only necessary if you have special needs for the location of
 *  the computed solutions.
 *
 *  ## Command line options and arguments
 *
 *  The syntax to run the PFSS computation suite is
 *
 *  	pfss --option0 [argument0] --option1 [argument1] ...
 *
 *  The following options and arguments are supported. Default values for not specified options can be adjusted in the Makefile.
 *
 *   \-\- **config** _filename_
 *
 *  _filename_ is path to configuration file [default: \-\-config ../config/config]
 *
 *  \-\- **compute** _filename_
 *
 *	Invokes the PFSS solver for given photospheric magnetogram at path _filename_. For additional arguments see below.
 *
 *	\-\- **map** _filename_
 *
 *	Computes the magnetic configuration at specified height [default: photosphere and source surface]. _filename_ points to the configuration file of the computed solution (ending in \*config.cfg). For additional options see below.
 *
 *	\-\- **batchcompute** _directory_
 *
 *	Invokes the PFSS solver for all magnetic magnetograms found in _directory_ (non-recursive).
 *
 *	\-\- **batchmap**
 *
 *	Invokes mapper for all solutions computed in data directory.
 *
 *	**Additional options for compute**
 *
 *	\-\- **rss** _value_		source surface height (multiples of solar radius), _value_ is floating point\n
 *	\-\- **ell** _value_		ellipticity of source surface, _value_ is floating point, [default: 1.0 (spheric)]\n
 *	\-\- **resCompR** _value_	computational grid resolution in radial direction, other directions are determined automatically,
 *								_value_ is unsigned integer\n
 *	\-\- **method** _value_		solution method to be used for PFSS computation. _value_ is either 'shc' for the classic spherical harmonic coefficient approach or 'numeric' for the finite difference solver\n
 *	\-\- **order** _value_		maximum principal order to be used with the SHC approach, _value_ is unsigned integer
 *
 *	**Additional options for map**
 *
 *	\-\- **resMapTheta** _value_		resolution of mapping in meridional direction, _value_ is unsigned integer\n
 *	\-\- **resMapPhi** _value_			resolution of mapping in zonal direction, _value_ is unsigned integer\n
 *	\-\- **height** _value_			height between phot. and source surface to be mapped, _value_ is multiple of solar radius\n
 *
 *	# Supported photospheric magnetograms
 *
 *	The following photospheric synoptic magnetogram sources are supported. The PFSS computation suite identifies the source instrument and necessary pre-processing steps by filename, meaning you
 *	cannot alter them or the program will not be able to handle the containing data.
 *
 *	| observatory name				| filename example 				| resolution	| URL 										|
 *	|----------------------------	|------------------				|------------	|-----										|
 *	| WSO							| WSO.2066.F.txt				| 72 x 30		| http://wso.stanford.edu/synopticl.html	|
 *	| SOHO-MDI						| synop_Ml_0.2066.fits			| 3600 x 1080	| http://sun.stanford.edu/synop/ 			|
 *	| SDO-HMI						| hmi.Synoptic_Ml.2100.fits		| 3600 x 1440	| http://hmi.stanford.edu/data/synoptic.html|
 *	| NSO-GONG						| mrmqs080208t0128c2066_000.fits| 360 x 180		| https://gong.nso.edu/data/magmap/crmap.html|
 *
 *	# Example use cases
 *
 *	All examples are executed from the PFSS/bin-directory, so change there:
 *
 *		cd PFSS/bin
 *
 *  To perform a PFSS model evaluation with default parameters for the synoptic photospheric magnetogram found at PFSS/data/input/synop_Ml_0.2066.fits (not included in the repository, you need to download the magnetogram and place it at that location):
 *
 *  	./pfss --compute ../data/input/synop_Ml_0.2066.fits
 *
 *  To generate magnetic mappings at the photosphere and source surface of the computed solution in the previous example:
 *
 *		./pfss --map ../data/2066/2066_MDI_Kiel_PFSS2.50_GridSph_35x87x175_config.cfg
 *
 *	To generate a magnetic mappings of the same solution at height r=2.4 r_sol with a resolution of 130 x 200 pixels:
 *
 *		./pfss --map ../data/2066/2066_MDI_Kiel_PFSS2.50_GridSph_35x87x175_config.cfg --height 2.4 \
 *		--resMapTheta 130 --resMapPhi 200
 *
 *	To evaluate the PFSS model for all photospheric magnetic maps found in directory PFSS/data/batchinput with a radial grid resolution of 40 grid points and a source surface radius of r=3.1 r_sol:
 *
 *		./pfss --batchcompute ../data/input/ --resCompR 40 --rss 3.1
 */
