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


bool setConfig(const string &filename)
{
	bool retval = true;

	if(!doesFileExist(filename))
	{
		cerr << __FILE__ << ":" << __LINE__ << ": --config option given, but file '" << filename << "'does not exist.\n";
		return false;
	}

	std::ifstream config(filename);
	string line = "";
    while(getline(config, line))
    {
    	string option 		= "";
    	string argument 	= "";
    	stringstream stream(line);
    	stream >> option >> argument;
    	if(option == "dirData")			dirData			= argument;
    	if(option == "dirConfig")		dirConfig		= argument;
    }
    config.close();

    if(dirData == "")
    {
    	cerr << __FILE__ << ":" << __LINE__ << ": Configuration file '" <<  filename << "' does not contain dirData field.\n";
    	retval = false;
    }
    else if(!directoryExists(dirData)) retval = createDir(dirData);

    if(dirConfig == "")
	{
		cerr << __FILE__ << ":" << __LINE__ << ": Configuration file '" <<  filename << "' does not contain dirConfig field.\n";
		retval = false;
	}
	else retval = directoryExists(dirConfig);

	return retval;
}

bool parse(int argc, char **argv)
{
    int c;
    int digit_optind 	= 0;

    bool batchMap		= false;
    bool loadAndMap 	= false;
    uint resMapTheta	= 100;
    uint resMapPhi		= 200;
    bool inter			= false;

    bool batchCompute	= false;
    bool compute		= false;
    uint resCompR		= 15;
    hcFloat rss			= 2.5*r_sol;
    hcFloat ellipticity	= 1.0;
    methodID method		= METH_NUMERIC;

    string filename		= "";

	while (1)
	{
		int this_option_optind 	= optind ? optind : 1;
		int option_index 		= 0;
		static struct option long_options[] = {
			{"map",  		required_argument, 	0,  0},
			{"compute", 	required_argument, 	0,  0},
			{"config",  	required_argument, 	0,  0},
			{"ell",			required_argument, 	0,  0},
			{"batchcompute",required_argument, 	0,  0},
			{"batchmap",	no_argument, 		0,  0},
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
				else if(command == "map"){			loadAndMap 		= true; filename = optarg;}
				else if(command == "compute"){		compute 		= true; filename = optarg;}
				else if(command == "batchcompute"){	batchCompute 	= true; filename = optarg;}
				else if(command == "batchmap"){		batchMap	 	= true;}
				else if(command == "resCompR"){		resCompR 		= atoi(optarg);}
				else if(command == "resMapTheta"){	resMapTheta 	= atoi(optarg);}
				else if(command == "resMapPhi"){	resMapPhi 		= atoi(optarg);}
				else if(command == "ell")
				{
					hcFloat val = strtof(optarg, NULL);
					if(val < 1.0)
					{
						printErrMess(__FILE__, __LINE__, "ellipticity must be >= 1.0, you supplied --ell " + string(optarg));
						return false;
					}
					else ellipticity = val;
				}
				else if(command == "rss")
				{
					hcFloat val = strtof(optarg, NULL);
					if(val <= 1.0)
					{
						printErrMess(__FILE__, __LINE__, "source surface radius must be > 1.0, you supplied --rss " + string(optarg));
						return false;
					}
					else rss = val * r_sol;
				}

				else
				{
					cout << "option " << long_options[option_index].name;
					if (optarg)	cout << " with arg " << optarg;
					cout << "\n";
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
		cout << "non-option ARGV-elements: ";
		while (optind < argc) cout << argv[optind++] << " ";
		cout << "\n";
	}

	if(dirData == "")	// no config supplied via cmd-line, search for default config pfss/config/config
	{
		string stdFilename = "../config/config";
		if(!setConfig(stdFilename))
		{
			printErrMess(__FILE__, __LINE__, "configuration file supplied (--config filename) and standard config '" + stdFilename + "' invalid");
			return false;
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

	PFSSsolution sol;
	if(batchCompute)	sol.batchKielGrid(filename, rss, resCompR, 100, 200, false, ellipticity);
	if(batchMap)		sol.batchMap(resMapTheta, resMapPhi, inter);
	if(loadAndMap)		sol.loadAndMapKielGrid(filename, resMapTheta, resMapPhi, false);
	if(compute)			sol.computeKielGrid(filename, "", rss, resCompR, ellipticity);

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
 *  \brief This PFSS computation suite can compute the PFSS solution using photospheric magnetograms from MDI, HMI, GONG, and WSO. It can also create
 *  image-like maps of the magnetic configuration and the expansion factor at arbitrary heights between photsphere and source surface.
 *
 *  Commands are given via command line arguments to the PFSS suite, e.g. pfss \-\-compute ./synop2066.fits for computing the PFSS model with standard
 *  parameters for the Carrington rotation given in synop2066.fits. Several solar observatories are supported and automatically recognized.
 *
 *  Commands and their options are:
 *
 *  \-\-config filename
 *
 *  filename is path to configuration file, [default: \-\-config ../config/config]
 *
 *  \-\-compute filename
 *
 *	Invokes the PFSS solver for given photospheric magnetogram found at path filename (for supported observatories see above).
 *
 *	Support arguments:
 *
 *	\-\-rss value			source surface height (in solar radii), value is floating point\n
 *	\-\-ell value			ellipticity of source surface, value is floating point, [default: 1.0 (spheric)]\n
 *	\-\-resCompR value	computational grid resolution in radial direction, other directions are determined automatically,
 *							value is unsigned integer
 *
 *	\-\-map filename
 *
 *	Computes the magnetic configuration at specified height (default: photosphere and source surface)
 *
 *	Support arguments:
 *
 *	\-\-resMapTheta value		resolution of mapping in meridional direction, value is unsigned integer\n
 *	\-\-resMapPhi value		resolution of mapping in zonal direction, value is unsigned integer\n
 *
 *	\-\-batchcompute directory
 *
 *	Invokes the PFSS solver for all magnetic magnetograms found in directory (non-recursive)
 *
 *	Support arguments: see \-\-compute filename
 *
 *	\-\-batchmap
 *
 *	Invokes mapper for all solutions computed in data directory.
 */
