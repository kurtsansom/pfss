#include "src/carRotInfo.h"
#include "src/filenames.h"

#include "engine/hcTools.h"

#include <boost/regex.hpp>

#include <iomanip>
#include <limits>
#include <cstring>

uint crInfoList::numCRinList            = 0;
int *crInfoList::listCRNum              = NULL;
long int *crInfoList::listJulianDayNum  = NULL;
double *crInfoList::listJulianDayFrac   = NULL;
hcDate *crInfoList::dateStart           = NULL;

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//
//                          CRinfo
//
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

crInfo::crInfo(uint crnum)
{
	initNULL();
	init(crnum);
}

crInfo::crInfo(const crInfo &other)
{
	initNULL();
	operator=(other);
}

crInfo::~crInfo(){

	clear();
}

crInfo &crInfo::operator=(const crInfo &other)
{
	if(this == &other)
		return *this;

	CRnum		= other.CRnum;
	timeStart 	= other.timeStart;
	timeStop	= other.timeStop;

	return *this;
}

bool crInfo::operator==(const crInfo &other)
{
	return this->CRnum == other.CRnum;
}

bool crInfo::operator>(const crInfo &other)
{
	return this->CRnum > other.CRnum;
}

bool crInfo::operator>=(const crInfo &other)
{
	return this->CRnum >= other.CRnum;
}

bool crInfo::operator<(const crInfo &other)
{
	return this->CRnum < other.CRnum;
}

bool crInfo::operator<=(const crInfo &other)
{
	return this->CRnum <= other.CRnum;
}

void crInfo::initNULL()
{
	CRnum		= 0;
	timeStart	= hcDate();
	timeStop	= hcDate();
}

void crInfo::clear()
{
	initNULL();
}

void crInfo::init(uint CRnum)
{
	clear();
	this->CRnum	= CRnum;
	timeStart 	= crInfoList::getStartDate(CRnum);
	timeStop 	= crInfoList::getStopDate(CRnum);
}

void crInfo::dump(uint indent) const
{
	stringstream ind;
	if(indent > 0) ind << setw(indent) << setfill(' ') << " ";
	cout << ind.str() << "Dumping crInfo:\n";
	cout << left;
	cout << ind.str() << setw(20) << setfill(' ') << "Carrington rot #: " 	<< CRnum						<< "\n";
	cout << ind.str() << setw(20) << setfill(' ') << "start time:" 			<< timeStart.toSpiceString() 	<< "\n";
	cout << ind.str() << setw(20) << setfill(' ') << "end time:" 			<< timeStop.toSpiceString()		<< "\n";
}

bool crInfo::exportBinary(std::ofstream &stream)
{
	bool retval = true;
	stream.write(reinterpret_cast<char*>(&CRnum), sizeof(uint));
	retval &= timeStart.exportBinary(stream);
	retval &= timeStop.exportBinary(stream);
	return retval;
}

bool crInfo::importBinary(std::ifstream &stream)
{
	bool retval = true;
	stream.read(reinterpret_cast<char*>(&CRnum), sizeof(uint));
	retval &= timeStart.importBinary(stream);
	retval &= timeStop.importBinary(stream);
	return retval;
}


//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//
//                          crInfloList
//
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

crInfoList::crInfoList()
{
    initNULL();
}

crInfoList::crInfoList(const crInfoList &other)
{
	initNULL();
    printf("ERROR! crInfoList cpy constructor not implemented!\n");
    exit(1);
}

crInfoList::~crInfoList()
{
    clear();
}

crInfoList &crInfoList::operator=(const crInfoList &other)
{
    if(this == &other)
        return *this;

    printf("ERROR crInfoList assignment op not implemented!\n");
    exit(1);

    return *this;
}

void crInfoList::clear()
{
    if(first != NULL)
    {
        crListElement *element0 = first;
        crListElement *element1;

        while(element0 != NULL)
        {
            element1 = element0->next;
            delete element0;
            element0 = element1;
        }
    }
    initNULL();
}

void crInfoList::initNULL()
{
    first   = NULL;
    last    = NULL;
}

void crInfoList::init()
{
    clear();
}

crListElement *crInfoList::getcrListElement(uint carRotNum)
{
    crListElement *element = first;

    while(element != NULL && element->CRnum < carRotNum)
        element = element->next;

    if(element == NULL)		        return NULL;
    if(element->CRnum == carRotNum) return element;
    return NULL;
}

bool crInfoList::appendObservation(uint carRotNum, originID origin, int pos)
{
    crListElement *element = getcrListElement(carRotNum);

    if(element != NULL)
        return element->appendObservatory(origin, pos);

    element = new crListElement(carRotNum);
    if(!element->appendObservatory(origin, pos))
    {
        delete element;
        return false;
    }

    if(first == NULL)
    {
        first   = element;
        last    = element;
        return true;
    }

    if(first->CRnum > carRotNum)
    {
        element->next   = first;
        first->prev     = element;
        first           = element;        
        return true;
    }

    crListElement *oldElement = first;

    while(true)
    {

        if(oldElement   == NULL)
        {
            element->prev       = last;
            element->prev->next = element;
            last                = element;
            return true;
        }

        if(oldElement->CRnum > carRotNum)
        {
            element->prev       = oldElement->prev;
            element->prev->next = element;
            element->next       = oldElement;
            oldElement->prev    = element;
            return true;

        }
        oldElement = oldElement->next;
    }

    return false;
}

bool crInfoList::initStaticMembers()
{
    clearStaticMembers();

    string filename = getFilename_crlist();

    if(!doesFileExist(filename))
    {
    	cerr << __FILE__ << ":" << __LINE__ << ": File '" << filename << "' does not exist.\n";
    	return false;
    }

    std::ifstream in(filename);

    char instrings[6][100];
    uint counter = 0;

    char line[1000];
    line[0] = '\t';
    while(line[0] != '\0' && line[0] != '\n')
    {
        line[0] = '\0';
        in.getline(line, 1000);
        ++counter;
    }

    in.close();
    numCRinList         = counter-1;
    listCRNum           = new int[counter];
    listJulianDayNum    = new long int[counter];
    listJulianDayFrac   = new double[counter];
    dateStart           = new hcDate[counter];
    in.open(filename);

    counter = 0;
    while(in >> instrings[0] >> instrings[1] >> instrings[2] >> instrings[3] >> instrings[4] >> instrings[5])
    {
        listCRNum[counter]  = atoi(instrings[0]);
        uint i = 0;
        char a = instrings[5][0];
        char tempStr[100] = "0\0";

        while(a != '\0')
        {
            if(a == '.')
                break;
            ++i;
            a = instrings[5][i];
        }
        strncpy(tempStr, instrings[5], i);
        listJulianDayNum[counter]   = atol(tempStr);

        instrings[5][i-1] = '0';
        strncpy(tempStr, instrings[5]+i-1, 10);
        listJulianDayFrac[counter]  = atof(tempStr);

        dateStart[counter].setFromJD(listJulianDayNum[counter], listJulianDayFrac[counter]);
        ++counter;
    }
    in.close();

    return true;
}

void crInfoList::clearStaticMembers(){

    delete [] listCRNum;
    delete [] listJulianDayNum;
    delete [] listJulianDayFrac;
    delete [] dateStart;
}

hcDate crInfoList::getStartDate(int carRotNum)
{
    hcDate retval;

    if(numCRinList == 0)
	{
    	printf("ERROR! crInfoList::getStartDate: This structure is not initialized! numCRinList==0\n");
    	return retval;
	}

    if(carRotNum < listCRNum[0])
    {
        retval.set(-1000000000, 0, 0, 0, 0, 0);
        return retval;
    }

    if(carRotNum > listCRNum[numCRinList-1])
    {
        retval.set(1000000000, 0, 0, 0, 0, 0);
        return retval;
    }

    for(uint i=0; i<numCRinList; ++i)
    {
        if(carRotNum == listCRNum[i])
        {
            retval = dateStart[i];
            return retval;
        }
    }

    printf("ERROR! crInfoList::getStartDate(): No date found!\n");
    return retval;
}

hcDate crInfoList::getStopDate(int carRotNum)
{
    hcDate retval;

    if(numCRinList == 0)
	{
		printf("ERROR! crInfoList::getStopDate: This structure is not initialized! numCRinList==0\n");
		return retval;
	}

    if(carRotNum < listCRNum[0])
    {
        retval.set(-1000000000, 0, 0, 0, 0, 0);
        return retval;
    }

    if(carRotNum > listCRNum[numCRinList-1])
    {
        retval.set(1000000000, 0, 0, 0, 0, 0);
        return retval;
    }

    for(uint i=0; i<numCRinList; ++i)
    {
        if(carRotNum == listCRNum[i])
        {
            if(i < numCRinList - 1)
                retval = dateStart[i+1];
            else
                retval = dateStart[i] + 25.38 * 86400 * hcDate::facSec;
            return retval;
        }
    }

    printf("ERROR! crInfoList::getStopDate(): No date found!\n");
    return hcDate();
}

int crInfoList::getCRnumber(const hcDate &date)
{
    if(date - getStartDate(listCRNum[0]) < 0)
        return std::numeric_limits<int>::min();

    if(getStartDate(listCRNum[numCRinList-1]) - date < 0)
        return std::numeric_limits<int>::max();

    uint i = 0;

    while(getStopDate(listCRNum[i]) - date < 0)
        ++i;

    return listCRNum[i];
}

void crInfoList::dump() const
{
    uint i=0;
    printf("Dumping crInfoList:\n");
    crListElement *element = first;

    while(element != NULL)
    {
        printf("\nElement %u:\n", ++i);
        element->dump();
        element = element->next;
    }

}

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//
//                          crListElement
//
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

crListElement::crListElement(int carRotNum)
{
    initNULL();
    init(carRotNum);
}

crListElement::crListElement(const crListElement &other)
{
	initNULL();
    printf("ERROR! crListElement cpy-constructor not implemented!\n");
    exit(1);
}

crListElement::~crListElement()
{
    clear();
}

crListElement &crListElement::operator=(const crListElement &other)
{
    if(this == &other)
        return *this;

    printf("ERROR! crListElement assignment op not implemented!\n");
    exit(1);

    return *this;
}

void crListElement::clear()
{
    delete [] WSO_solution;
    delete [] KPVT_solution;
    delete [] MDI_solution;
    delete [] MDIDAILY_solution;
    delete [] GONG_solution;
    delete [] HMI_solution;
    delete [] OWN_solution;
    initNULL();
}

void crListElement::initNULL()
{
    prev            	= NULL;
    next            	= NULL;
    WSO_solution    	= NULL;
    KPVT_solution   	= NULL;
    MDI_solution    	= NULL;
    MDIDAILY_solution   = NULL;
    GONG_solution   	= NULL;
    HMI_solution    	= NULL;
    OWN_solution		= NULL;
    numWSO          	= 0;
    numKPVT         	= 0;
    numMDI          	= 0;
    numMDIDAILY     	= 0;
    numGONG         	= 0;
    numHMI          	= 0;
    numOWN				= 0;
}

void crListElement::init(int carRotNum)
{
    clear();
    ((crInfo*)this)->init(carRotNum);
}

void crListElement::addElementToArray(int **array, uint &numEntries, uint newValue)
{
    int *oldarray   = *array;
    *array          = new int[numEntries + 1];
    for(uint i=0; i<numEntries; ++i)
        (*array)[i] = oldarray[i];
    (*array)[numEntries] = newValue;

    if(numEntries > 0)
        delete [] oldarray;

    ++numEntries;
}

bool crListElement::appendObservatory(originID origin, int pos)
{
    if(origin == ORIGIN_WSO)		        addElementToArray(&WSO_solution, 		numWSO, 		pos);
    else if(origin == ORIGIN_NSOKPVT)		addElementToArray(&KPVT_solution, 		numKPVT, 		pos);
    else if(origin == ORIGIN_SOHOMDI)		addElementToArray(&MDI_solution, 		numMDI, 		pos);
    else if(origin == ORIGIN_SOHOMDI_DAILY)	addElementToArray(&MDIDAILY_solution, 	numMDIDAILY, 	pos);
    else if(origin == ORIGIN_NSOGONG)		addElementToArray(&GONG_solution, 		numGONG, 		pos);
    else if(origin == ORIGIN_SDOHMI)		addElementToArray(&HMI_solution, 		numHMI, 		pos);
    else if(origin == ORIGIN_OWN)			addElementToArray(&OWN_solution, 		numHMI, 		pos);
    else
    {
        printf("ERROR! crListElement::appendObservatory: originID (%u) not supported!\n", origin);
        printf("WSO: %u, KPVT: %u, MDI: %u, UNKNOWN: %u\n", ORIGIN_WSO, ORIGIN_NSOKPVT, ORIGIN_SOHOMDI, ORIGIN_UNKNOWN);
        return false;
    }
    return true;
}

void crListElement::dump()
{
    printf("Dumping crListElement %i\n", CRnum);
    if(prev == NULL)        printf("Previous: NULL");
    else			        printf("Previous: %i", prev->CRnum);
    if(next == NULL)        printf("\tNext: NULL\n");
    else			        printf("\tNext: %i\n", next->CRnum);

    printf("WSO (%u):\t", numWSO);
    for(uint i=0; i<numWSO; ++i)	    printf("%i\t", WSO_solution[i]);

    printf("\nKPVT (%u):\t", numKPVT);
    for(uint i=0; i<numKPVT; ++i)		printf("%i\t", KPVT_solution[i]);

    printf("\nMDI (%u):\t", numMDI);
    for(uint i=0; i<numMDI; ++i)		printf("%i\t", MDI_solution[i]);

    printf("\nMDIDAILY (%u):\t", numMDIDAILY);
	for(uint i=0; i<numMDIDAILY; ++i)	printf("%i\t", MDIDAILY_solution[i]);

    printf("\nGONG (%u):\t", numGONG);
    for(uint i=0; i<numGONG; ++i)       printf("%i\t", GONG_solution[i]);

    printf("\nHMI (%u):\t", numHMI);
    for(uint i=0; i<numHMI; ++i)		printf("%i\t", HMI_solution[i]);

    printf("\nOWN (%u):\t", numOWN);
	for(uint i=0; i<numOWN; ++i)		printf("%i\t", OWN_solution[i]);

    printf("\n\n");
}

