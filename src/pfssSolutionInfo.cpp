#include "src/pfssSolutionInfo.h"
#include "src/ellipticalGrid.h"

#include <iostream>
#include <iomanip>

#include "boost/regex.hpp"

using namespace boost;

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//
//                          PFSSsolutionInfo
//
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

PFSSsolutionInfo::PFSSsolutionInfo()
{
    initNULL();
}

PFSSsolutionInfo::PFSSsolutionInfo(const PFSSsolutionInfo &other) :
		SynopticInfo(other)
{
	initNULL();
	operator =(other);
}

PFSSsolutionInfo::PFSSsolutionInfo(SynopticInfo synInf, uint sizeofFloat, modelID model, methodID method, groupID group,
									hcFloat rss, hcFloat ell, uint orderSHC, uint numR, uint numTheta, uint numPhi)
{
	initNULL();
	init(synInf, sizeofFloat, model, method, group, rss, ell, orderSHC, numR, numTheta, numPhi);
}

PFSSsolutionInfo &PFSSsolutionInfo::operator=(const PFSSsolutionInfo &other)
{
    if(this == &other)
        return *this;

    this->sizeofFloat		= other.sizeofFloat;
    this->model				= other.model;
    this->method			= other.method;
    this->group				= other.group;
    this->rss				= other.rss;
    this->ell				= other.ell;
    this->orderSHC			= other.orderSHC;
    this->numR				= other.numR;
    this->numTheta			= other.numTheta;
    this->numPhi			= other.numPhi;
    this->dateComputed		= other.dateComputed;
    this->computationTime 	= other.computationTime;

    ((SynopticInfo*)this)->operator=(other);

    return *this;
}

bool PFSSsolutionInfo::operator==(const PFSSsolutionInfo &other)
{
	bool retval = true;
	retval &= SynopticInfo::operator==(other);
	retval &= rss			== other.rss;
	retval &= ell			== other.ell;
	retval &= orderSHC		== other.orderSHC;
	retval &= numR			== other.numR;
	retval &= numTheta		== other.numTheta;
	retval &= numPhi		== other.numPhi;
	return retval;
}

bool PFSSsolutionInfo::operator>(const PFSSsolutionInfo &other)
{
	if(		CRnum	> other.CRnum)		return true;
	else if(CRnum	< other.CRnum)		return false;

	if(		rss > other.rss)			return true;
	else if(rss < other.rss)			return false;

	if(		ell > other.ell)			return true;
	else if(ell < other.ell)			return false;

	if(		method	> other.method)		return true;
	else if(method	< other.method)		return false;

	if(		orderSHC > other.orderSHC)	return true;
	else if(orderSHC < other.orderSHC)	return false;

	if(		instrument > other.instrument)	return true;
	else if(instrument < other.instrument)	return false;

	if(		maxSinLat > other.maxSinLat)	return true;
	else if(maxSinLat < other.maxSinLat)	return false;

	if(		dailyID > other.dailyID)	return true;
	else if(dailyID < other.dailyID)	return false;

	if(		numR > other.numR)			return true;
	else if(numR < other.numR)			return false;

	if(		numTheta > other.numTheta)	return true;
	else if(numTheta < other.numTheta)	return false;

	if(		numPhi > other.numPhi)		return true;
	else if(numPhi < other.numPhi)		return false;

	return false;
}

bool PFSSsolutionInfo::operator<(const PFSSsolutionInfo &other)
{
	return !(operator>(other) || operator==(other));
}

void PFSSsolutionInfo::initNULL()
{
    method				= METH_NUMERIC;
    sizeofFloat			= 0;
    rss					= 0.0;
    ell					= 0.0;
    orderSHC			= 0;
    numR                = 0;
    numTheta            = 0;
    numPhi              = 0;
    dateComputed		= hcDate(1970, 0, 0, 0, 0, 0);
    computationTime		= 0;
    ((SynopticInfo*)this)->initNULL();
}

void PFSSsolutionInfo::clear()
{
	initNULL();
}

void PFSSsolutionInfo::init(SynopticInfo synInf, uint sizeofFloat, modelID model, methodID method, groupID group,
							hcFloat rss, hcFloat ell, uint orderSHC, uint numR, uint numTheta, uint numPhi)
{
	((SynopticInfo*)this)->operator =(synInf);
	this->sizeofFloat	= sizeofFloat;
	this->model			= model;
	this->method		= method;
	this->group			= group;
	this->rss			= rss;
	this->ell			= ell;
	this->orderSHC		= orderSHC;
	this->numR			= numR;
	this->numTheta		= numTheta;
	this->numPhi		= numPhi;
}

bool PFSSsolutionInfo::exportBinary(std::ofstream &stream)
{
	bool retval = true;
	stream.write(reinterpret_cast<char*>(&sizeofFloat),		sizeof(uint));
	SynopticInfo::exportBinary(stream);
	stream.write(reinterpret_cast<char*>(&model),			sizeof(int));
	stream.write(reinterpret_cast<char*>(&method),			sizeof(int));
	stream.write(reinterpret_cast<char*>(&group),			sizeof(int));
	stream.write(reinterpret_cast<char*>(&rss),				sizeofFloat);
	stream.write(reinterpret_cast<char*>(&ell),				sizeofFloat);
	stream.write(reinterpret_cast<char*>(&orderSHC), 		sizeof(uint));
	stream.write(reinterpret_cast<char*>(&numR), 			sizeof(uint));
	stream.write(reinterpret_cast<char*>(&numTheta), 		sizeof(uint));
	stream.write(reinterpret_cast<char*>(&numPhi), 			sizeof(uint));
	retval &= dateComputed.exportBinary(stream);
	stream.write(reinterpret_cast<char*>(&computationTime), sizeof(uint));
	stream.write(reinterpret_cast<char*>(&solutionSteps), 	sizeof(uint));
	return retval;
}

bool PFSSsolutionInfo::importBinary(std::ifstream &stream)
{
	bool retval = true;
	uint sizeoffloat	= 0;
	stream.read(reinterpret_cast<char*>(&sizeoffloat),	sizeof(uint));
	if(sizeoffloat != 4 && sizeoffloat != 8) return false;
	char *tempFloat		= sizeoffloat == 4 ? reinterpret_cast<char*>(new float()): reinterpret_cast<char*>(new double());
	SynopticInfo::importBinary(stream);
	stream.read(reinterpret_cast<char*>(&model),			sizeof(int));
	stream.read(reinterpret_cast<char*>(&method),			sizeof(int));
	stream.read(reinterpret_cast<char*>(&group),			sizeof(int));
	stream.read(reinterpret_cast<char*>(tempFloat),	sizeoffloat); rss = sizeoffloat==4 ? *(reinterpret_cast<float*>(tempFloat)) : *(reinterpret_cast<double*>(tempFloat));
	stream.read(reinterpret_cast<char*>(tempFloat),	sizeoffloat); ell = sizeoffloat==4 ? *(reinterpret_cast<float*>(tempFloat)) : *(reinterpret_cast<double*>(tempFloat));
	stream.read(reinterpret_cast<char*>(&orderSHC), 		sizeof(uint));
	stream.read(reinterpret_cast<char*>(&numR), 			sizeof(uint));
	stream.read(reinterpret_cast<char*>(&numTheta), 		sizeof(uint));
	stream.read(reinterpret_cast<char*>(&numPhi), 			sizeof(uint));
	retval &= dateComputed.importBinary(stream);
	stream.read(reinterpret_cast<char*>(&computationTime), sizeof(uint));
	stream.read(reinterpret_cast<char*>(&solutionSteps), 	sizeof(uint));
	delete tempFloat;
	return retval;
}

string PFSSsolutionInfo::toString() const
{
	stringstream retval;
	retval.precision(4);
	retval << setw(4) 	<< setfill('0') << CRnum						 		<< " ";
	retval << setw(10) 	<< setfill(' ') << getStringFromMethodID(method) 		<< " ";
	retval << setw(10) 	<< setfill(' ') << getStringFromOriginID(instrument)	<< " ";
	retval << fixed 	<< rss/r_sol											<< " ";
	retval << ell 																<< " ";
	retval << setw(4)	<< orderSHC 											<< " ";
	retval << setw(4) 	<< setfill(' ') << numR									<< " ";
	retval << setw(4) 	<< setfill(' ') << numTheta								<< " ";
	retval << setw(4) 	<< setfill(' ') << numPhi								<< " ";
	retval << maxSinLat															<< " ";
	return retval.str();
}

void PFSSsolutionInfo::dump(uint indent) const
{
	stringstream ind;
	if(indent > 0) ind << setw(indent) << setfill(' ') << " ";
	cout << ind.str() << "Dumping PFSSsolutionInfo:\n";
	SynopticInfo::dump(indent+1);
	cout << left;
    cout << ind.str() << setw(20) << setfill(' ') << "float:" <<(sizeofFloat==4?"single":(sizeofFloat==8?"double":(sizeofFloat==16?"long double":"unknown"))) << "\n";
    cout << ind.str() << setw(20) << setfill(' ') << "Model:"			<< getStringFromModelID(model)		<< "\n";
    cout << ind.str() << setw(20) << setfill(' ') << "Group:"			<< getStringFromGroupID(group)		<< "\n";
    cout << ind.str() << setw(20) << setfill(' ') << "Method:"			<< getStringFromMethodID(method)	<< "\n";
    cout << ind.str() << setw(20) << setfill(' ') << "rss:"				<< rss / r_sol						<< "\n";
    cout << ind.str() << setw(20) << setfill(' ') << "ell:"				<< ell								<< "\n";
    cout << ind.str() << setw(20) << setfill(' ') << "orderSHC:"		<< orderSHC							<< "\n";
    cout << ind.str() << setw(20) << setfill(' ') << "numR:"			<< numR								<< "\n";
    cout << ind.str() << setw(20) << setfill(' ') << "numT:"			<< numTheta							<< "\n";
    cout << ind.str() << setw(20) << setfill(' ') << "numP:"			<< numPhi							<< "\n";
    cout << ind.str() << setw(20) << setfill(' ') << "dailyID:"			<< dailyID							<< "\n";
    cout << ind.str() << setw(20) << setfill(' ') << "compTime:"		<< computationTime					<< "\n";
    cout << ind.str() << setw(20) << setfill(' ') << "timestamp:"		<< dateComputed.toSpiceString()		<< "\n";
    fflush(stdout);
}
