#ifndef PFSSSOLUTIONINFO_H
#define PFSSSOLUTIONINFO_H

#include "src/enum.h"
#include "src/synPhotMagfield.h"
#include "src/laplaceSolver.h"

typedef unsigned int uint;

class PFSSsolutionInfo : public SynopticInfo{
public:

    uint sizeofFloat;							/*!< \brief floating point type used for computation 4=float, 8=double		*/
    modelID model;								/*!< \brief model employed for magnetic computation							*/
    methodID method;							/*!< \brief method applied for (PFSS) computation 							*/
    groupID group;								/*!< \brief working group which computed this solution						*/
    hcFloat rss;								/*!< \brief source surface radius 										(m)	*/
    hcFloat ell;								/*!< \brief ellipticity of source surface (if method==METH_ELLIPTICAL)		*/
    uint orderSHC;								/*!< \brief maximum principal order of SHC solution (if method==METH_SHC)	*/
    uint numR;									/*!< \brief number of grid points in radial direction						*/
    uint numTheta;								/*!< \brief number of grid points in meridional direction					*/
    uint numPhi;								/*!< \brief number of grid points in zonal direction						*/
    hcDate dateComputed;						/*!< \brief date when the solution was computed								*/
    uint computationTime;						/*!< \brief time in seconds it took to compute								*/
    uint solutionSteps;							/*!< \brief number of solution steps to fall below accuracy threshold		*/

    PFSSsolutionInfo();
    	/*!< \brief std constructor         																				*/

    PFSSsolutionInfo(const PFSSsolutionInfo &other);
    	/*!< \brief cpy constructor         																				*/

    PFSSsolutionInfo(	SynopticInfo synInf, uint sizeofFloat, modelID model, methodID method, groupID group,
        				hcFloat rss, hcFloat ell, uint orderSHC, uint numR, uint numTheta, uint numPhi);
    	/*!< \brief constructor																								*/

    virtual ~PFSSsolutionInfo(){}
    	/*!< \brief destructor																								*/

    PFSSsolutionInfo &operator=(const PFSSsolutionInfo &other);
    	/*!< \brief assignment operator																						*/

    bool operator==(const PFSSsolutionInfo &other);
    	/*!< \brief comparison operator																						*/

    bool operator>(const PFSSsolutionInfo &other);
    	/*!< \brief comparison operator																						*/

    bool operator<(const PFSSsolutionInfo &other);
    	/*!< \brief comparison operator																						*/

    void initNULL();
    void clear();
    void init(	SynopticInfo synInf, uint sizeofFloat, modelID model, methodID method, groupID group,
    			hcFloat rss, hcFloat ell, uint orderSHC, uint numR, uint numTheta, uint numPhi);

    bool exportBinary(ofstream &stream);
    	/*!< \brief export instance to binary stream																		*/

    bool importBinary(ifstream &stream);
    	/*!< \brief imports instance from binary stream																		*/

    string toString() const;
    	/*!< \brief returns information on this instance in string															*/

    void dump(uint indent=0) const;
    	/*!< \brief dumps information on this instance to stdout with optional indendtation									*/
};

#endif
