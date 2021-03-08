#ifndef CARROTINFO_H
#define CARROTINFO_H

#include "engine/hcTime.h"

#include "src/enum.h"


/*! \brief information on specific Carrington rotation
 */
class crInfo{
public:

	uint CRnum;				/*!< \brief number of Carrington Rotation							*/
	hcDate timeStart;		/*!< \brief start date of Carrington Rotation						*/
	hcDate timeStop;		/*!< \brief stop date of Carrington Rotation						*/

	crInfo(uint crNum = 0);			/*!< \brief std constructor									*/
	crInfo(const crInfo &other);	/*!< \brief cpy constructor									*/
	virtual ~crInfo();				/*!< \brief destructor										*/

	crInfo &operator=(const crInfo &other);		/*!< \brief assignment operator					*/
	bool operator==(const crInfo &other);		/*!< \brief comparison operator					*/
	bool operator>(const crInfo &other);		/*!< \brief comparison operator					*/
	bool operator<(const crInfo &other);		/*!< \brief comparison operator					*/
	bool operator>=(const crInfo &other);		/*!< \brief comparison operator					*/
	bool operator<=(const crInfo &other);		/*!< \brief comparison operator					*/

	void initNULL();
	void clear();
	void init(uint CRnum);

	void dump() const;
};

class crListElement;

/*! \brief list of carrington rotation numbers and corresponding start times
 */
class crInfoList{
public:
    static uint numCRinList;            /*!< \brief number of carrington rotations stored in list           */
    static int *listCRNum;              /*!< \brief list of carrington rotation numbers                     */
    static long int *listJulianDayNum;  /*!< \brief list of corresponding start times (JD integer part)     */
    static double *listJulianDayFrac;   /*!< \brief list of corresponding start times (JD fractional part)  */
    static hcDate *dateStart;           /*!< \brief list of start dates                                     */
    crListElement *first;				/*!< \brief first element in list									*/
	crListElement *last;				/*!< \brief last element in list									*/


    crInfoList();
    crInfoList(const crInfoList &other);
    ~crInfoList();

    crInfoList &operator=(const crInfoList &other);

    void clear();
    void initNULL();
    void init();

    crListElement *getcrListElement(int carRotNum);
        /*!< \brief returns the list element representing the specified carRotNum or NULL, if non-existent  */

    bool appendObservation(int carRotNum, originID origin, int pos);


    static void initStaticMembers(const char *filename);
    static void clearStaticMembers();

    static hcDate getStartDate(int carRotNum);
        /*!< \brief computes start date of specific Carrington rotation                                     */

    static hcDate getStopDate(int carRotNum);
        /*!< \brief computes stop date of specific Carrington rotation                                      */

    static int getCRnumber(const hcDate &date);
        /*!< \brief computes the Carrington rotation number of a specific date                              */

    void dump() const ;
};

/*! information on which instrument has been used to compute the magnetic field of Carrington Rotation carRotNum
 *
 *  the integer values XXX_solution are either -1 if no solution has been computed using that instrument
 *  or the positional number in the solutions-set of HelioMagfield
 */
class crListElement : public crInfo{
public:

    crListElement *prev;
    crListElement *next;

    int *WSO_solution;
    int *KPVT_solution;
    int *MDI_solution;
    int *MDIDAILY_solution;
    int *GONG_solution;
    int *HMI_solution;
    int *OWN_solution;

    uint numWSO;
    uint numKPVT;
    uint numMDI;
    uint numMDIDAILY;
    uint numGONG;
    uint numHMI;
    uint numOWN;

    crListElement(int carRotNum=0);             /*!< \brief std-constructor         */
    crListElement(const crListElement &other);  /*!< \brief cpy constructor         */
    ~crListElement();

    crListElement &operator=(const crListElement &other);

    void clear();
    void initNULL();
    void init(int carRotNum);

    void addElementToArray(int **array, uint &numEntries, uint newValue);
        /*!< \brief expands array by one and appends newValue at the end                        */

    bool appendObservatory(originID origin, int pos);
        /*!< \brief appends the position in solutions-set of the new observation to element     */

    void dump();

};

#endif // CARROTINFO_H
