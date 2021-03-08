#ifndef CONTAINER_H
#define CONTAINER_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>

typedef unsigned int uint;

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//
//                          hcSet<T>
//
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

template <class T>
class hcSet{
public:

    uint numMaxSlots;
    uint numOccupiedSlots;
    uint numFreeSlots;
    uint firstFreeSlot;
    uint lastOccSlot;
    T **elements;
    bool *occupied;

    hcSet(uint numMaxSlots = 4);
    hcSet(const hcSet<T> &other);
    virtual ~hcSet();

    hcSet<T> &operator=(const hcSet<T> &other);
    T &operator[](uint i) const;

    void clear();
    void initNULL();
    void init(uint numMaxSlots);

    virtual void expand(uint numNewElements);
    virtual int findElement(const T &object);
    virtual bool hasElement(const T &object);


    virtual bool isSlotOccupied(uint n) const;
    	/*!< \brief is slot n occupied?														*/

    int appendElement(T &object);
    	/*!< \brief appends given object to this set										*/

    bool appendElementAtPos(T &object, uint pos);
        /*!< \brief appends element at a specific position if possible                      */

    virtual int removeElement(T &object);
    	/*!< \brief remove given object from this set										*/

    T &elementAt(uint num);
        /*!< \brief returns reference to element at position num                            */

    void dump() const;

};


template <class T>
hcSet<T>::hcSet(uint numMaxSlots){

    initNULL();
    init(numMaxSlots);
}

template <class T>
hcSet<T>::hcSet(const hcSet<T> &other){

    initNULL();
    *this = other;
}

template <class T>
hcSet<T>::~hcSet(){

    clear();
}

template <class T>
hcSet<T> &hcSet<T>::operator=(const hcSet<T> &other)
{
	if(this == &other)
		return *this;

	clear();
	init(other.numMaxSlots);

	for(uint i=0;i<numMaxSlots;++i)
	{
		this->occupied[i] = other.occupied[i];
		this->elements[i] = other.elements[i];
	}

	this->numOccupiedSlots  = other.numOccupiedSlots;
	this->numFreeSlots      = other.numFreeSlots;
	this->firstFreeSlot     = other.firstFreeSlot;
	this->lastOccSlot       = other.lastOccSlot;

	return *this;
}

template <class T>
T &hcSet<T>::operator[](uint n) const
{
	if(n<numMaxSlots && occupied[n])
		return *elements[n];

	// TODO: This is a really ugly hack! Get some better ideas you dumbass!!
    std::cerr << "hcSet<T>::operator[]: n >= numMaxSlots (" << n << ">=" << numMaxSlots << ") or slot["<< n << "] not occupied!\nMemory leak introduced!\n";
    exit(1);
	T *newT = new T();
	return *newT;
}


template <class T>
void hcSet<T>::clear()
{
    if(elements != NULL)
    {
        for(uint i=0; i<numMaxSlots; ++i)
            if(occupied[i])
                removeElement(*elements[i]);
        delete [] elements;
        delete [] occupied;
    }

    initNULL();
}

template <class T>
void hcSet<T>::initNULL()
{
	numMaxSlots			= 0;
	numOccupiedSlots	= 0;
	numFreeSlots		= 0;
	firstFreeSlot		= 0;
	lastOccSlot			= 0;

    elements = NULL;
    occupied = NULL;
}

template <class T>
void hcSet<T>::init(uint numMaxSlots){

    clear();
    elements	= new T*[numMaxSlots];
    occupied	= new bool[numMaxSlots];

    for(uint i=0;i<numMaxSlots;++i)
    {
        elements[i] = NULL;
        occupied[i] = false;
    }

    this->numMaxSlots       = numMaxSlots;
    this->numOccupiedSlots  = 0;
    this->numFreeSlots      = numMaxSlots;
    this->firstFreeSlot     = 0;
    this->lastOccSlot       = 0;
}

template <class T>
void hcSet<T>::expand(uint numNewElements)
{
#ifdef TODO
    printf("INFO: hcSet::expand: change set size from %u to %u.\n", numMaxSlots, numMaxSlots + numNewElements);
    printf("INFO: If you see this message too often, your program is probably in need of optimization\n");
    printf("INFO: meaning that this procedure is inefficient)\n");
#endif

    T **tempElementStorage  = new T*[numMaxSlots + numNewElements];
	bool *tempOccupied      = new bool[numMaxSlots + numNewElements];


    for(uint i=0;i<numMaxSlots;++i)
    {
        tempElementStorage[i]   = elements[i];
        tempOccupied[i]         = occupied[i];
    }

    for(uint i=numMaxSlots;i<numMaxSlots + numNewElements;++i)
    {
        tempElementStorage[i]   = NULL;
        tempOccupied[i]         = false;
    }

    numMaxSlots     += numNewElements;
    numFreeSlots    += numNewElements;

    delete [] elements;
    delete [] occupied;

    elements    = tempElementStorage;
    occupied    = tempOccupied;

    while(occupied[++firstFreeSlot] && firstFreeSlot < numMaxSlots);

    if(firstFreeSlot > lastOccSlot)
    {
        lastOccSlot = firstFreeSlot;
        if(!occupied[lastOccSlot])
            while(lastOccSlot > 0 && !occupied[--lastOccSlot]);
    }
}

template <class T>
int hcSet<T>::findElement(const T &object){

    for(uint i=0;i <= lastOccSlot;++i)
        if(occupied[i])
        {
            if(elements[i] == &object)
                return i;
        }
    return -1;
}

template <class T>
bool hcSet<T>::hasElement(const T &object){

    if(findElement(object) > -1)
        return true;

    return false;
}

template <class T>
bool hcSet<T>::isSlotOccupied(uint n) const{

	if(n >= numMaxSlots)
	{
        std::cerr << "hcSet<T>::isSlotOccupied: n > numMaxSlots (" << n << ">" << numMaxSlots <<")\n";
		return false;
	}

	return occupied[n];
}

/*!	\returns position where object is stored
 *
 */
template <class T>
int hcSet<T>::appendElement(T &object)
{
    uint retval 		= firstFreeSlot;
    elements[retval] 	= &object;
    occupied[retval] 	= true;

    if ((numFreeSlots > 0)  && --numFreeSlots > 0)
        while(occupied[++firstFreeSlot] && firstFreeSlot < numMaxSlots);
    else
    	expand(numMaxSlots);	// TODO: catch memory errors

    if(firstFreeSlot > lastOccSlot)
    {
        lastOccSlot = firstFreeSlot;
        if(!occupied[lastOccSlot])
            while(lastOccSlot > 0 && !occupied[--lastOccSlot]);
    }

    ++numOccupiedSlots;

    return retval;
}

#include <iostream>

template <class T>
bool hcSet<T>::appendElementAtPos(T &object, uint pos)
{
    if(pos >= numMaxSlots)
    {
        printf("ERROR! hcSet::appendElementAtPos: requested position (%u) out of bounds (%u)!\n", pos, numMaxSlots);
        return false;
    }

    if(occupied[pos])
    {
        printf("ERROR! hcSet::appendElementAtPos: reqested position (%u) already occupied!\n", pos);
        return false;
    }

    elements[pos] = &object;
    occupied[pos] = true;

    if(pos == firstFreeSlot)
    {
        if ((numFreeSlots > 0)  && --numFreeSlots > 0)
            while(occupied[++firstFreeSlot] && firstFreeSlot < numMaxSlots);

        if(firstFreeSlot > lastOccSlot)
        {
            lastOccSlot = firstFreeSlot;
            if(!occupied[lastOccSlot])
                while(lastOccSlot > 0 && !occupied[--lastOccSlot]);
        }
    }
    else
    {
        --numFreeSlots;
        if(pos > lastOccSlot)
            lastOccSlot = pos;
    }

    ++numOccupiedSlots;
    return true;
}

/*! \returns position where object has been found and removed or -1 on failure
 *
 */
template <class T>
int hcSet<T>::removeElement(T &object){

    for(uint i=0;i<=lastOccSlot;++i)
    {
        if(occupied[i])
        {
            if(&object == elements[i])
            {                
            	//void *q = (void*)this->elements[i];
            	//uintptr_t s = (uintptr_t)q;
            	//printf("hcSet::removeElement %llx\n", s);fflush(stdout);
                elements[i] = NULL;
                occupied[i] = false;
                ++numFreeSlots;
                --numOccupiedSlots;
                if(i < firstFreeSlot)
                    firstFreeSlot = i;

                if(lastOccSlot == i)
                    while(lastOccSlot > 0 && !occupied[--lastOccSlot]);


                return i;
            }
        }
    }
    return -1;
}

template <class T>
T &hcSet<T>::elementAt(uint num){

    return this->operator[](num);
}

template <class T>
void hcSet<T>::dump() const{

    printf("--- Dumping hcSet/hcSetStorage:\n");
    printf("# of slots: %u\n", this->numMaxSlots);
    printf("# of occupied slots: %u\n", this->numOccupiedSlots);
    printf("lastOccupiedSlot: %u\n", this->lastOccSlot);
}

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//
//                          hcSetStorage<T>
//
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

template <class T>
class hcSetStorage : public hcSet<T>{
public:

	hcSetStorage(uint numMaxSlots = 1);			/*!< \brief std constructor										*/
    hcSetStorage(const hcSetStorage &other);	/*!< \brief cpy constructor										*/
    virtual ~hcSetStorage();

    hcSetStorage<T> &operator=(const hcSetStorage<T> &other);

    virtual int removeElement(T &object);       /*!< \brief removes and destroys an element             */
    bool extractElement(T &object);              /*!< \brief removes an element and transfers ownership  */
};

template <class T>
hcSetStorage<T>::hcSetStorage(uint numMaxSlots) : hcSet<T>(numMaxSlots){}

template <class T>
hcSetStorage<T>::hcSetStorage(const hcSetStorage<T> &other) :
    hcSet<T>(other.numMaxSlots)
{
    *this = other;
}


template <class T>
hcSetStorage<T>::~hcSetStorage(){

    this->clear();
}

template <class T>
hcSetStorage<T> &hcSetStorage<T>::operator=(const hcSetStorage<T> &other){

	if(this == &other)
		return *this;

	hcSet<T>::clear();
	this->init(other.numMaxSlots);

	for(uint i=0;i<this->numMaxSlots;++i)
	{
		this->occupied[i] = other.occupied[i];
        if(this->occupied[i])
            this->elements[i] = new T(*other.elements[i]);
        else
            this->elements[i] = NULL;
	}

	this->numOccupiedSlots  = other.numOccupiedSlots;
	this->numFreeSlots      = other.numFreeSlots;
	this->firstFreeSlot     = other.firstFreeSlot;
	this->lastOccSlot       = other.lastOccSlot;

	return *this;
}

template <class T>
int hcSetStorage<T>::removeElement(T &object)
{
    for(uint i=0;i<=this->lastOccSlot;++i)
        if(this->occupied[i])
        {            
            if(&object == this->elements[i])
            {
            	//void *q = (void*)this->elements[i];
            	//uintptr_t s = (uintptr_t)q;
            	//printf("hcSetStorage::removeElement %llx\n", s);fflush(stdout);
                delete this->elements[i];

                this->elements[i] = NULL;
                this->occupied[i] = false;
                ++this->numFreeSlots;
                --this->numOccupiedSlots;
                if(i < this->firstFreeSlot)
                    this->firstFreeSlot = i;

                if(this->lastOccSlot == i)
                    while(hcSet<T>::lastOccSlot > 0 && !this->occupied[--hcSet<T>::lastOccSlot]);

                return i;
            }
        }

    return -1;
}

template <class T>
bool hcSetStorage<T>::extractElement(T &object){

    for(uint i=0;i<=this->lastOccSlot;++i)
        if(this->occupied[i] && &object == this->elements[i])
        {
            this->elements[i] = NULL;
            this->occupied[i] = false;
            ++this->numFreeSlots;
            --this->numOccupiedSlots;
            if(i < this->firstFreeSlot)
                this->firstFreeSlot = i;

            if(this->lastOccSlot == i)
                while(hcSet<T>::lastOccSlot > 0 && !this->occupied[--hcSet<T>::lastOccSlot]);

            return true;
        }

    return false;
}

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//
//                          hcDualContainer<T>
//
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

template <class T1, class T2>
class hcDualSet : public hcSet<T1>{
public:

    T2 **otherElements;

    hcDualSet(uint numMaxSlots = 1);
    ~hcDualSet();

    void clear();
    void initNULL();
    void init(uint numMaxSlots);

    T2* getOtherObject(const T1 &object);
    int appendElement(T1 &object, T2 &otherObject);
    int removeElement(T1 &object);
};

template <class T1, class T2>
hcDualSet<T1, T2>::hcDualSet(uint numMaxSlots) : hcSet<T1>(numMaxSlots)
{
    initNULL();
    init(numMaxSlots);
}

template <class T1, class T2>
hcDualSet<T1, T2>::~hcDualSet()
{
    clear();
}

template <class T1, class T2>
void hcDualSet<T1, T2>::clear()
{
    hcSet<T1>::clear();
    if(otherElements != NULL)
    	delete [] otherElements;

    initNULL();
}

template <class T1, class T2>
void hcDualSet<T1, T2>::initNULL()
{
    otherElements = NULL;
}

template <class T1, class T2>
void hcDualSet<T1, T2>::init(uint numMaxSlots)
{
    hcSet<T1>::init(numMaxSlots);
    otherElements = new T2*[numMaxSlots];

    for(uint i=0;i<numMaxSlots;++i)
        otherElements[i] = NULL;
}

template <class T1, class T2>
T2* hcDualSet<T1, T2>::getOtherObject(const T1 &object)
{
    int pos = hcSet<T1>::findElement(object);

    if(pos != -1)
        return otherElements[pos];

    return NULL;
}

template <class T1, class T2>
int hcDualSet<T1, T2>::appendElement(T1 &object, T2 &otherObject)
{
    int pos = hcSet<T1>::appendElement(object);

    if(pos == -1)
    {
        printf("ERROR! hcDualContainer:appendObject(): Container full!\n");
        return -1;
    }

    otherElements[pos] = &otherObject;

    return pos;
}

template <class T1, class T2>
int hcDualSet<T1, T2>::removeElement(T1 &object)
{
    int pos = removeElement(object);

    if(pos == -1)
    {
        printf("WARNING! hyDualContainer::removeObject(): Object not found!\n");
        return -1;
    }

    otherElements[pos] = NULL;

    return pos;
}

#endif // CONTAINER_H
