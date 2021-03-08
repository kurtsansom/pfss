#ifndef SORTEDLIST_H
#define SORTEDLIST_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
using namespace std;

typedef unsigned int uint;

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//
//                          hcSortedList<T>
//
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
template<class T>
class hcSortedListStorage;

template <class T>
class hcSortedList{
public:

    uint numElements;
    uint numMaxSlots;
    T **elements;

    hcSortedList(uint numMaxSlots = 128);
    hcSortedList(const hcSortedList<T> &other);
    virtual ~hcSortedList();

    hcSortedList<T> &operator=(const hcSortedList<T> &other);
    hcSortedList<T> &operator=(const hcSortedListStorage<T> &other);
    T &operator[](uint i) const;

    void clear();
    void initNULL();
    void init(uint numMaxSlots=1024);

    virtual void expand(uint numNewElements);
    virtual int findElement(const T &object);
    virtual bool hasElement(const T &object);

    uint insertElement(T &object);
    	/*!< \brief appends given object to this list										*/

    virtual int removeElement(T &object);
    	/*!< \brief remove given object from this set										*/

    bool sort();

    void dump() const;

};


template <class T>
hcSortedList<T>::hcSortedList(uint numMaxSlots)
{
    initNULL();
    init(numMaxSlots);
}

template <class T>
hcSortedList<T>::hcSortedList(const hcSortedList<T> &other)
{
    initNULL();
    *this = other;
}

template <class T>
hcSortedList<T>::~hcSortedList()
{
    clear();
}

template <class T>
hcSortedList<T> &hcSortedList<T>::operator=(const hcSortedList<T> &other)
{
	if(this == &other)
		return *this;

	clear();
	init(other.numMaxSlots);

	for(uint i=0;i<numElements;++i)
		this->elements[i] = other.elements[i];

	return *this;
}

template <class T>
hcSortedList<T> &hcSortedList<T>::operator=(const hcSortedListStorage<T> &other)
{
	clear();
	init(other.numMaxSlots);

	for(uint i=0;i<numElements;++i)
		this->elements[i] = other.elements[i];

	return *this;
}

template <class T>
T &hcSortedList<T>::operator[](uint n) const
{
	if(n<numElements)
		return *elements[n];

	// TODO: This is a really ugly hack! Get some better ideas you dumbass!!
    printf("ERROR! hcSortedList<T>::operator[]: n >= numElements (%u >= %u)!\nMemory leak introduced!\n", n, numElements);
    fflush(stdout);
    exit(1);
	T *newT = new T();
	return *newT;
}

template <class T>
void hcSortedList<T>::clear()
{
    if(elements != NULL)
    {
    	while(this->numElements>0)
			removeElement(*elements[0]);
        delete [] elements;
    }

    initNULL();
}

template <class T>
void hcSortedList<T>::initNULL()
{
	numMaxSlots	= 0;
	numElements	= 0;
    elements 	= NULL;
}

template <class T>
void hcSortedList<T>::init(uint numMaxSlots){

    clear();
    elements	= new T*[numMaxSlots];

    for(uint i=0;i<numMaxSlots;++i)
        elements[i] = NULL;

    this->numMaxSlots	= numMaxSlots;
    this->numElements 	= 0;
}

template <class T>
void hcSortedList<T>::expand(uint numNewElements)
{
#ifdef TODO
    printf("INFO: hcSortedList::expand: change set size from %u to %u.\n", numMaxSlots, numMaxSlots + numNewElements);
    printf("INFO: If you see this message too often, your program is probably in need of optimization\n");
    printf("INFO: meaning that this procedure is inefficient)\n");
#endif

    if(numNewElements==0)
    	numNewElements = 1024;

    T **tempElementStorage  = new T*[numMaxSlots + numNewElements];

    for(uint i=0; i<numMaxSlots + numNewElements; ++i)
    	tempElementStorage[i] = NULL;

    for(uint i=0;i<numElements;++i)
        tempElementStorage[i]   = elements[i];

    numMaxSlots     += numNewElements;

    delete [] elements;

    elements    = tempElementStorage;
}

template <class T>
int hcSortedList<T>::findElement(const T &object)
{
    for(uint i=0;i <numElements; ++i)
		if(elements[i] == &object)
			return i;

    return -1;
}

template <class T>
bool hcSortedList<T>::hasElement(const T &object)
{
    if(findElement(object) > -1)
        return true;

    return false;
}


/*!	\returns position where object is stored
 *
 */
template <class T>
uint hcSortedList<T>::insertElement(T &object)
{
    if(numElements == numMaxSlots)
    	expand(numMaxSlots);

    uint l = 0;
    uint m = numElements/2;
    uint h = numElements-1;

    while(true)
    {
    	if(numElements == 0)									    		break;
    	if(object <= *elements[m] && (m==0 || object >= *elements[m-1])) 	break;

    	if(object > *elements[m])
    	{
    		if(m==h-1 || m==h)
    		{
    			if(object > *elements[h])	m=h+1;
    			else	    				m=h;
    			break;
    		}
    		l = m;
    		m = (h-m)/2+l;
    	}
    	else
    	{
    		if(m==l) break;
    		h = m;
    		m = (m-l)/2+l;
    	}
    }

    T **tempElements = new T*[numMaxSlots];
    for(uint j=0; j<m; ++j)
		tempElements[j] = elements[j];
    tempElements[m] = &object;
    for(uint j=m; j<numElements; ++j)
		tempElements[j+1] = elements[j];

    delete [] elements;
    elements = tempElements;
    ++numElements;

    return m;
}

/*! \returns position where object has been found and removed or -1 on failure
 *
 */
template <class T>
int hcSortedList<T>::removeElement(T &object)
{
	int retval = -1;

    for(uint i=0;i<numElements;++i)
		if(&object == elements[i])
			retval = i;

    if(retval > -1)
    {
    	for(uint i=retval+1; i<numElements; ++i)
    		elements[i-1] = elements[i];

    	 elements[numElements-1] = NULL;
		--numElements;
    }

    return retval;
}

template <class T>
bool hcSortedList<T>::sort()
{
	bool switched 	= false;
	bool retval		= false;

	do
	{
		switched = false;
		for(uint i=1; i<numElements; ++i)
			if(*elements[i-1] > *elements[i])
			{
				T *temp 		= elements[i];
				elements[i]		= elements[i-1];
				elements[i-1]	= temp;
				switched 		= true;
				retval 			= true;
			}
	}while(switched);

	return retval;
}

template <class T>
void hcSortedList<T>::dump() const{

    printf("--- Dumping hcSortedList:\n");
    printf("# of slots: %u\n", this->numMaxSlots);
    printf("# of occupied slots: %u\n", this->numElements);
}

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//
//                          hcSortedListStorage<T>
//
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

template <class T>
class hcSortedListStorage : public hcSortedList<T>{
public:

	hcSortedListStorage(uint numMaxSlots = 128);			/*!< \brief std constructor										*/
	hcSortedListStorage(const hcSortedListStorage &other);	/*!< \brief cpy constructor										*/
    virtual ~hcSortedListStorage();

    hcSortedListStorage<T> &operator=(const hcSortedListStorage<T> &other);

    virtual int removeElement(T &object);       	/*!< \brief removes and destroys an element             */
    bool extractElement(T &object);              	/*!< \brief removes an element and transfers ownership  */
};

template <class T>
hcSortedListStorage<T>::hcSortedListStorage(uint numMaxSlots) : hcSortedList<T>(numMaxSlots){}

template <class T>
hcSortedListStorage<T>::hcSortedListStorage(const hcSortedListStorage<T> &other) :
    hcSortedList<T>(other.numMaxSlots)
{
    *this = other;
}


template <class T>
hcSortedListStorage<T>::~hcSortedListStorage()
{
    this->clear();
}

template <class T>
hcSortedListStorage<T> &hcSortedListStorage<T>::operator=(const hcSortedListStorage<T> &other){

	if(this == &other)
		return *this;

	hcSortedList<T>::clear();
	this->init(other.numMaxSlots);
	this->numElements = other.numElements;

	for(uint i=0;i<this->numElements;++i)
		this->elements[i] = new T(*other.elements[i]);

	return *this;
}

template <class T>
int hcSortedListStorage<T>::removeElement(T &object)
{
	int retval = -1;

	for(uint i=0;i<this->numElements;++i)
		if(&object == this->elements[i])
		{
        	//void *q = (void*)this->elements[i];
        	//uintptr_t s = (uintptr_t)q;
        	//printf("hcSortedListStorage::removeElement %llx\n", s);fflush(stdout);
			delete this->elements[i];
			retval = i;
			break;
		}

	if(retval > -1)
	{
		for(uint i=retval+1; i<this->numElements; ++i)
			this->elements[i-1] = this->elements[i];
	}

	this->elements[this->numElements-1] = NULL;
	--this->numElements;

	return retval;
}

template <class T>
bool hcSortedListStorage<T>::extractElement(T &object)
{
	int retval = -1;

	for(uint i=0;i<this->numElements;++i)
		if(&object == this->elements[i])
		{
			retval = i;
			break;
		}

	if(retval > -1)
	{
		for(uint i=retval+1; i<this->numElements; ++i)
			this->elements[i-1] = this->elements[i];

		this->elements[this->numElements-1] = NULL;
		--this->numElements;
	}

	return retval > -1;
}

#endif // SORTEDLIST_H
