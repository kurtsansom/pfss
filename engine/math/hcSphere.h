#ifndef SPHERES_H
#define SPHERES_H
#include "engine/math/hcCircle.h"
#include "engine/math/hcVec.h"

/*! \brief mathematical model and functions for nD-spheres
 */
template<uint dim>
class hcSphere : public hcCircle<dim>{
public:

#ifdef __NVCC__
      __host__ __device__
#endif
      hcSphere();												/*!< \brief std constructor					*/

#ifdef __NVCC__
      __host__ __device__
#endif
      hcSphere(const hcSphere<dim>& sphere);					/*!< \brief cpy constructor					*/

#ifdef __NVCC__
      __host__ __device__
#endif
      hcSphere(const Vec<dim, hcFloat> &pos, float radius);		/*!< \brief constructor						*/

#ifdef __NVCC__
      __host__ __device__
#endif
      ~hcSphere();												/*!< \brief destructor						*/

#ifdef __NVCC__
      __host__ __device__
#endif
      hcSphere<dim> &operator=(const hcSphere<dim> &other);		/*!< \brief assignment operator				*/

#ifdef __NVCC__
      __host__ __device__
#endif
      void init(const Vec<dim, hcFloat> &pos, float radius);

#ifdef __NVCC__
      __host__ __device__
#endif
      void dump() const;
};

template<uint dim>
hcSphere<dim>::hcSphere(){}

template<uint dim>
hcSphere<dim>::hcSphere(const hcSphere<dim> &sphere) : hcCircle<dim>(sphere){

}

template<uint dim>
hcSphere<dim>::hcSphere(const Vec<dim, hcFloat> &pos, float radius){

    init(pos, radius);
}

template<uint dim>
hcSphere<dim>::~hcSphere(){}

template<uint dim>
hcSphere<dim> &hcSphere<dim>::operator=(const hcSphere<dim> &other){

	if(this == &other)
		return *this;

	hcCircle<dim>::operator=(other);

	return *this;
}

template<uint dim>
void hcSphere<dim>::init(const Vec<dim, hcFloat> &pos, float radius){

    hcCircle<dim>::init(pos, radius);
}

template<uint dim>
void hcSphere<dim>::dump() const{

    printf("Dumping hcSphere:\n");
    hcCircle<dim>::dump();
}

#endif
