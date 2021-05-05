#-----------------------------------------------------------------------------
# PFSS
#-----------------------------------------------------------------------------

MAKEFILE      = Makefile
VERSION_MAJOR	= 0						# Major Version Number
VERSION_MINOR	= 1						# Minor Version Number

####### Compiler, tools and options

CUDA			= 1
DEFINES       	=  -D_DEBUG
DEFINES			+= -DhcFloat=float		# floating point precision to be used
DEFINES			+= -DNUMTHREADS=12		# number of thread to be used in multi-threaded computations
DEFINES			+= -DSPRO=1 -DMPRO=2 -DMCUDA=3 -DMAGMAPPING=MPRO # magnetic mapping to be computed by single thread (SPRO), multiple threads (MPRO) or CUDA (do not use that last one, its not efficient)
#DEFINES		+= -DTODO				# shows hints where improvements are due
#DEFINES			+= -DPRESENTATION		# better visibility for presentation
DEFINES			+= -DRSCALE				# scales radial coordinates by a factor of r_sol for numeric stability
DEFINES			+= -DREMOVEMONOPOLE		# removes artificial magnetic monopole from magnetograms (WSO does not remove monopole from their data before SHC computation, see their coefficient files)
DEFINES			+= -DVERBOSE			# gives more output
DEFINES			+=-DPFSS_VER_MAJ=${VERSION_MAJOR}
DEFINES			+=-DPFSS_VER_MIN=${VERSION_MINOR}


ifeq ($(CUDA), 1)
	DEFINES		+= -DCUDA
	CUDAARCH	=-arch=sm_30
	ifeq ($(CUDAARCH), -arch=sm_35)
	DEFINES		+= -DCUDACC=35
	else ifeq ($(CUDAARCH), -arch=sm_30)
	DEFINES		+= -DCUDACC=30
	else ifeq ($(CUDAARCH), -arch=sm_21)
	DEFINES		+= -DCUDACC=21
	else ifeq ($(CUDAARCH), -arch=sm_20)
	DEFINES		+= -DCUDACC=20
	endif
endif

CC            	= gcc
CXX           	= g++
CFLAGS        	= -pipe -g -fPIC $(DEFINES)
CXXFLAGS      	= -pipe -g -std=c++0x -fPIC -O0 $(DEFINES) 
INCPATH       	= -I. $(CUDAINCPATH)
INCTEST			= -Itest
LINK          	= g++
LFLAGS        	= -O0	
LIBS			+= -lpthread 											# Multiprocessing library
LIBS			+= -lcfitsio 											# FITS library
LIBS			+= -lgsl -lgslcblas										# Gnu Science Library 
LIBS			+= -lboost_system -lboost_filesystem -lboost_regex		# Boost libraries
LIBS			+= -lfftw3												# Fast Fourier Transform library
LIBS			+= -lfreeimage											# Image processing library
LIBS			+= -l:cspice.a -lm										# NAIF spice library
DEL_FILE      	= rm -f
DEL_DIR       	= rmdir

ifeq ($(CUDA), 1)
	NVCC		= nvcc
	LIBS		+= -lcuda -lcudart
	NVCCFLAGS	= -g -G -O0 --ptxas-options=-v --maxrregcount=32 --machine 64 -x cu -rdc=true -Xcompiler -fPIC $(CUDAARCH) $(DEFINES)
endif 

####### Output directory

DST_DIR   = bin/

####### Files

OBJECTS	= main.o \
			hcConstants.o \
			hcTools.o \
			hcImage.o \
			carRotInfo.o \
			hcTime.o \
			hcFunction.o \
			pfssSolution.o \
			pfssSolutionInfo.o \
			pfssSolution_batches.o \
			laplaceSolver.o \
			hcImageFITS.o \
			synPhotMagfield.o \
			pfss.o \
			filenames.o \
			imageStatistics.o \
			hcLine.o \
			hcPlane.o
			
				
ifeq ($(CUDA), 1)
CUDA_OBJECTS 	= cuda_interface_cuda.o \
					cuda_grids_cuda.o \
					cuda_laplaceSolver_cuda.o \
					hcVec.o \
					hcMatrix.o \
					magMapping.o \
					magline.o \
					grids.o \
					ellipticalGrid.o
					
else
OBJECTS	+= hcVec.o \
			hcMatrix.o \
			magMapping.o \
			magline.o \
			grids.o \
			ellipticalGrid.o
endif

OBJECTS_IN_DIR = $(addprefix $(DST_DIR), $(OBJECTS))
CUDAOBJ_IN_DIR = $(addprefix $(DST_DIR), $(CUDA_OBJECTS))

TARGET        = pfss


first: all

####### Build rules

all: $(TARGET)

$(TARGET):  $(OBJECTS_IN_DIR) $(CUDAOBJ_IN_DIR)	
ifeq ($(CUDA), 1)
	$(NVCC) $(CUDAARCH) -Xcompiler -fPIC -dlink $(CUDAOBJ_IN_DIR) -o $(DST_DIR)cudalink.o
	$(LINK) $(LFLAGS) -o $(DST_DIR)$(TARGET) $(OBJECTS_IN_DIR) $(CUDAOBJ_IN_DIR) $(DST_DIR)cudalink.o $(LIBS)
else
	$(LINK) $(LFLAGS) -o $(DST_DIR)$(TARGET) $(OBJECTS_IN_DIR) $(LIBS)
endif

clean:
	-$(DEL_FILE) $(OBJECTS_IN_DIR) $(CUDAOBJ_IN_DIR) $(DST_DIR)cudalink.o	
	-$(DEL_FILE) *~ core *.core
	-$(DEL_FILE) $(DST_DIR)$(TARGET)
	
$(DST_DIR)main.o: main.cpp
	$(CXX) $(CXXFLAGS) $(INCPATH) -c -o $(DST_DIR)main.o main.cpp
	
$(DST_DIR)hcConstants.o: engine/hcConstants.cpp engine/hcConstants.h
	$(CXX) $(CXXFLAGS) $(INCPATH) -c -o $(DST_DIR)hcConstants.o engine/hcConstants.cpp
	
$(DST_DIR)coordTransform.o: engine/math/coordTransform.cpp engine/math/coordTransform.h
	$(CXX) $(CXXFLAGS) $(INCPATH) -c -o $(DST_DIR)coordTransform.o engine/math/coordTransform.cpp
	
$(DST_DIR)cuda_interface_cuda.o: src/cuda_interface.cu
	$(NVCC) $(NVCCFLAGS) $(INCPATH) -c -o $(DST_DIR)cuda_interface_cuda.o src/cuda_interface.cu

$(DST_DIR)cuda_grids_cuda.o: src/grids.h \
		src/cuda_interface.h \
		src/cuda_interface.cu \
		src/cuda_grids.cu
	$(NVCC) $(NVCCFLAGS) $(INCPATH) -c -o $(DST_DIR)cuda_grids_cuda.o src/cuda_grids.cu	 
	
$(DST_DIR)cuda_laplaceSolver_cuda.o: src/cuda_interface.h \
		src/cuda_interface.cu \
		src/cuda_grids.cu \
		src/cuda_laplaceSolver.cu
	$(NVCC) $(NVCCFLAGS) $(INCPATH) -c -o $(DST_DIR)cuda_laplaceSolver_cuda.o src/cuda_laplaceSolver.cu
	
$(DST_DIR)cuda_magMapping_cuda.o: src/cuda_magMapping.cu
	$(NVCC) $(NVCCFLAGS) $(INCPATH) -c -o $(DST_DIR)cuda_magMapping_cuda.o src/cuda_magMapping.cu
	
$(DST_DIR)cuda_magline_cuda.o: src/cuda_magline.cu
	$(NVCC) $(NVCCFLAGS) $(INCPATH) -c -o $(DST_DIR)cuda_magline_cuda.o src/cuda_magline.cu
	
$(DST_DIR)cuda_testEnv_cuda.o: engine/cuda_testEnv.cu	
	$(NVCC) $(NVCCFLAGS) $(INCPATH) -c -o $(DST_DIR)cuda_testEnv_cuda.o engine/cuda_testEnv.cu
	
$(DST_DIR)hcVec.o: engine/math/hcVec.cpp engine/math/hcVec.h
ifeq ($(CUDA), 1)
	$(NVCC) $(NVCCFLAGS) $(INCPATH) -c -o $(DST_DIR)hcVec.o engine/math/hcVec.cpp
else
	$(CXX) $(CXXFLAGS) $(INCPATH) -c -o $(DST_DIR)hcVec.o engine/math/hcVec.cpp
endif

$(DST_DIR)hcMatrix.o: engine/math/hcMatrix.cpp engine/math/hcMatrix.h
ifeq ($(CUDA), 1)
	$(NVCC) $(NVCCFLAGS) $(INCPATH) -c -o $(DST_DIR)hcMatrix.o engine/math/hcMatrix.cpp
else
	$(CXX) $(CXXFLAGS) $(INCPATH) -c -o $(DST_DIR)hcMatrix.o engine/math/hcMatrix.cpp
endif	

$(DST_DIR)grids.o: src/grids.cpp src/grids.h
ifeq ($(CUDA), 1)
	$(NVCC) $(NVCCFLAGS) $(INCPATH) -c -o $(DST_DIR)grids.o src/grids.cpp
else
	$(CXX) $(CXXFLAGS) $(INCPATH) -c -o $(DST_DIR)grids.o src/grids.cpp
endif

$(DST_DIR)ellipticalGrid.o: src/ellipticalGrid.cpp src/ellipticalGrid.h src/grids.h
ifeq ($(CUDA), 1)
	$(NVCC) $(NVCCFLAGS) $(INCPATH) -c -o $(DST_DIR)ellipticalGrid.o src/ellipticalGrid.cpp
else
	$(CXX) $(CXXFLAGS) $(INCPATH) -c -o $(DST_DIR)ellipticalGrid.o src/ellipticalGrid.cpp
endif

$(DST_DIR)pfss.o: src/pfss.cpp src/pfss.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o $(DST_DIR)pfss.o src/pfss.cpp
	
$(DST_DIR)magMapping.o: src/magMapping.cpp src/magMapping.h \
		src/magline.h \
		engine/hcImage.h \
		src/laplaceSolver.h
ifeq ($(CUDA), 1)
	$(NVCC) -c $(NVCCFLAGS) $(INCPATH) -o $(DST_DIR)magMapping.o src/magMapping.cpp
else
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o $(DST_DIR)magMapping.o src/magMapping.cpp
endif

$(DST_DIR)magline.o: src/magline.cpp src/magline.h
ifeq ($(CUDA), 1)
	$(NVCC) -c $(NVCCFLAGS) $(INCPATH) -o $(DST_DIR)magline.o src/magline.cpp
else
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o $(DST_DIR)magline.o src/magline.cpp
endif

$(DST_DIR)hcPlane.o: engine/math/hcPlane.cpp engine/math/hcPlane.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o $(DST_DIR)hcPlane.o engine/math/hcPlane.cpp
	
$(DST_DIR)hcImage.o: engine/hcImage.cpp engine/hcImage.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o $(DST_DIR)hcImage.o engine/hcImage.cpp
	
$(DST_DIR)hcImageFITS.o: engine/hcImageFITS.cpp engine/hcImageFITS.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o $(DST_DIR)hcImageFITS.o engine/hcImageFITS.cpp
	
$(DST_DIR)carRotInfo.o: src/carRotInfo.cpp src/carRotInfo.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o $(DST_DIR)carRotInfo.o src/carRotInfo.cpp
	
$(DST_DIR)synPhotMagfield.o: src/synPhotMagfield.cpp src/synPhotMagfield.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o $(DST_DIR)synPhotMagfield.o src/synPhotMagfield.cpp
	
$(DST_DIR)hcTools.o: engine/hcTools.cpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o $(DST_DIR)hcTools.o engine/hcTools.cpp	
	
$(DST_DIR)hcTime.o: engine/hcTime.cpp engine/hcTime.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o $(DST_DIR)hcTime.o engine/hcTime.cpp
	
$(DST_DIR)hcFunction.o: engine/math/hcFunction.cpp engine/math/hcFunction.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o $(DST_DIR)hcFunction.o engine/math/hcFunction.cpp
	
$(DST_DIR)pfssSolution.o: src/pfssSolution.cpp src/pfssSolution.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o $(DST_DIR)pfssSolution.o src/pfssSolution.cpp
	
$(DST_DIR)pfssSolutionInfo.o:  src/pfssSolutionInfo.cpp
	$(CXX) $(CXXFLAGS) $(INCPATH) -c -o $(DST_DIR)pfssSolutionInfo.o src/pfssSolutionInfo.cpp
	
$(DST_DIR)pfssSolution_batches.o:  src/pfssSolution_batches.cpp
	$(CXX) $(CXXFLAGS) $(INCPATH) -c -o $(DST_DIR)pfssSolution_batches.o src/pfssSolution_batches.cpp
	
$(DST_DIR)laplaceSolver.o: src/laplaceSolver.cpp src/laplaceSolver.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o $(DST_DIR)laplaceSolver.o src/laplaceSolver.cpp
	
$(DST_DIR)filenames.o:  src/filenames.cpp
	$(CXX) $(CXXFLAGS) $(INCPATH) -c -o $(DST_DIR)filenames.o src/filenames.cpp
	
$(DST_DIR)imageStatistics.o: src/imageStatistics.cpp src/imageStatistics.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o $(DST_DIR)imageStatistics.o src/imageStatistics.cpp
	
$(DST_DIR)hcLine.o: engine/math/hcLine.cpp engine/math/hcLine.h
	$(CXX) $(CXXFLAGS) $(INCPATH) -c -o $(DST_DIR)hcLine.o engine/math/hcLine.cpp

