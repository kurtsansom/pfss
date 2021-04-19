#ifndef HCIMAGE_H
#define HCIMAGE_H

#ifdef GUI
#include "extern/soil/SOIL.h"
#endif

#include "engine/hcTools.h"
#include "engine/math/hcVec.h"

#include "FreeImage.h"
#include "stdio.h"
#include <cstring>

typedef unsigned int uint;



//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//
//                          hcImage
//
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

template <class T>
class hcImage{
protected:
	T *data;

public:

    uint width;
    uint height;

    hcImage();                          /*!< \brief std constructor         */
    hcImage(const hcImage<T> &other);   /*!< \brief cpy constructor         */
    hcImage(uint width, uint height);
    virtual ~hcImage();

    hcImage &operator=(const hcImage &other);	/*!< \brief assignment operator											*/

    T &operator()(uint x, uint y);		/*!< \brief returns manipulatable reference to contents							*/

    T content(uint x, uint y) const;	/*!< \brief returns copy of contents											*/

    void initNULL();
    void clear();
    void init(uint width, uint height);
        /*!< \brief allocates memory for an image of the according dimenstions                                  		*/

    //TODO: following not working? return value that will be destroyed at end of function?
    hcImage extractSubimage(int extractPosX, int extractPosY, uint width, uint height);
        /*!< \brief extracts subimage of (width * height) from the specified position                           		*/

    bool extractSubimage(int extractPosX, int extractPosY, hcImage<T> &retval);
        /*!< \brief see above                                                                                   		*/

    void shiftXdirection(int numPixels);
    	/*!< \brief shifts image numPixels in positive x-direction, pixels leaving right boundary enter from the left	*/

    virtual bool load(const string &filename){return true;}
    	/*!< \brief loads image from file																				*/

    virtual bool save(const string &filename){return true;}
    	/*!< \brief saves image to file																					*/

    void pushToArray(T **array, uint &width, uint &height);
        /*!< \brief saves data to array and reports back width and height                                       		*/

    T *getData(){return data;}	// TODO get rid of this

    void dump() const;
};

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//
//                          hcImageRGBA
//
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

class hcImageRGBA : public hcImage<uint>{
public:

    hcImageRGBA();                          /*!< \brief std constructor         */
    hcImageRGBA(const hcImageRGBA &other);  /*!< \brief cpy constructor         */
    hcImageRGBA(uint width, uint height);
    hcImageRGBA(uint width, uint height, uint color);
    virtual ~hcImageRGBA();

    hcImageRGBA &operator=(const hcImageRGBA &other);

    void init(uint width, uint height);

    void init(uint width, uint height, uint color);
        /*!< \brief allocates memory for an image of the specified dimensions and paints it according to color  */

    virtual bool save(const string &filename);
        /*!< \brief saves to filename                                                                           */

#ifdef GUI
    virtual bool load(const string &filename);
        /*!< \brief loads an image from filename                                                                */
#endif

    uint numNotBlackPixels() const;

    void magDiff(hcImage &img1, hcImage &img2);

    void setBackgroundColor(uint bgColor, int insertPosX=0, int insertPosY=0, uint subWidth=0, uint subHeight=0);
        /*!< \brief paints the entire image or parts of it in bgColor                                                          */

    bool insertSubimage(hcImage &src, int insertPosX, int insertPosY);
        /*!< \brief inserts a image into this one at the specified location                                     */

    bool loadFromArray(uint *array, uint width, uint height);
        /*!< \brief loads image data from array                                                                 */

    bool interpolateRectangularImage(uint width_out, uint height_out);
        /*!< \brief interpolates image to the requested dimension                                               */

    void dump() const;
};

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//
//                          hcImageFloat
//
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

class hcImageFloat : public hcImage<hcFloat>{
public:

    hcImageFloat();                             /*!< \brief std constructor         */
    hcImageFloat(const hcImageFloat &other);    /*!< \brief cpy constructor         */
    hcImageFloat(uint width, uint height);
    hcImageFloat(uint width, uint height, float color);
    virtual ~hcImageFloat();

    hcImageFloat &operator=(const hcImageFloat &other);

    void init(uint width, uint height);

    void init(uint width, uint height, float color);
        /*!< \brief allocates memory for an image of the specified dimensions and paints it according to color  				*/

    virtual bool load(const string &filename);
		/*!< \brief loads an image from filename                                                                				*/

    virtual bool save(const string &filename);
        /*!< \brief saves to filename                                                                           				*/

    void reldiffImage(const hcImageFloat &img0, const hcImageFloat &img1);
    	/*!< \brief computes pixelwise relative difference between img0 and img1												*/

    void diffImage(const hcImageFloat &img0, const hcImageFloat &img1);
    	/*!< \brief compute pixelwise difference between img0 and img1															*/

    float absMean();
    	/*!< \brief sums up all absolute pixel values and divides by number of pixels											*/

    void setBackgroundColor(float bgColor, int insertPosX=0, int insertPosY=0, uint subWidth=0, uint subHeight=0);
        /*!< \brief paints the entire image or parts of it in bgColor                                                          	*/

    bool insertSubimage(const hcImageFloat &src, int insertPosX, int insertPosY);
        /*!< \brief inserts a image into this one at the specified location                                     				*/

    bool loadFromArray(float *array, uint width, uint height);
        /*!< \brief loads image data from array                                                                 				*/

    virtual void addHomWhiteNoise(float sigma, uint seed);
    	/*!< \brief adds white noise by adding gaussian distributed value														*/

    virtual void addPixelNoise(float fraction, uint seed);
    	/*!< \brief adds noise to each pixel seperately by adding normal distributed values with sigma = fraction * pixelValue	*/

    virtual void addSignedFractionNoise(float fraction, uint seed);
    	/*!< \brief adds +/- fraction * pixel value to each pixel																*/

    hcFloat meanSquaredDiff(hcImage &other);

    void meanFilter(uint windowWidth, bool circular);
    	/*!< \brief applies mean filter to picture, x-coordinate may be circular												*/

    void medianFilter(uint windowWidth, bool circular);
    	/*!< \brief applies median filter to picture, x-coordinate may be circular												*/

    //void rescale(uint newWidth, uint newHeight);
    	/*!< \brief rescales image to new dimensions																			*/

    void dump() const;
};


//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//
//                          hcImageVec3D
//
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

class hcImageVec3D : public hcImage<Vec3D>{
public:

	hcImageVec3D();                          /*!< \brief std constructor         */
	hcImageVec3D(const hcImageVec3D &other);  /*!< \brief cpy constructor         */
	hcImageVec3D(uint width, uint height);
    virtual ~hcImageVec3D();

    hcImageVec3D &operator=(const hcImageVec3D &other);

    void init(uint width, uint height);

    void dump() const;
};

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//
//                          hcImageInt
//
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

class hcImageInt : public hcImage<int>{
public:

	hcImageInt();                          /*!< \brief std constructor         */
	hcImageInt(const hcImageInt &other);  /*!< \brief cpy constructor         */
	hcImageInt(uint width, uint height);
    virtual ~hcImageInt();

    hcImageInt &operator=(const hcImageInt &other);

    void init(uint width, uint height);

    void dump() const;
};

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//
//                          hcImageBool
//
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

class hcImageBool : public hcImage<bool>{
public:

	hcImageBool();                          /*!< \brief std constructor         */
	hcImageBool(const hcImageBool &other);  /*!< \brief cpy constructor         */
	hcImageBool(uint width, uint height);
    virtual ~hcImageBool();

    hcImageBool &operator=(const hcImageBool &other);

    void init(uint width, uint height);

    void dump() const;
};

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//
//                          hcImage
//
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

template <class T>
hcImage<T>::hcImage()
{
    initNULL();
}

template <class T>
hcImage<T>::hcImage(const hcImage<T> &other)
{
    initNULL();
    *this = other;
}

template <class T>
hcImage<T>::hcImage(uint width, uint height)
{
    initNULL();
    init(width, height);
}

template <class T>
hcImage<T>::~hcImage()
{
    clear();
}

template <class T>
hcImage<T> &hcImage<T>::operator=(const hcImage<T> &other)
{
    if(this == &other)
        return *this;

    init(other.width, other.height);

    memcpy((void*)this->data, (void*)(other.data), sizeof(T) * width * height);

    return *this;
}

template <class T>
T &hcImage<T>::operator()(uint x, uint y)
{
    if(x >= width || y >= height)
    {
        cerr << __FILE__ << ":" << __LINE__ << ": Indices out of bounds (x = " << x << ", y = " << y;
        cerr << ", width = " << width << ", height: " << height << ")\n";
        return data[0];
    }

    return data[y * width + x];
}

template <class T>
T hcImage<T>::content(uint x, uint y) const
{
    if(x >= width || y >= height)
    {
    	cerr << __FILE__ << ":" << __LINE__ << ": Indices out of bounds (x = " << x << ", y = " << y;
		cerr << ", width = " << width << ", height: " << height << ")\n";
		return data[0];
    }

    return data[y * width + x];
}

template <class T>
void hcImage<T>::initNULL()
{
    data    = NULL;
    width   = 0;
    height  = 0;
}

template <class T>
void hcImage<T>::clear()
{
    delete [] data;
    initNULL();
}

template <class T>
void hcImage<T>::init(uint width, uint height)
{
    clear();
    this->width     = width;
    this->height    = height;
    uint numElements= width*height;

    if(numElements > 0)				 	data 	= new T[numElements];
    for(uint i=0; i<numElements; ++i)  	data[i]	= T();
}

template <class T>
hcImage<T> hcImage<T>::extractSubimage(int extractPosX, int extractPosY, uint width, uint height)
{
    hcImage<T> retval;
    retval.init(width, height);
    for(int x = extractPosX; x < extractPosX+width; ++x)
        for(int y = extractPosY; y < extractPosY+height; ++y)
            if(!((x < 0) || (x >= width) || (y < 0) || (y >= height)))
                retval(x-extractPosX, y-extractPosY) = this->operator()(x, y);

    return retval;
}

template <class T>
bool hcImage<T>::extractSubimage(int extractPosX, int extractPosY, hcImage<T> &retval)
{
    for(int x = extractPosX; x < extractPosX+retval.width; ++x)
        for(int y = extractPosY; y < extractPosY+retval.height; ++y)
            if(!((x < 0) || (x >= width) || (y < 0) || (y >= height)))
                retval(x-extractPosX, y-extractPosY) = this->operator()(x, y);

    return true;
}

template <class T>
void hcImage<T>::shiftXdirection(int numPixels)
{
	if(abs(numPixels) > width)
	{
		int sign 	= (numPixels < 0 ? -1 : 1);
		uint n		= abs(numPixels) / width;
		numPixels -=  sign * n * width;
	}

	T *newData = new T[width*height];

	for(uint x=0; x<width; ++x)
		for(uint y=0; y<height; ++y)
		{
			int oldX = (int)x - numPixels;
			oldX += (oldX < 0 				? width : 0);
			oldX -= (oldX > (int)width-1	? width : 0);

			newData[y * width + x] = data[y * width + oldX];
		}

	delete [] data;
	data = newData;
}

template <class T>
void hcImage<T>::pushToArray(T **array, uint &width, uint &height)
{
    width   = this->width;
    height  = this->height;
    *array  = new T[width * height];

    for(uint x=0; x<width; ++x)
        for(uint y=0; y<height; ++y)
            (*array)[y * width + x] = data[y * width + x];
}

template <class T>
void hcImage<T>::dump() const
{
    printf("--- Dumping hcImage<T>:\n");
    printf("Widht:\t%u\n", width);
    printf("Height:\t%u\n", height);
}

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//
//                          hcScribble
//
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

/*!< \brief implements objects (lines, points, crosses, whatever) for scribbeling in hcImages
 */
class hcScribble : public hcImageRGBA{
public:

    hcImageRGBA *canvas;        /*!< \brief not-owning pointer to the image to be scribbled in                      */
    hcImageRGBA canvasOrig;     /*!< \brief canvas content that has been overdrawn by this scribble (TODO: this allows scribbles to think other scribbles are the original image) */
    hcImageRGBA original;       /*!< \brief not scaled image to be scribbled in canvas                              */

    int insertPosX;             /*!< \brief x-pos of this in canvas                                                 */
    int insertPosY;             /*!< \brief y-pos of this in canvas                                                 */

    double relPosX;             /*!< \brief relative (-1 < x < 1) x-pos of this in canvas                           */
    double relPosY;             /*!< \brief relative (-1 < x < 1) y-pos of this in canvas                           */

    double newRelPosX;          /*!< \brief where to draw when draw-function is called next time                    */
    double newRelPosY;          /*!< \brief where to draw when draw-function is called next time                    */

    double relWidth;            /*!< \brief relative width of this scribble in canvas                               */

    hcScribble();
    hcScribble(const hcScribble &other);
    virtual ~hcScribble();

    hcScribble &operator=(const hcScribble &other);

    void initNULL();
    virtual void clear();
    void init(double relWidth, double relPosX, double relPosY, hcImageRGBA *canvas);

    void setCanvas(hcImageRGBA *canvas);
    	/*!< \brief sets canvas to be painted upon																	*/

    void reinit();
        /*!< \brief this function empties canvasOrig so that a new canvas can be loaded without parts of the old one be dumped in it */

    bool updateScaledScribble();
        /*!< \brief adjusts this hcImageRGBA according to relative position / width                                 */

    bool updatePosition();
        /*!< \brief recomputes the insertPos of the scribble according to relPos                                    */

    bool draw();
        /*!< \brief draws this scribble at the specified position                                                   */

    bool redrawOrig();
        /*!< \brief redraws the original part of the image where the scribble was drawn                             */

    void dump() const;
};

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//
//                          hcScribbleVertLine
//
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------


class hcScribbleVertLine : public hcScribble{
public:

    hcScribbleVertLine();                                   /*!< \brief std constructor                             */
    hcScribbleVertLine(const hcScribbleVertLine &other);    /*!< \brief cpy constructor                             */
    hcScribbleVertLine(uint color, double relWidth, double relPosX, double relPosY, hcImageRGBA *canvas);
    virtual ~hcScribbleVertLine();

    hcScribbleVertLine &operator=(const hcScribbleVertLine &other);

    void initNULL();
    virtual void clear();
    void init(uint color, double relWidth, double relPosX, double relPosY, hcImageRGBA *canvas);

    void dump() const;
};

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//
//                          hcScribbleDot
//
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------


class hcScribbleDot : public hcScribble{
public:

    hcScribbleDot();                                    /*!< \brief std constructor                             */
    hcScribbleDot(const hcScribbleDot &other);          /*!< \brief cpy constructor                             */
    hcScribbleDot(uint color, double relWidth, double relPosX, double relPosY, hcImageRGBA *canvas);
    virtual ~hcScribbleDot();

    hcScribbleDot &operator=(const hcScribbleDot &other);

    void initNULL();
    virtual void clear();
    void init(uint color, double relWidth, double relPosX, double relPosY, hcImageRGBA *canvas);

    void dump() const;
};

#endif

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//
//                          other functions
//
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------


hcFloat getMaxSinLat(uint numPixels, hcFloat northPixelBorder);
	/*!< \brief returns the highest possible value for north-most pixel center so that its border will not exceed northPixelBorder 	*/

void bordertest(uint numPixels, hcFloat minSinLat, hcFloat maxSinLat);
	/*!< \brief test procedure to see metrics for sinLatitude magnetograms															*/

void setBackgroundColor(hcImageRGBA &img, uint bgColor);
    /*!< \brief fills image with bgColor (in RGBA8-format) */

uint insertSubimage(
        uint *dest, uint destWidth, uint destHeight,
        uint *src,  uint srcWidth,  uint srcHeight,
        uint insertPosX, uint insertPosY);
    /*!< \brief inserts a smaller image (src) into a bigger one (dst) at position (insertPosX, insertPosY) */

void interpolateRectangularImage_orig(uint *image_in, uint *image_out, uint width_in, uint height_in, uint width_out, uint height_out);

void createLogPlot(float *image_in, uint width, uint height, uint *image_out);
void createOneSidedLogPlot(float *image_in, uint width, uint height, uint *image_out);
