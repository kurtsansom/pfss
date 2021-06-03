#include "engine/hcImage.h"

#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_cdf.h"

/*! @param numPixels 			number of pixels in meridional direction
 *  @param northPixelBorder		border of north-most pixel (in sin(latitude)), 1.0 for pixels stretching all the way to the north pole
 */
hcFloat getMaxSinLat(uint numPixels, hcFloat northPixelBorder)
{
	return (northPixelBorder * (numPixels-1) / numPixels);
}

/*! gives information on pixels for different coordinate schemes
 *
 * 	case 1 - Hoeksema
 *
 * 				image stretches in meridional direction from maxSinLat to minSinLat, where
 * 				the center coordinates of the outermost pixels are maxSinLat and minSinLat,
 * 				the border of the pixels lie even further out
 *
 * 				------
 * 				|    |
 * 				|    |   --- maxSinLat
 * 				|    |
 * 				------
 *
 * 	case 2 - mine
 *
 * 				------   --- maxSinLat
 * 				|    |
 * 				|    |
 * 				|    |
 * 				------
 */
void bordertest(uint numPixels, hcFloat minSinLat, hcFloat maxSinLat)
{
	hcFloat dSinLat1	= (maxSinLat - minSinLat) / (numPixels-1);
	hcFloat dSinLat2	= (maxSinLat - minSinLat) / numPixels;

	hcFloat dist1[numPixels];	// distances between pixel centers in rad (case 1)
	hcFloat dist2[numPixels];	// distances between pixel centers in rad (case 2)
	hcFloat pos1[numPixels];	// center positions of pixels in rad (case 1)
	hcFloat pos2[numPixels];	// center positions of pixels in rad (case 2)
	hcFloat size1[numPixels];	// pixel size in rad (case 1)
	hcFloat size2[numPixels];	// pixel size in rad (case 2)
	hcFloat res1[numPixels];	// resolution in numPixels/rad (case 1)
	hcFloat res2[numPixels];	// resolution in numPixels/rad (case 2)

	printf("\nnumPixels:\t%u\n", numPixels);
	printf("maxSinLat:\t%E\n", maxSinLat);
	printf("minSinLat:\t%E\n", minSinLat);

	printf("\nAngular center position of pixels:\n");
	printf("pixel#\ttheta/rad(c1)\ttheta/rad(c2)\n");
	for(uint i=0; i<numPixels; ++i)
	{
		pos1[i]		= PI/2.0 - asin(maxSinLat - dSinLat1 * i);
		pos2[i]		= PI/2.0 - asin(maxSinLat - dSinLat2 * (i+1.0/2.0));

		size1[i]	= asin(maxSinLat - dSinLat1 * (i - 1.0/2.0)) - asin(maxSinLat - dSinLat1 * (i + 1.0/2.0));
		size2[i]	= asin(maxSinLat - dSinLat2 * i)       		 - asin(maxSinLat - dSinLat2 * (i + 1.0));
		printf("%u\t%E\t%E\n", i, pos1[i], pos2[i]);
	}

	printf("\nAngular size of pixels:\n");
	printf("pixel#\tsize/rad(c1)\tsize/rad(c2)\n");
	for(uint i=0; i<numPixels; ++i)
	{
		size1[i]	= asin(maxSinLat - dSinLat1 * (i - 1.0/2.0)) - asin(maxSinLat - dSinLat1 * (i + 1.0/2.0));
		size2[i]	= asin(maxSinLat - dSinLat2 * i)       		 - asin(maxSinLat - dSinLat2 * (i + 1.0));
		printf("%u\t%E\t%E\n", i, size1[i], size2[i]);
	}

	printf("\nAngular resolution in numPixels/rad:\n");
	printf("pixel#\tcase 1\tcase 2\n");
	for(uint i=0; i<numPixels; ++i)
	{
		printf("%u\t%E\t%E\n", i, 1/size1[i], 1/size2[i]);
	}

	printf("\nAngular distance to previous pixels\n");
	printf("pixel#\tdist/rad(c1)\tdist/rad(c2)\n");
	for(uint i=0; i<numPixels; ++i)
	{
		if(i==0)
		{
			dist1[i] = 0.0;
			dist2[i] = 0.0;
		}
		else
		{
			dist1[i] = pos1[i] - pos1[i-1];
			dist2[i] = pos2[i] - pos2[i-1];
		}

		printf("%u\t%E\t%E\n", i, dist1[i], dist2[i]);
	}

	hcFloat lat 	= asin(maxSinLat+1.0/2.0*dSinLat1);
	hcFloat theta	= PI/2.0 - lat;
	printf("\nCase 1 north border:\nsin(lat)\tlat/rad\t\tlat/deg\t\ttheta/rad\ttheta/deg\n%E\t%E\t%E\t%E\t%E\n",
				maxSinLat+1.0/2.0*dSinLat1, lat, lat*rad2deg, theta, theta*rad2deg);

	lat 	= asin(maxSinLat);
	theta	= PI/2.0 - lat;
	printf("\nCase 1 north border:\nsin(lat)\tlat/rad\t\tlat/deg\t\ttheta/rad\ttheta/deg\n%E\t%E\t%E\t%E\t%E\n",
					maxSinLat, lat, lat*rad2deg, theta, theta*rad2deg);
}

void setBackgroundColor(hcImageRGBA &img, uint bgColor)
{
    uint x,y;
    for(y=0;y<img.height;++y)
        for(x=0;x<img.width;++x)
            img(x,y) = bgColor;
}

uint insertSubimage(uint *dest, uint destWidth, uint destHeight, uint *src, uint srcWidth, uint srcHeight, uint insertPosX, uint insertPosY){

    if(srcWidth + insertPosX > destWidth)
    {
        printf("ERROR! insertSubimage: Subimage does not fit in destination image!\n(destWidth: %u, srcWidth: %u, insertPosX: %u)\n", destWidth, srcWidth, insertPosX);
        return 0;
    }
    if(srcHeight + insertPosY > destHeight)
    {
        printf("ERROR! insertSubimage: Subimage does not fit in destination image!\n(destHeight: %u, srcHeight: %u, insertPosY: %u)\n", destHeight, srcHeight, insertPosY);
        return 0;
    }

    uint x, y;
    float alpha;
    uint colorDest, colorSrc;
    unsigned char rDst, gDst, bDst, aDst;
    unsigned char rSrc, gSrc, bSrc, aSrc;

    for(x=0;x<srcWidth;++x)
        for(y=0;y<srcHeight;++y)
        {
            colorDest 	= dest[(y + insertPosY) * destWidth + (x + insertPosX)];
            colorSrc 	= src[y * srcWidth + x];

            RGBA82char(colorSrc, rSrc, gSrc, bSrc, aSrc);
            RGBA82char(colorDest, rDst, gDst, bDst, aDst);

            alpha = aSrc / 255.0;
            uint finalAlpha = aDst > aSrc ? aDst : aSrc;

            dest[(y + insertPosY) * destWidth + (x + insertPosX)] = \
                    char2RGBA8( (1 - alpha) * rDst + alpha * rSrc,\
                                (1 - alpha) * gDst + alpha * gSrc,\
                                (1 - alpha) * bDst + alpha * bSrc,\
                                finalAlpha
                                );
        }
    return 1;
}

uint nearestPowerOfTwo(uint number){

    uint i = 0;
    while((unsigned)(1<<++i) < number);
    --i;
    if(fabs((signed)((1<<i) - (signed)number)) <= fabs((signed)((1<<(i+1)) - (signed)number)))
        return i;
    else
        return ++i;
}




/*! \brief resizes a rectangular image via bilinear interpolation
 *
 */
void interpolateRectangularImage_orig(uint *image_in, uint *image_out, uint width_in, uint height_in, uint width_out, uint height_out)
{
    uint i_out, j_out;
    uint *tempImage 		= NULL;
    uint *tempImage2 		= NULL;

    float scaleFactorWidth 	= (float)width_in / width_out;
    float scaleFactorHeight = (float)height_in / height_out;

    // if input image is twice as wide / tall, it has to be treated special to not
    // lose energy in image space (due to unsampled pixels in input image)
    if(scaleFactorWidth >= 2 || scaleFactorHeight >= 2)
    {

        uint newWidthIn 	= 1 << nearestPowerOfTwo(width_in);
        uint newHeightIn 	= 1 << nearestPowerOfTwo(height_in);

        tempImage = new uint[newWidthIn * newHeightIn];

        interpolateRectangularImage_orig(image_in, tempImage, width_in, height_in, newWidthIn, newHeightIn);

        scaleFactorWidth 	= (float)newWidthIn  / width_out;
        scaleFactorHeight 	= (float)newHeightIn / height_out;

        tempImage2 = new uint[newWidthIn * newHeightIn];

        while(scaleFactorWidth >= 2 || scaleFactorHeight >= 2)
        {
            unsigned char r1 = 0;
            unsigned char r2 = 0;
            unsigned char g1 = 0;
            unsigned char g2 = 0;
            unsigned char b1 = 0;
            unsigned char b2 = 0;
            unsigned char a1 = 0;
            unsigned char a2 = 0;

            bool signal = false;

            if(scaleFactorWidth >= 2)
            {
                signal = true;
                for(uint i=0;i<newWidthIn/2;++i)
                    for(uint j=0;j<newHeightIn;++j)
                    {
                        RGBA82char(tempImage[j * newWidthIn + 2*i],     r1, g1, b1, a1);
                        RGBA82char(tempImage[j * newWidthIn + 2*i + 1], r2, g2, b2, a2);
                        if(a2 > a1)
                        	a1 = a2;
                        tempImage2[j * newWidthIn/2 + i] = char2RGBA8(0.5*(r1+r2), 0.5*(g1+g2), 0.5*(b1+b2), a1);
                    }

                newWidthIn = newWidthIn / 2;
                scaleFactorWidth = (float) newWidthIn / width_out;
            }

            if(scaleFactorHeight >= 2)
            {
                if(signal)
                {
                    for(uint i=0;i<newWidthIn;++i)
                        for(uint j=0;j<newHeightIn;++j)
                            tempImage[j*newWidthIn+i] = tempImage2[j*newWidthIn+i];
                }
                for(uint i=0;i<newWidthIn;++i)
                    for(uint j=0;j<newHeightIn/2;++j)
                    {
                        RGBA82char(tempImage[(2*j) * newWidthIn + i],   r1, g1, b1, a1);
                        RGBA82char(tempImage[(2*j+1) * newWidthIn + i], r2, g2, b2, a2);
                        if(a2 > a1)
							a1 = a2;
                        tempImage2[j * newWidthIn + i] = char2RGBA8(0.5*(r1+r2), 0.5*(g1+g2), 0.5*(b1+b2), a1);
                    }

                newHeightIn = newHeightIn / 2;
                scaleFactorHeight = (float) newHeightIn / height_out;
            }

            // image still more than twice as wide/tall ? copy new image and repeat
            // : set new width_in / height_in, repoint image_in and let the
            // algorithm work its magic
            if (scaleFactorWidth >= 2 || scaleFactorHeight >= 2)
            {
                for(uint i=0;i<newWidthIn;++i)
                    for(uint j=0;j<newHeightIn;++j)
                         tempImage[j*newWidthIn+i] = tempImage2[j*newWidthIn+i];
            }
            else
            {
                image_in = tempImage2;
                width_in = newWidthIn;
                height_in = newHeightIn;
            }
        }
    }

    for(i_out = 0; i_out < width_out; ++i_out)
        for(j_out = 0; j_out < height_out; ++j_out)
        {
            int i  = floor((scaleFactorWidth   * (i_out + 0.5)) - 0.5);
            int j  = floor((scaleFactorHeight  * (j_out + 0.5)) - 0.5);
            if(i<0)
                i=0;
            if(j<0)
                j=0;

            float xOut  = (i_out    + 0.5) / width_out;
            float yOut  = (j_out    + 0.5) / height_out;
            float xi    = (i        + 0.5) / width_in;
            float xip   = (i + 1    + 0.5) / width_in;
            float yi    = (j        + 0.5) / height_in;
            float yip   = (j + 1    + 0.5) / height_in;

            unsigned char r_ij      = 0;
            unsigned char r_ipj     = 0;
            unsigned char r_ijp     = 0;
            unsigned char r_ipjp    = 0;

            unsigned char g_ij      = 0;
            unsigned char g_ipj     = 0;
            unsigned char g_ijp     = 0;
            unsigned char g_ipjp    = 0;

            unsigned char b_ij      = 0;
            unsigned char b_ipj     = 0;
            unsigned char b_ijp     = 0;
            unsigned char b_ipjp    = 0;

            unsigned char a_ij      = 0;
            unsigned char a_ipj     = 0;
            unsigned char a_ijp     = 0;
            unsigned char a_ipjp    = 0;


            RGBA82char(image_in[j   * width_in  + i],       r_ij, g_ij, b_ij, a_ij);
        if(i<(signed)width_in-1)
            RGBA82char(image_in[j   * width_in  + i + 1],   r_ipj, g_ipj, b_ipj, a_ipj);
        if(j<(signed)height_in-1)
            RGBA82char(image_in[(j+1)* width_in + i],       r_ijp, g_ijp, b_ijp, a_ijp);
        if(i<(signed)width_in-1 && j<(signed)height_in-1)
            RGBA82char(image_in[(j+1)* width_in + i + 1],   r_ipjp, g_ipjp, b_ipjp, a_ipjp);

            float a1 = (xip - xOut) / (xip - xi);
            float a2 = 1 - a1;

            float rIntL = a1 * r_ij + a2 * r_ipj;
            float gIntL = a1 * g_ij + a2 * g_ipj;
            float bIntL = a1 * b_ij + a2 * b_ipj;
            float aIntL = a_ij > a_ipj ? a_ij : a_ipj;
            //max(a_ij, a_ipj);

            float rIntU = a1 * r_ijp + a2 * r_ipjp;
            float gIntU = a1 * g_ijp + a2 * g_ipjp;
            float bIntU = a1 * b_ijp + a2 * b_ipjp;
            float aIntU = a_ijp > a_ipjp ? a_ijp : a_ipjp;
            //max(a_ijp, a_ipjp);

            a1 = (yip - yOut) / (yip - yi);
            a2 = 1 - a1;

            float red   = a1 * rIntL + a2 * rIntU;
            float green = a1 * gIntL + a2 * gIntU;
            float blue  = a1 * bIntL + a2 * bIntU;
            float alpha = a1 * aIntL + a2 * aIntU;

            image_out[j_out * width_out + i_out] = char2RGBA8(red, green, blue, alpha);
        }
    if(tempImage != NULL)
    {
        delete [] tempImage;
        delete [] tempImage2;
    }
}

void createLogPlot(float *image_in, uint width, uint height, uint *image_out)
{
    float maxValue = 4;
    float minValue = 2;

    for(uint x=0;x<width;++x)
        for(uint y=0;y<height;++y)
        {
            uint ind 	= y*width+x;
            float value = image_in[ind];
            int sign 	= value > 0 ? 1 : -1;
            value  		= (log10((fabs(value)*pow(10,minValue) + 1))) / (maxValue+minValue);
            if (value > 1.0)
                value = 1.0;

            char r,g,b;

            if(sign > 0)
                r = 192 * value + 64;
            else
                r = 64 * (1-value);

            if(sign < 0)
                b = 192 * value + 64;
            else
                b = 64 * (1-value);

            g = 64 * (1-value);

            image_out[ind] = char2RGBA8(r, g, b, 255);
        }
}

void createOneSidedLogPlot(float *image_in, uint width, uint height, uint *image_out)
{
    uint x,y;

    float maxValue = 2;
    float minValue = 4;

    for(x=0;x<width;++x)
        for(y=0;y<height;++y)
        {
            uint ind 	= y*width+x;
            float value = image_in[ind];
            int sign 	= value > 0 ? 1 : -1;
            value  		= (log10((fabs(value)*pow(10,minValue) + 1))) / (maxValue+minValue);
            if (value > 1.0)
                value = 1.0;

            char r,g,b;

            r = 255 * value;
            b = 0;
            g = 0;

            image_out[ind] = sign>0 ? char2RGBA8(r, g, b, 255) : char2RGBA8(0, 0, 0, 255);
        }
}

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//
//                          hcImageRGBA
//
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

hcImageRGBA::hcImageRGBA()
{
    initNULL();
}

hcImageRGBA::hcImageRGBA(const hcImageRGBA &other)
{
    initNULL();
    hcImage<uint>::operator=(other);
}

hcImageRGBA::hcImageRGBA(uint width, uint height)
{
    initNULL();
    init(width, height);
}

hcImageRGBA::hcImageRGBA(uint width, uint height, uint color)
{
    initNULL();
    init(width, height, color);
}

hcImageRGBA::~hcImageRGBA(){

    clear();
}

hcImageRGBA &hcImageRGBA::operator=(const hcImageRGBA &other)
{
    if(this == &other)
        return *this;

    hcImage<uint>::operator=(other);

    return *this;
}

void hcImageRGBA::init(uint width, uint height)
{
    hcImage<uint>::init(width, height);
}

void hcImageRGBA::init(uint width, uint height, uint color)
{
    clear();
    init(width, height);
    setBackgroundColor(color);
}

uint hcImageRGBA::numNotBlackPixels() const
{
	uint retval = 0;

	for(uint x=0; x<width; ++x)
		for(uint y=0; y<height; ++y)
		{
			unsigned char r, g, b, a;
			uint color = data[y*width+x];

			RGBA82char(color, r, g, b, a);

			if(r != 0 || g != 0 || b != 0)
				++retval;
		}

	return retval;
}

void hcImageRGBA::magDiff(hcImage &img1, hcImage &img2)
{

	if(img1.width != img2.width || img1.height != img2.height)
	{
		printf("ERROR! hcImageRGBA::magDiff: Dimensions do not match! (w1: %u, w2: %u, h1: %u, h2: %u)\n", img1.width, img2.width, img1.height, img2.height);
		return;
	}

	uint color1, color2, color;
	unsigned char r1, r2, g1, g2, b1, b2, a1, a2;
	unsigned char r, g, b, a;
	uint width 	= img1.width;
	uint height	= img1.height;

	init(width, height);
	uint count = 0;

	for(uint x=0; x<width; ++x)
		for(uint y=0; y<height; ++y)
		{
			color1 = img1(x,y);
			color2 = img2(x,y);

			RGBA82char(color1, r1, g1, b1, a1);
			RGBA82char(color2, r2, g2, b2, a2);

			if(r1==r2 && g1==g2 && b1==b2)
			{
				r = 0;
				g = 0;
				b = 0;
				a = 255;
			}
			else
			{
				r = r1;
				g = g1;
				b = b1;
				a = 255;
				++count;
			}

			color = char2RGBA8(r, g, b, a);

			this->operator ()(x, y) = color;
		}

	//printf("hcImageRGBA::magDiff: %u pixels differ (= %E of entire image)\n", count, (hcFloat)count / (width*height));
}

void hcImageRGBA::setBackgroundColor(uint bgColor, int insertPosX, int insertPosY, uint subWidth, uint subHeight)
{
    if(subWidth == 0)
        subWidth = width;

    if(subHeight == 0)
        subHeight = height;

    uint heightTo = insertPosY + subHeight;
    uint widthTo  = insertPosX + subWidth;

    if(heightTo > height)
        heightTo = height;

    if(widthTo > width)
        widthTo = width;

    for(int y=insertPosY; y<heightTo; ++y)
        for(int x=insertPosX; x<widthTo; ++x)
            if(!((y >= height) || (y < 0) || (x >= width) || (x < 0)))
                this->data[y * width + x] = bgColor;
}

bool hcImageRGBA::insertSubimage(hcImage &src, int insertPosX, int insertPosY)
{
    float alpha;
    uint colorDst, colorSrc;
    unsigned char rDst, gDst, bDst, aDst;
    unsigned char rSrc, gSrc, bSrc, aSrc;

    for(int x=0; x<src.width; ++x)
        for(int y=0; y<src.height; ++y)
        {
            // if position outside of destination image boundaries, skip
            if((y + insertPosY >= height) || (y + insertPosY < 0) || (x + insertPosX >= width) || (x + insertPosX < 0))
                continue;

            colorDst   	= data[(y + insertPosY) * width + (x + insertPosX)];
            colorSrc	= src(x,y);

            RGBA82char(colorSrc, rSrc, gSrc, bSrc, aSrc);
            RGBA82char(colorDst, rDst, gDst, bDst, aDst);

            alpha           = aSrc / 255.0;
            uint finalAlpha = aDst > aSrc ? aDst : aSrc;

            data[(y + insertPosY) * width + (x + insertPosX)] = \
                    char2RGBA8( (1 - alpha) * rDst + alpha * rSrc,\
                                (1 - alpha) * gDst + alpha * gSrc,\
                                (1 - alpha) * bDst + alpha * bSrc,\
                                finalAlpha
                                );
        }

    return true;
}

/*! \brief resizes a rectangular image via bilinear interpolation TODO: redo this function!
 *
 */
bool hcImageRGBA::interpolateRectangularImage(uint width_out, uint height_out)
{
    if(data == NULL || width == 0 || height == 0)
    {
        printf("ERROR! hcImageRGBA::interpolateRectangularImage(): Image is empty!");
        printf("Width: %u, height: %u\n", width, height);
        return false;
    }

    uint *image_out         = new uint[width_out * height_out];

    uint i_out, j_out;

    uint newWidthIn     	= 1 << nearestPowerOfTwo(width);
	uint newHeightIn    	= 1 << nearestPowerOfTwo(height);

    uint *tempImage         = new uint[newWidthIn * newHeightIn];
    uint *tempImage2        = new uint[newWidthIn * newHeightIn];

    float scaleFactorWidth  = (float)width  / width_out;
    float scaleFactorHeight = (float)height / height_out;

    // if input image is twice as wide / tall, it has to be treated special to not
    // lose energy in image space (due to unsampled pixels in input image)
    if(scaleFactorWidth >= 2 || scaleFactorHeight >= 2)
    {
        interpolateRectangularImage_orig(data, tempImage, width, height, newWidthIn, newHeightIn);

        scaleFactorWidth    = (float)newWidthIn 	/ width_out;
        scaleFactorHeight   = (float)newHeightIn 	/ height_out;

        while(scaleFactorWidth >= 2 || scaleFactorHeight >= 2)
        {
            unsigned char r1 = 0;
            unsigned char r2 = 0;
            unsigned char g1 = 0;
            unsigned char g2 = 0;
            unsigned char b1 = 0;
            unsigned char b2 = 0;
            unsigned char a1 = 0;
            unsigned char a2 = 0;

            bool signal = false;

            if(scaleFactorWidth >= 2)
            {
                signal = true;

                for(uint i=0;i<newWidthIn/2;++i)
                    for(uint j=0;j<newHeightIn;++j)
                    {
                        RGBA82char(tempImage[j * newWidthIn + 2*i],     r1, g1, b1, a1);
                        RGBA82char(tempImage[j * newWidthIn + 2*i + 1], r2, g2, b2, a2);
                        if(a2 > a1)
							a1 = a2;
                        tempImage2[j * newWidthIn/2 + i] = char2RGBA8(0.5*(r1+r2), 0.5*(g1+g2), 0.5*(b1+b2), a1);
                    }

                newWidthIn 			= newWidthIn / 2;
                scaleFactorWidth 	= (float) newWidthIn / width_out;
            }

            if(scaleFactorHeight >= 2)
            {
                if(signal)
                {
                    for(uint i=0;i<newWidthIn;++i)
                        for(uint j=0;j<newHeightIn;++j)
                            tempImage[j*newWidthIn+i] = tempImage2[j*newWidthIn+i];
                }

                for(uint i=0;i<newWidthIn;++i)
                    for(uint j=0;j<newHeightIn/2;++j)
                    {
                        RGBA82char(tempImage[(2*j) * newWidthIn + i],   r1, g1, b1, a1);
                        RGBA82char(tempImage[(2*j+1) * newWidthIn + i], r2, g2, b2, a2);
                        if(a2 > a1)
							a1 = a2;
                        tempImage2[j * newWidthIn + i] = char2RGBA8(0.5*(r1+r2), 0.5*(g1+g2), 0.5*(b1+b2), a1);
                    }

                newHeightIn = newHeightIn / 2;
                scaleFactorHeight = (float) newHeightIn / height_out;
            }

            // image still more than twice as wide/tall ? copy new image and repeat
            // : set new width / height, repoint data and let the
            // algorithm work its magic
            if (scaleFactorWidth >= 2 || scaleFactorHeight >= 2)
            {
                for(uint i=0;i<newWidthIn;++i)
                    for(uint j=0;j<newHeightIn;++j)
                         tempImage[j*newWidthIn+i] = tempImage2[j*newWidthIn+i];
            }
            else
            {
            	init(newWidthIn, newHeightIn);

            	for(uint i=0;i<newWidthIn;++i)
					for(uint j=0;j<newHeightIn;++j)
						data[j*newWidthIn+i] = tempImage2[j*newWidthIn+i];
            }
        }
    }

    delete [] tempImage;
	delete [] tempImage2;

#define max($1, $2) $1 > $2 ? $1 : $2

    for(i_out = 0; i_out < width_out; ++i_out)
        for(j_out = 0; j_out < height_out; ++j_out)
        {
            int i  = floor((scaleFactorWidth   * (i_out + 0.5)) - 0.5);
            int j  = floor((scaleFactorHeight  * (j_out + 0.5)) - 0.5);
            if(i<0)
                i=0;
            if(j<0)
                j=0;

            int ind_i   = i;
            int ind_ip  = i+1 < (signed)width  ? i+1 : i;
            int ind_j   = j;
            int ind_jp  = j+1 < (signed)height ? j+1 : j;


            float xOut  = (i_out    + 0.5) / width_out;
            float yOut  = (j_out    + 0.5) / height_out;
            float xi    = (i        + 0.5) / width;
            float xip   = (i + 1    + 0.5) / width;
            float yi    = (j        + 0.5) / height;
            float yip   = (j + 1    + 0.5) / height;

            unsigned char r_ij      = 0;
            unsigned char r_ipj     = 0;
            unsigned char r_ijp     = 0;
            unsigned char r_ipjp    = 0;

            unsigned char g_ij      = 0;
            unsigned char g_ipj     = 0;
            unsigned char g_ijp     = 0;
            unsigned char g_ipjp    = 0;

            unsigned char b_ij      = 0;
            unsigned char b_ipj     = 0;
            unsigned char b_ijp     = 0;
            unsigned char b_ipjp    = 0;

            unsigned char a_ij      = 0;
            unsigned char a_ipj     = 0;
            unsigned char a_ijp     = 0;
            unsigned char a_ipjp    = 0;

            RGBA82char(data[ind_j   * width  + ind_i],  r_ij,   g_ij,   b_ij,   a_ij);
            RGBA82char(data[ind_j   * width  + ind_ip], r_ipj,  g_ipj,  b_ipj,  a_ipj);
            RGBA82char(data[ind_jp  * width  + ind_i],  r_ijp,  g_ijp,  b_ijp,  a_ijp);
            RGBA82char(data[ind_jp  * width  + ind_ip], r_ipjp, g_ipjp, b_ipjp, a_ipjp);

            float a1 	= (xip - xOut) / (xip - xi);
            float a2 	= 1 - a1;

            float rIntL = a1 * r_ij + a2 * r_ipj;
            float gIntL = a1 * g_ij + a2 * g_ipj;
            float bIntL = a1 * b_ij + a2 * b_ipj;
            float aIntL = max(a_ij, a_ipj);

            float rIntU = a1 * r_ijp + a2 * r_ipjp;
            float gIntU = a1 * g_ijp + a2 * g_ipjp;
            float bIntU = a1 * b_ijp + a2 * b_ipjp;
            float aIntU = max(a_ijp, a_ipjp);

            a1 			= (yip - yOut) / (yip - yi);
            a2 			= 1 - a1;

            float red   = a1 * rIntL + a2 * rIntU;
            float green = a1 * gIntL + a2 * gIntU;
            float blue  = a1 * bIntL + a2 * bIntU;
            float alpha = a1 * aIntL + a2 * aIntU;

            image_out[j_out * width_out + i_out] = char2RGBA8(red, green, blue, alpha);
        }

    init(width_out, height_out);
    for(uint i=0;i<width_out;++i)
		for(uint j=0;j<height_out;++j)
			data[j*width_out+i] = image_out[j*width_out+i];

    delete [] image_out;

    return true;
}


bool hcImageRGBA::save(const string &filename)
{
	FIBITMAP *bitmap	= FreeImage_Allocate(width, height, 32);
	uint inRowBytes		= width * sizeof(uint);

	for (int y=0; y<height; ++y)
		for(int x=0; x<width; ++x)
		{
			RGBQUAD val;
			uint d = operator()(x,y);
			unsigned char r, g, b, a;
			RGBA82char(d, r, g, b, a);
			val.rgbRed 		= r;
			val.rgbGreen 	= g;
			val.rgbBlue 	= b;
			FreeImage_SetPixelColor(bitmap, x, y, &val);
		}

	FreeImage_Save(FIF_BMP, bitmap, filename.data());
	FreeImage_Unload(bitmap);

    return true;
}

#ifdef GUI
bool hcImageRGBA::load(const string &filename) //TODO get rid of SOIL and use FreeImage
{
    if(!checkFileEx(filename, "hcImageRGBA::load"))
    {
    	cerr << __FILE__ << ":" << __LINE__ << ": File '" << filename << "' does not exist.\n";
        return false;
    }

    int channels, width_in, height_in;

    unsigned char *imageData   = SOIL_load_image(filename.data(), &width_in, &height_in, &channels, SOIL_LOAD_AUTO );
    init(width_in, height_in);

    if(channels != 4 && channels != 3)
    {
        printf("ERROR! hcImageRGBA::load(): number of channels in image not 4 or 3 (is %i)!\n", channels);
        free(imageData);
        return false;
    }

    for(uint x=0; x<width; ++x)
        for(uint y=0; y<height; ++y)
        {
        	if(channels == 4)
				data[(height - y - 1) * width + x] =
						char2RGBA8(imageData[(y*width+x)*4 + 0], imageData[(y*width+x)*4 + 1],
								   imageData[(y*width+x)*4 + 2], imageData[(y*width+x)*4 + 3]);
			else if(channels == 3)
				data[(height - y - 1) * width + x] =
						char2RGBA8(imageData[(y*width+x)*3 + 0], imageData[(y*width+x)*3 + 1],
								   imageData[(y*width+x)*3 + 2], 255);
            else
            {
            	printf("ERROR! hcImageRGBA::load: Number of channels (%i) not supported!\n", channels);
            	return false;
            }
        }

    free(imageData);

    return true;
}
#endif

bool hcImageRGBA::loadFromArray(uint *array, uint width, uint height)
{
    init(width, height);

    for(uint x=0; x<width; ++x)
        for(uint y=0; y<height; ++y)
            data[y * width + x] = array[y * width + x];
    return true;
}


void hcImageRGBA::dump() const
{
    printf("--- Dumping hcImageRGBA:\n");
    hcImage<uint>::dump();
}

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//
//                          hcImageFloat
//
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------


hcImageFloat::hcImageFloat()
{
    initNULL();
}

hcImageFloat::hcImageFloat(const hcImageFloat &other)
{
    initNULL();
    hcImage<hcFloat>::operator=(other);
}

hcImageFloat::hcImageFloat(uint width, uint height)
{
    initNULL();
    init(width, height);
}

hcImageFloat::hcImageFloat(uint width, uint height, float color)
{
    initNULL();
    init(width, height, color);
}

hcImageFloat::~hcImageFloat(){

    clear();
}

hcImageFloat &hcImageFloat::operator=(const hcImageFloat &other)
{
    if(this == &other)
        return *this;

    hcImage<hcFloat>::operator=(other);

    return *this;
}

void hcImageFloat::init(uint width, uint height)
{
    hcImage<hcFloat>::init(width, height);
}

void hcImageFloat::init(uint width, uint height, float color)
{
    clear();
    init(width, height);
    setBackgroundColor(color);
}

void hcImageFloat::reldiffImage(const hcImageFloat &img0, const hcImageFloat &img1)
{
	if(img0.width != img1.width || img0.height != img1.height)
	{
		printf("ERROR! hcImageFloat::diffImage: Dimensions do not match!\nimg0.width: %u, img1.width: %u, img0.height: %u, img1.height: %u\n", img0.width, img1.width, img0.height, img1.height);
		return;
	}

	uint width 	= img0.width;
	uint height	= img0.height;

	init(width, height);

	for(uint x=0; x<width; ++x)
		for(uint y=0; y<height; ++y)
			this->operator ()(x,y) = fabs(img0.data[y * width + x] - img1.data[y * width + x])/fabs(img0.data[y * width + x]);
}

void hcImageFloat::diffImage(const hcImageFloat &img0, const hcImageFloat &img1)
{
	if(img0.width != img1.width || img0.height != img1.height)
	{
		printf("ERROR! hcImageFloat::diffImage: Dimensions do not match!\nimg0.width: %u, img1.width: %u, img0.height: %u, img1.height: %u\n", img0.width, img1.width, img0.height, img1.height);
		return;
	}

	uint width 	= img0.width;
	uint height	= img0.height;

	init(width, height);

	for(uint x=0; x<width; ++x)
		for(uint y=0; y<height; ++y)
			this->operator ()(x,y) = img0.data[y * width + x] - img1.data[y * width + x];
}

float hcImageFloat::absMean()
{
	float retval = 0.0;

	for(uint x=0; x<width; ++x)
		for(uint y=0; y<height; ++y)
			retval += fabs(this->operator ()(x,y));

	retval /= width*height;

	return retval;
}

void hcImageFloat::setBackgroundColor(float bgColor, int insertPosX, int insertPosY, uint subWidth, uint subHeight)
{
    if(subWidth == 0)	subWidth = width;
    if(subHeight == 0)	subHeight = height;

    int heightTo = insertPosY + subHeight;
    int widthTo  = insertPosX + subWidth;

    if(heightTo > height)	heightTo 	= height;
    if(widthTo  > width)	widthTo 	= width;

    for(int y=insertPosY; y<heightTo; ++y)
        for(int x=insertPosX; x<widthTo; ++x)
            if(!((y >= height) || (y < 0) || (x >= width) || (x < 0)))
                this->data[y * width + x] = bgColor;
}//*/

/*! inserts image supplied by src at the supplied positions by replacing the data stored
 *  at these positions in this (opaque, no transperency)
 */
bool hcImageFloat::insertSubimage(const hcImageFloat &src, int insertPosX, int insertPosY)
{
    for(int x=0; x<src.width; ++x)
        for(int y=0; y<src.height; ++y)
        {
            // if position outside of destination image boundaries, skip
            if((y + insertPosY >= height) || (y + insertPosY < 0) || (x + insertPosX >= width) || (x + insertPosX < 0))
                continue;

            data[(y + insertPosY) * width + (x + insertPosX)] = src.data[y * src.width + x];
        }

    return true;
}


bool hcImageFloat::save(const string &filename)
{
    printf("ERROR! hcImageFloat::save: Not implemented yet!\n");
    //SOIL_save_image(filename, SOIL_SAVE_TYPE_BMP, width, height, 4, reinterpret_cast<unsigned char*>(data));

    return true;
}

bool hcImageFloat::load(const string &filename)
{
    if(!checkFileEx(filename, "hcImageFloat::load"))
        return false;

    printf("ERROR! hcImageFloat::load: Not implemented yet!\n");

    return true;
}

bool hcImageFloat::loadFromArray(float *array, uint width, uint height)
{
    init(width, height);

    for(uint x=0; x<width; ++x)
        for(uint y=0; y<height; ++y)
            data[y * width + x] = array[y * width + x];
    return true;
}

#ifdef GUI
void hcImageFloat::addHomWhiteNoise(float sigma, uint seed)
{
	gsl_rng *rng = gsl_rng_alloc(gsl_rng_taus2);
	gsl_rng_set(rng, seed);

	for(uint x=0; x<width; ++x)
		for(uint y=0;y<height; ++y)
			data[y * width + x] += gsl_ran_gaussian(rng, sigma);

	gsl_rng_free(rng);
}

void hcImageFloat::addPixelNoise(float fraction, uint seed)
{
	gsl_rng *rng = gsl_rng_alloc(gsl_rng_taus2);
	gsl_rng_set(rng, seed);

	for(uint x=0; x<width; ++x)
		for(uint y=0;y<height; ++y)
		{
			uint ind 	= y * width + x;
			float sigma = fabs(data[ind]) * fraction;
			data[ind] 	+= gsl_ran_gaussian(rng, sigma);
		}

	gsl_rng_free(rng);
}

void hcImageFloat::addSignedFractionNoise(float fraction, uint seed)
{
	gsl_rng *rng = gsl_rng_alloc(gsl_rng_taus2);
	gsl_rng_set(rng, seed);

	for(uint x=0; x<width; ++x)
		for(uint y=0;y<height; ++y)
		{
			uint ind 		= y * width + x;
			double pixVal	= data[ind];
			double rand 	= gsl_rng_uniform(rng);
			int sign		= rand < 0.5 ? -1 : +1;
			double retval	= sign * pixVal * fraction;
			data[ind]		+= retval;
		}

	gsl_rng_free(rng);
}
#endif

hcFloat hcImageFloat::meanSquaredDiff(hcImage &other)
{
	double difference = 0.0;

	for(uint x=0; x<width; ++x)
		for(uint y=0; y<height; ++y)
		{
			uint ind 	= y*width+x;
			difference 	+= (operator()(x,y) - other(x,y)) * (operator()(x,y) - other(x,y));
		}

	difference /= width * height;
	return difference;
}

void hcImageFloat::meanFilter(uint windowWidth, bool circular)
{
	hcImageFloat temp(width, height);

	for(int x=0; x<width; ++x)
		for(int y=0; y<height; ++y)
		{
			float mean = 0.0;
			for(int xx=-(int)windowWidth/2; xx<=(int)windowWidth/2; ++xx)
				for(int yy=-(int)windowWidth/2; yy<=(int)windowWidth/2; ++yy)
				{
					uint xInd 	= (x + xx < 0 		? (circular ? x + xx + width  : 0       ) :
								  (x + xx >= width 	? (circular ? x + xx - width  : width-1 ) : x + xx));
					uint yInd	= (y + yy < 0		? 0 									  :
								  (y + yy >=height  ? height-1 								  : y + yy));

					uint ind	= yInd * width + xInd;
					mean 		+= data[ind];
				}
			uint numPixels 	= (2 * windowWidth/2 + 1) * (2 * windowWidth/2 + 1);
			mean 			/= numPixels;
			temp(x,y) 		= mean;
		}

	this->operator =(temp);
}

void sortBubble(hcFloat *array, uint numElements, bool ascending=true)
{
	if(numElements <= 1)
		return;

	bool changed = true;

	while(changed)
	{
		changed = false;
		for(uint i=1; i<numElements; ++i)
		{
			if(ascending ? array[i-1] > array[i] : array[i-1] < array[i])
			{
				hcFloat temp 	= array[i-1];
				array[i-1]		= array[i];
				array[i]		= temp;
				changed 		= true;
			}
		}
	}
}

void hcImageFloat::medianFilter(uint windowWidth, bool circular)
{
	hcImageFloat temp(width, height);

	uint numPixels 	= (windowWidth + 1) * (windowWidth + 1);
	hcFloat sortArray[numPixels];

	for(int x=0; x<width; ++x)
		for(int y=0; y<height; ++y)
		{
			uint i = 0;
			for(int xx=-(int)windowWidth/2; xx<=(int)windowWidth/2; ++xx)
				for(int yy=-(int)windowWidth/2; yy<=(int)windowWidth/2; ++yy)
				{
					uint xInd 	= (x + xx < 0 		? (circular ? x + xx + width  : 0       ) :
								  (x + xx >= width 	? (circular ? x + xx - width  : width-1 ) : x + xx));
					uint yInd	= (y + yy < 0		? 0 									  :
								  (y + yy >=height  ? height-1 								  : y + yy));

					uint ind	= yInd * width + xInd;

					sortArray[i++] = data[ind];
				}

			sortBubble(sortArray, numPixels);

			temp(x,y) 		= numPixels % 2 == 0 ? sortArray[numPixels/2] : 0.5 * (sortArray[numPixels/2] + sortArray[numPixels/2+1]);
		}

	this->operator =(temp);
}

/*
void hcImageFloat::rescale(uint newWidth, uint newHeight)
{
	hcFloat *dataNew = new hcFloat[newWidth*newHeight];

	rescaleImage(data, width, height, dataNew, newWidth, newHeight, FILTER_LANCZOS3, true);
	init(newWidth, newHeight);
	delete [] data;
	data = dataNew;
}//*/


void hcImageFloat::dump() const
{
    printf("--- Dumping hcImageFloat:\n");
    hcImage<hcFloat>::dump();
}

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//
//                          hcImageVec3D
//
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------


hcImageVec3D::hcImageVec3D()
{
    initNULL();
}

hcImageVec3D::hcImageVec3D(const hcImageVec3D &other)
{
    initNULL();
    hcImage<Vec3D>::operator=(other);
}

hcImageVec3D::hcImageVec3D(uint width, uint height)
{
    initNULL();
    init(width, height);
}

hcImageVec3D::~hcImageVec3D()
{
    clear();
}

hcImageVec3D &hcImageVec3D::operator=(const hcImageVec3D &other)
{
    if(this == &other)
        return *this;

    hcImage<Vec3D>::operator=(other);

    return *this;
}

void hcImageVec3D::init(uint width, uint height)
{
    hcImage<Vec3D>::init(width, height);
}

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//
//                          hcImageInt
//
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------


hcImageInt::hcImageInt()
{
    initNULL();
}

hcImageInt::hcImageInt(const hcImageInt &other)
{
    initNULL();
    hcImage<int>::operator=(other);
}

hcImageInt::hcImageInt(uint width, uint height)
{
    initNULL();
    init(width, height);
}

hcImageInt::~hcImageInt()
{
    clear();
}

hcImageInt &hcImageInt::operator=(const hcImageInt &other)
{
    if(this == &other)
        return *this;

    hcImage<int>::operator=(other);

    return *this;
}

void hcImageInt::init(uint width, uint height)
{
    hcImage<int>::init(width, height);
}

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//
//                          hcImageBool
//
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------


hcImageBool::hcImageBool()
{
    initNULL();
}

hcImageBool::hcImageBool(const hcImageBool &other)
{
    initNULL();
    hcImage<bool>::operator=(other);
}

hcImageBool::hcImageBool(uint width, uint height)
{
    initNULL();
    init(width, height);
}

hcImageBool::~hcImageBool()
{
    clear();
}

hcImageBool &hcImageBool::operator=(const hcImageBool &other)
{
    if(this == &other)
        return *this;

    hcImage<bool>::operator=(other);

    return *this;
}

void hcImageBool::init(uint width, uint height)
{
    hcImage<bool>::init(width, height);
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

hcScribble::hcScribble()
{
    initNULL();
}

hcScribble::hcScribble(const hcScribble &other)
{
    initNULL();
    *this = other;
}

hcScribble::~hcScribble()
{
	clear();
}

hcScribble &hcScribble::operator=(const hcScribble &other)
{
    if(this == &other)
        return *this;

    hcImageRGBA::operator=(other) ;
    original            = other;
    this->insertPosX    = other.insertPosX;
    this->insertPosY    = other.insertPosY;
    this->relWidth      = other.relWidth;
    this->newRelPosX    = other.newRelPosX;
    this->newRelPosY    = other.newRelPosY;

    return *this;
}

void hcScribble::initNULL()
{
    canvas          = NULL;
    insertPosX      = 0;
    insertPosY      = 0;
    relWidth        = 0.0;
    relPosX         = 0.0;
    relPosY         = 0.0;
    newRelPosX      = 0.0;
    newRelPosY      = 0.0;
    hcImageRGBA::init(0, 0);
    original.init(0, 0);
    canvasOrig.init(0, 0);
}

void hcScribble::clear()
{
    initNULL();
}

void hcScribble::init(double relWidth, double relPosX, double relPosY, hcImageRGBA *canvas)
{
    clear();
    this->canvas        = canvas;
    this->relWidth      = relWidth;
    this->relPosX       = -2.0;
    this->relPosY       = -2.0;
    this->newRelPosX    = relPosX;
    this->newRelPosY    = relPosY;
}

void hcScribble::setCanvas(hcImageRGBA *canvas)
{
	this->canvas	= canvas;
}

void hcScribble::reinit()
{
    canvasOrig.init(0, 0);
}

bool hcScribble::updateScaledScribble()
{
    if(canvas == NULL)
        return false;

    double aspectRatio  = (double) original.width / original.height;
    uint pxWidth        = relWidth  * canvas->width  / 2.0;
    uint pxHeight       = pxWidth / aspectRatio;
    uint posX           = posTransform(-1, 1, 0, canvas->width,  relPosX);
    uint posY           = posTransform(-1, 1, 0, canvas->height, relPosY);
    insertPosX          = posX - pxWidth  / 2.0;
    insertPosY          = posY - pxHeight / 2.0;

    ((hcImageRGBA*)this)->operator=(original);

    canvasOrig.init(pxWidth, pxHeight);
    *((hcImageRGBA*)(this)) = original;
    this->interpolateRectangularImage(pxWidth, pxHeight);

    return true;
}

bool hcScribble::updatePosition()
{
    if(canvas == NULL)
        return false;

    double aspectRatio  = (double) original.width / original.height;
    uint pxWidth        = relWidth  * canvas->width  / 2.0;
    uint pxHeight       = pxWidth / aspectRatio;
    uint posX           = posTransform(-1, 1, 0, canvas->width,  relPosX);
    uint posY           = posTransform(-1, 1, 0, canvas->height, relPosY);
    insertPosX          = posX - pxWidth  / 2.0;
    insertPosY          = posY - pxHeight / 2.0;

    return true;
}

bool hcScribble::draw()
{
	if(canvas == NULL)
		return false;

	bool posBool 	= newRelPosX == relPosX && newRelPosY == relPosY;
	bool initiated	= canvasOrig.width != 0;

	if(posBool && initiated)
        return false;

	initiated = redrawOrig();

    this->relPosX = newRelPosX;
    this->relPosY = newRelPosY;

    if(!initiated)
        updateScaledScribble();
    else
        updatePosition();

    canvas->extractSubimage(insertPosX, insertPosY, canvasOrig);
    canvas->insertSubimage(*this, insertPosX, insertPosY);

    return true;
}

bool hcScribble::redrawOrig()
{
    if(canvasOrig.width == 0)
        return false;

    canvas->insertSubimage(canvasOrig, insertPosX, insertPosY);
    return true;
}

void hcScribble::dump() const
{
    printf("--- Dumping hcScribble:\n");
    printf("insertPosX/Y:\t%u/%u\n", insertPosX, insertPosY);
    printf("relPosX/Y:\t\t%E / %E\n", relPosX, relPosY);
    printf("newRelPosX/Y:\t%E / %E\n", newRelPosX, newRelPosY);
    printf("relWidth:\t\t%E\n", relWidth);
    //printf("Original:\n");
    //original.dump();
    //printf("Rescaled:\n");
    //hcImageRGBA::dump();
}

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//
//                          hcScribbleVertLine
//
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

hcScribbleVertLine::hcScribbleVertLine()
{
    initNULL();
    init(char2RGBA8(0, 0, 0, 255), 0.1, 0.0, 0.0, NULL);
}

hcScribbleVertLine::hcScribbleVertLine(const hcScribbleVertLine &other)
{
    initNULL();
    *this = other;
}

hcScribbleVertLine::hcScribbleVertLine(uint color, double relWidth, double relPosX, double relPosY, hcImageRGBA *canvas)
{
    initNULL();
    init(color, relWidth, relPosX, relPosY, canvas);
}

hcScribbleVertLine::~hcScribbleVertLine()
{
    clear();
}

hcScribbleVertLine &hcScribbleVertLine::operator=(const hcScribbleVertLine &other)
{
    if(this == &other)
        return *this;

    hcScribble::operator=(other);

    return *this;
}

void hcScribbleVertLine::initNULL(){}

void hcScribbleVertLine::clear()
{
    initNULL();
}

void hcScribbleVertLine::init(uint color, double relWidth, double relPosX, double relPosY, hcImageRGBA *canvas)
{
    hcScribble::init(relWidth, relPosX, relPosY, canvas);
    original.init(1, 100, color);
}

void hcScribbleVertLine::dump() const
{
    printf("--- Dumping hcScribbleVertLine:\n");
    hcScribble::dump();
}

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//
//                          hcScribbleDot
//
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

hcScribbleDot::hcScribbleDot()
{
    initNULL();
    init(char2RGBA8(0, 0, 0, 255), 0.1, 0.0, 0.0, NULL);
}

hcScribbleDot::hcScribbleDot(const hcScribbleDot &other)
{
    initNULL();
    *this = other;
}

hcScribbleDot::hcScribbleDot(uint color, double relWidth, double relPosX, double relPosY, hcImageRGBA *canvas)
{
    initNULL();
    init(color, relWidth, relPosX, relPosY, canvas);
}

hcScribbleDot::~hcScribbleDot()
{
    clear();
}

hcScribbleDot &hcScribbleDot::operator=(const hcScribbleDot &other)
{
    if(this == &other)
        return *this;

    hcScribble::operator=(other);

    return *this;
}

void hcScribbleDot::initNULL(){}

void hcScribbleDot::clear()
{
    initNULL();
}

void hcScribbleDot::init(uint color, double relWidth, double relPosX, double relPosY, hcImageRGBA *canvas)
{
    hcScribble::init(relWidth, relPosX, relPosY, canvas);
    original.init(1, 1, color);
}

void hcScribbleDot::dump() const
{
    printf("--- Dumping hcScribbleDot:\n");
    hcScribble::dump();
}
