//To compile and run:
// On os x:
// g++ coordTransform.cpp -g -I ./halide/include -I stack/DarwinX86/daf_base/11.0/include/ -I stack/DarwinX86/ndarray/10.1+24/include/ -I stack/DarwinX86/pex_exceptions/10.0+37/include/ -I stack/DarwinX86/pex_policy/11*/include/ -I stack/DarwinX86/daf_persistence/11.0-1-g7418c06/include/ -I stack/DarwinX86/utils/11*/include/ -I stack/DarwinX86/daf_base/10.1+2/include/ -I stack/DarwinX86/base/11*/include/ -I stack/DarwinX86/ndarray/10.1+2/include/ -I stack/DarwinX86/pex_exceptions/11*/include/ -I stack/DarwinX86/eigen/3.2.0/include/ -I stack/DarwinX86/afw/11*/include -L ./halide/bin -L stack/DarwinX86/afw/11*/lib -L stack/DarwinX86/daf_persistence/11.0-1-g7418c06/lib/ -L stack/DarwinX86/boost/1.59.lsst3/lib/ -L stack/DarwinX86/pex_exceptions/10.0+37/lib/ -L stack/DarwinX86/daf_base/11.0/lib/ -L stack/DarwinX86/pex_policy/11.0/lib/ -lHalide -lafw -ldaf_base -ldaf_persistence -lboost_system -o coordTransform -std=c++11
//
// DYLD_LIBRARY_PATH=./halide/bin:stack/DarwinX86/afw/11.0-7-g974e250/lib:stack/DarwinX86/daf_persistence/11.0-1-g7418c06/lib/:stack/DarwinX86/boost/1.59.lsst3/lib/:stack/DarwinX86/pex_exceptions/10.0+37/lib/:stack/DarwinX86/daf_base/11.0/lib/:stack/DarwinX86/pex_policy/11.0/lib/:stack/DarwinX86/pex_logging/11.0/lib/:stack/DarwinX86/utils/11.0-1-g47edd16/lib/:stack/DarwinX86/base/11.0/lib/ ./coordTransform

#include "lsst/afw/image.h"
#include "lsst/afw/geom.h"

namespace afwGeom = lsst::afw::geom;
namespace afwImage = lsst::afw::image;
namespace afwMath  = lsst::afw::math;

#include <stdio.h>
#include "Halide.h"
#include <bitset>
#include "clock.h"
using namespace std;
using namespace Halide;

int main(int argc, char *argv[]) {
    //Use the LSST .fits reader to create LSST MaskedImage types of the source
    //and destination images
    auto destImage = afwImage::MaskedImage<float>("./images/calexp-0289820_24.fits");
    auto srcImage = afwImage::MaskedImage<float>("./images/calexp-0288976_24.fits");

    //print out the size of the source image
    int width = srcImage.getWidth(), height = srcImage.getHeight();
    printf("Loaded: %d x %d\n", width, height);

    //Read destination image values into a Halide image
    Image<float> destImageVal(width, height);
    for (int y = 0; y < height; y++) {
        afwImage::MaskedImage<float, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>::x_iterator inPtr = destImage.x_at(0, y);
        for (int x = 0; x < width; x++){
            destImageVal(x, y) = (*inPtr).image();
            inPtr++;
        }
    }


    //Use LSST methods to extract coordinate information from the source
    //and destination images
    auto destFitsHeader = afwImage::readMetadata("./images/calexp-0289820_24.fits",1);
    auto destWcs = afwImage::makeWcs(destFitsHeader);

    auto srcFitsHeader = afwImage::readMetadata("./images/calexp-0288976_24.fits",1);
    auto srcWcs = afwImage::makeWcs(srcFitsHeader);
  
    //Time computing the coordinate transform for every pixel in the image
    //to show it takes a while
    clock_t t1 = clock();
    for(double destCol = 0; destCol < 2048; destCol++){
        for(double destRow = 0; destRow < 1489; destRow++){
            //The next 3 lines are where the coordinate transform actually happens
            afwGeom::Angle sky1, sky2;
            destWcs->pixelToSky(destCol, destRow, sky1, sky2);
            afwGeom::Point2D srcPosition = srcWcs->skyToPixel(sky1, sky2);
        }
    }
    clock_t t2 = clock();
    printf("Calculating coordinate transform for every point in the image took %f seconds\n", (float)(t2-t1)/CLOCKS_PER_SEC);


    //Sometimes images may be shifted such that their (0, 0) pixel location does
    //not line up with the (0, 0) locaiton of their coordinate system.
    //I think this is primarily used when computing on an image that is a subregion
    //of a larger image (so the cropped image is actually in the middle of the
    //larger image's coordinate system)
    //Let's check our images.  It turns out they aren't shifted.
    afwGeom::Point2D const destXY0(destImage.getXY0());
    afwGeom::Point2D const srcXY0(srcImage.getXY0());
    printf("destXY0[0] = %f, destXY0[1] = %f\n", destXY0[0], destXY0[1]);
    printf("srcXY0[0] = %f, srcXY0[1] = %f\n", srcXY0[0], srcXY0[1]);


    //Now we'll do the coordinate transform for a particular point in the destination
    //image that is in the center of a star to make sure everything works correctly. 
    double destCol = destXY0[0] + 237;
    double destRow = destXY0[1] + 67;
    afwGeom::Angle sky1, sky2;
    destWcs->pixelToSky(destCol, destRow, sky1, sky2);
    afwGeom::Point2D srcPosition = srcWcs->skyToPixel(sky1, sky2);
    printf("dest location (%f, %f) converted to source location (%f, %f)\n",
            destCol, destRow,
            srcPosition[0] - srcXY0[0], srcPosition[1] - srcXY0[1]);
    printf("destImage(237, 67) = %f\n", destImageVal(237, 67));
    //You can see that the converted source location falls in the center of the
    //same star in the source image.  Also, running testWarp.py, you can see
    //that in the warped output image this star has roughly the same position as in 
    //the destination image and maximum intensity as in the source image. 

}

