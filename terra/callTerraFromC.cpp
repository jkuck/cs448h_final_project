//To Run:
//$ terra blurExampleTerraStandalone.t
//$ terra-Linux-x86_64-84bbb0b/bin/terra blurExampleTerraStandalone.t
//$ g++ -c callTerraFromC.cpp
//$ g++ blurExampleTerraStandalone.o callTerraFromC.o
//$ ./a.out
#include "blurExampleTerraStandalone.h"
#include <math.h>       /* cos */
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */


//#define TEST_VS_LSST
//#define PRINT_CORNERS_OF_OUTPUT_AND_REFERENCE //only define if TEST_VS_LSST is defined
//#define USE_5_KERNELS 

//to test double precision kernel define this
//also define USE_5_KERNELS and 
//#define DOUBLE_PRECISION_KERNEL

//can be one of the following strings:
//"gaussian_contains_denormal_numbers"
//"gaussian_denormals_zeroed" 
//"random" //output will not match LSST, but for checking performance
#define KERNEL_INFO "random"

//smallest kernel value if gaussian_denormals_zeroed
#define MIN_KERNEL_VAL pow(10, -30) 

//count kernel vals below this number as denormals
#define COUNT_AS_DENORMAL pow(10, -38)

float kernel1(int i, int j){
    float sigmaX1 = 2.0;
    float sigmaY1 = 2.0;
    float theta1 = 0.0; //rotation of sigmaX axis
    return (exp(-((i*cos(theta1) +j*sin(theta1))*(i*cos(theta1) +j*sin(theta1)))
                    /(2*sigmaX1*sigmaX1)) / (sqrt(2*M_PI)*sigmaX1))
                    *(exp(-((j*cos(theta1) - i*sin(theta1))*(j*cos(theta1) - i*sin(theta1)))
                    /(2*sigmaY1*sigmaY1)) / (sqrt(2*M_PI)*sigmaY1));
}



float kernel2(int i, int j){
    float sigmaX2 = 0.5;
    float sigmaY2 = 4.0;
    float theta2 = 0.0; //rotation of sigmaX axis
    return (exp(-((i*cos(theta2) +j*sin(theta2))*(i*cos(theta2) +j*sin(theta2)))
                    /(2*sigmaX2*sigmaX2)) / (sqrt(2*M_PI)*sigmaX2))
                    *(exp(-((j*cos(theta2) - i*sin(theta2))*(j*cos(theta2) - i*sin(theta2)))
                    /(2*sigmaY2*sigmaY2)) / (sqrt(2*M_PI)*sigmaY2));
}

float kernel3(int i, int j){
    float sigmaX3 = 0.5;
    float sigmaY3 = 4.0;
    float theta3 = 3.14159/4; //rotation of sigmaX axis
    return (exp(-((i*cos(theta3) +j*sin(theta3))*(i*cos(theta3) +j*sin(theta3)))
                    /(2*sigmaX3*sigmaX3)) / (sqrt(2*M_PI)*sigmaX3))
                    *(exp(-((j*cos(theta3) - i*sin(theta3))*(j*cos(theta3) - i*sin(theta3)))
                    /(2*sigmaY3*sigmaY3)) / (sqrt(2*M_PI)*sigmaY3));
}

float kernel4(int i, int j){
    float sigmaX4 = 0.5;
    float sigmaY4 = 4.0;
    float theta4 = 3.14159/2; //rotation of sigmaX axis
    return (exp(-((i*cos(theta4) +j*sin(theta4))*(i*cos(theta4) +j*sin(theta4)))
                    /(2*sigmaX4*sigmaX4)) / (sqrt(2*M_PI)*sigmaX4))
                    *(exp(-((j*cos(theta4) - i*sin(theta4))*(j*cos(theta4) - i*sin(theta4)))
                    /(2*sigmaY4*sigmaY4)) / (sqrt(2*M_PI)*sigmaY4));
}


float kernel5(int i, int j){
    float sigmaX5 = 4.0;
    float sigmaY5 = 4.0;
    float theta5 = 0.0; //rotation of sigmaX axis
    return (exp(-((i*cos(theta5) +j*sin(theta5))*(i*cos(theta5) +j*sin(theta5)))
                    /(2*sigmaX5*sigmaX5)) / (sqrt(2*M_PI)*sigmaX5))
                    *(exp(-((j*cos(theta5) - i*sin(theta5))*(j*cos(theta5) - i*sin(theta5)))
                    /(2*sigmaY5*sigmaY5)) / (sqrt(2*M_PI)*sigmaY5));
}


float kernel(int i, int j, float sigmaX, float sigmaY, float theta){

    return (exp(-((i*cos(theta) +j*sin(theta))*(i*cos(theta) +j*sin(theta)))
                    /(2*sigmaX*sigmaX)) / (sqrt(2*M_PI)*sigmaX))
                    *(exp(-((j*cos(theta) - i*sin(theta))*(j*cos(theta) - i*sin(theta)))
                    /(2*sigmaY*sigmaY)) / (sqrt(2*M_PI)*sigmaY));
}


double kernel1Double(int i, int j){
    double sigmaX1 = 2.0;
    double sigmaY1 = 2.0;
    double theta1 = 0.0; //rotation of sigmaX axis
    return (exp(-((i*cos(theta1) +j*sin(theta1))*(i*cos(theta1) +j*sin(theta1)))
                    /(2*sigmaX1*sigmaX1)) / (sqrt(2*M_PI)*sigmaX1))
                    *(exp(-((j*cos(theta1) - i*sin(theta1))*(j*cos(theta1) - i*sin(theta1)))
                    /(2*sigmaY1*sigmaY1)) / (sqrt(2*M_PI)*sigmaY1));
}



double kernel2Double(int i, int j){
    double sigmaX2 = 0.5;
    double sigmaY2 = 4.0;
    double theta2 = 0.0; //rotation of sigmaX axis
    return (exp(-((i*cos(theta2) +j*sin(theta2))*(i*cos(theta2) +j*sin(theta2)))
                    /(2*sigmaX2*sigmaX2)) / (sqrt(2*M_PI)*sigmaX2))
                    *(exp(-((j*cos(theta2) - i*sin(theta2))*(j*cos(theta2) - i*sin(theta2)))
                    /(2*sigmaY2*sigmaY2)) / (sqrt(2*M_PI)*sigmaY2));
}

double kernel3Double(int i, int j){
    double sigmaX3 = 0.5;
    double sigmaY3 = 4.0;
    double theta3 = 3.14159/4; //rotation of sigmaX axis
    return (exp(-((i*cos(theta3) +j*sin(theta3))*(i*cos(theta3) +j*sin(theta3)))
                    /(2*sigmaX3*sigmaX3)) / (sqrt(2*M_PI)*sigmaX3))
                    *(exp(-((j*cos(theta3) - i*sin(theta3))*(j*cos(theta3) - i*sin(theta3)))
                    /(2*sigmaY3*sigmaY3)) / (sqrt(2*M_PI)*sigmaY3));
}

double kernel4Double(int i, int j){
    double sigmaX4 = 0.5;
    double sigmaY4 = 4.0;
    double theta4 = 3.14159/2; //rotation of sigmaX axis
    return (exp(-((i*cos(theta4) +j*sin(theta4))*(i*cos(theta4) +j*sin(theta4)))
                    /(2*sigmaX4*sigmaX4)) / (sqrt(2*M_PI)*sigmaX4))
                    *(exp(-((j*cos(theta4) - i*sin(theta4))*(j*cos(theta4) - i*sin(theta4)))
                    /(2*sigmaY4*sigmaY4)) / (sqrt(2*M_PI)*sigmaY4));
}


double kernel5Double(int i, int j){
    double sigmaX5 = 4.0;
    double sigmaY5 = 4.0;
    double theta5 = 0.0; //rotation of sigmaX axis
    return (exp(-((i*cos(theta5) +j*sin(theta5))*(i*cos(theta5) +j*sin(theta5)))
                    /(2*sigmaX5*sigmaX5)) / (sqrt(2*M_PI)*sigmaX5))
                    *(exp(-((j*cos(theta5) - i*sin(theta5))*(j*cos(theta5) - i*sin(theta5)))
                    /(2*sigmaY5*sigmaY5)) / (sqrt(2*M_PI)*sigmaY5));
}

double kernelDouble(int i, int j, double sigmaX, double sigmaY, double theta){

    return (exp(-((i*cos(theta) +j*sin(theta))*(i*cos(theta) +j*sin(theta)))
                    /(2*sigmaX*sigmaX)) / (sqrt(2*M_PI)*sigmaX))
                    *(exp(-((j*cos(theta) - i*sin(theta))*(j*cos(theta) - i*sin(theta)))
                    /(2*sigmaY*sigmaY)) / (sqrt(2*M_PI)*sigmaY));
}


int main(int argc, char *argv[]) {
//    for(int kernelSize = 3; kernelSize <=27; kernelSize = kernelSize+2){
        int numberOfBasisKernels = 5;
        int numberOfFuncCoef = 10;
        int kernelSize = 23; //kernel width and height
        int kernelArea = kernelSize*kernelSize;
        int boundingBox = (kernelSize-1)/2;
        
        int imageWidth = 2048;
        int imageHeight = 1489;
        float* outputImg = (float*)malloc(imageWidth*imageHeight*sizeof(float));
        float* outputVar = (float*)malloc(imageWidth*imageHeight*sizeof(float));
        uint16_t* outputMask = (uint16_t*)malloc(imageWidth*imageHeight*sizeof(uint16_t));
    
        printf("kernel size = %d, ", kernelSize);

    #ifdef USE_5_KERNELS 
    //HARD CODED TO 5 BASIS KERNELS, 10 func coefficients
        numberOfBasisKernels = 5;
        #ifdef DOUBLE_PRECISION_KERNEL
            double** kernelArray = (double**)malloc(numberOfBasisKernels*sizeof(double*));
            double (**kernelFuncs)(int, int) = (double (**) (int, int))malloc(5*
                                                    sizeof(double (*)(int, int)));
            kernelFuncs[0] = &kernel1Double;
            kernelFuncs[1] = &kernel2Double;
            kernelFuncs[2] = &kernel3Double;
            kernelFuncs[3] = &kernel4Double;
            kernelFuncs[4] = &kernel5Double;
        #else
            float** kernelArray = (float**)malloc(numberOfBasisKernels*sizeof(float*));
            float (**kernelFuncs)(int, int) = (float (**) (int, int))malloc(5*
                                                    sizeof(float (*)(int, int)));
            kernelFuncs[0] = &kernel1;
            kernelFuncs[1] = &kernel2;
            kernelFuncs[2] = &kernel3;
            kernelFuncs[3] = &kernel4;
            kernelFuncs[4] = &kernel5;
        #endif
    
        int denormalValCount = 0;

        for(int ker = 0; ker < 5; ker++){
            #ifdef DOUBLE_PRECISION_KERNEL
                double * curKernel = (double*)malloc(kernelArea*sizeof(double));
            #else
                float * curKernel = (float*)malloc(kernelArea*sizeof(float));
            #endif
            for(int j = 0; j < kernelSize; j++){
                for(int i = 0; i < kernelSize; i++){
                    if(KERNEL_INFO == "gaussian_contains_denormal_numbers"){
                        #ifdef DOUBLE_PRECISION_KERNEL
                            double curVal = kernelFuncs[ker](i-boundingBox, j-boundingBox);
                            if(curVal < COUNT_AS_DENORMAL){
                                denormalValCount++;
                            }
                            curKernel[j*kernelSize+i] = curVal;
                        #else
                            float curVal = kernelFuncs[ker](i-boundingBox, j-boundingBox);
                            if(curVal < COUNT_AS_DENORMAL){
                                denormalValCount++;
                            }
                            curKernel[j*kernelSize+i] = curVal;
                        #endif
                    }
                    else if(KERNEL_INFO == "gaussian_denormals_zeroed"){
                        #ifdef DOUBLE_PRECISION_KERNEL
                            double curVal = kernelFuncs[ker](i-boundingBox, j-boundingBox);
                            if(curVal > MIN_KERNEL_VAL)
                                curKernel[j*kernelSize+i] = curVal;
                            else
                                curKernel[j*kernelSize+i] = 0.0f;
                        #else
                            float curVal = kernelFuncs[ker](i-boundingBox, j-boundingBox);
                            if(curVal > MIN_KERNEL_VAL)
                                curKernel[j*kernelSize+i] = curVal;
                            else
                                curKernel[j*kernelSize+i] = 0.0f;
                        #endif
                    }
                    else if(KERNEL_INFO == "random"){
                        #ifdef DOUBLE_PRECISION_KERNEL
                            curKernel[j*kernelSize+i] = (double)(rand()) / (double)(RAND_MAX);
                        #else
                            curKernel[j*kernelSize+i] = (float)(rand()) / (float)(RAND_MAX);
                        #endif
                    }
                    else
                        printf("ERROR: invalid KERNEL_INFO");
                }
            }
            kernelArray[ker] = curKernel;
        }
        if(KERNEL_INFO == "gaussian_contains_denormal_numbers"){
            printf("#denormal Kernel values = %d, ", denormalValCount);
        }
    #else
    //FLEXIBLE NUMBER OF BASIS KERNELS
        float** kernelArray = (float**)malloc(numberOfBasisKernels*sizeof(float*));

        int denormalValCount = 0;
        for(int ker = 0; ker < numberOfBasisKernels; ker++){
            float * curKernel = (float*)malloc(kernelSize*kernelSize*sizeof(float));
            for(int j = 0; j < kernelSize; j++){
                for(int i = 0; i < kernelSize; i++){
                    if(KERNEL_INFO == "gaussian_contains_denormal_numbers"){
                        float curVal = curKernel[j*kernelSize+i] = kernel(i-boundingBox , j-boundingBox,
                                                    .12f*(ker+1), .34f*(ker+1), .56f*(ker+1));
                        if(curVal < COUNT_AS_DENORMAL){
                            denormalValCount++;
                        }
                            curKernel[j*kernelSize+i] = curVal;

                    }
                    else if(KERNEL_INFO == "gaussian_denormals_zeroed"){
                        float curVal = curKernel[j*kernelSize+i] = kernel(i-boundingBox , j-boundingBox,
                                                    .12f*(ker+1), .34f*(ker+1), .56f*(ker+1));
                        if(curVal > MIN_KERNEL_VAL)
                            curKernel[j*kernelSize+i] = curVal;
                        else
                            curKernel[j*kernelSize+i] = 0.0f;
                    }
                    else if(KERNEL_INFO == "random"){
                        curKernel[j*kernelSize+i] = (float)(rand()) / (float)(RAND_MAX);
                    }
                    else
                        printf("ERROR: invalid KERNEL_INFO");
    
                }
            }
            kernelArray[ker] = curKernel;
        }
        if(KERNEL_INFO == "gaussian_contains_denormal_numbers"){
            printf("#denormal Kernel values = %d", denormalValCount);
        }
    #endif


    #ifdef TEST_VS_LSST
    
        float* funcParams = (float*)malloc(5*10*sizeof(float));
        for(int i = 0; i < 5; i++){
            for(int j = 0; j < 10; j++){
                funcParams[i*10+j] = 1.0f*(i+1) + .001f*(j+1);
            }
        }
    
        //open the image file containing raw bytes representing floats
        float* inputImg = (float*)malloc(imageWidth*imageHeight*sizeof(float));
        FILE * rawImageFile;
        char imageLocation [100];
        sprintf (imageLocation, "./convolveRawData/inputImg.raw");
        rawImageFile = fopen ( imageLocation , "rb" );
        if (rawImageFile==NULL) {fputs ("File error",stderr); exit (1);}
        // copy the image file into the buffer:
        fread (inputImg,sizeof(float),imageWidth*imageHeight,rawImageFile);
    
    
        //open the variance file containing raw bytes representing floats
        float* inputVar = (float*)malloc(imageWidth*imageHeight*sizeof(float));
        FILE * rawVarianceFile;
        char varianceLocation [100];
        sprintf (varianceLocation, "./convolveRawData/inputVar.raw");
        rawVarianceFile = fopen ( varianceLocation , "rb" );
        if (rawVarianceFile==NULL) {fputs ("File error",stderr); exit (1);}
        // copy the variance file into the buffer:
        fread (inputVar,sizeof(float),imageWidth*imageHeight,rawVarianceFile);
    
        //open the file containing raw bytes representing uint16_t's
        uint16_t* inputMask = (uint16_t*)malloc(imageWidth*imageHeight*sizeof(uint16_t));
        FILE * rawMaskFile;
        char maskLocation [100];
        sprintf (maskLocation, "./convolveRawData/inputMask.raw");
        rawMaskFile = fopen ( maskLocation , "rb" );
        if (rawMaskFile==NULL) {fputs ("File error",stderr); exit (1);}
        // copy the file into the buffer:
        fread (inputMask,sizeof(uint16_t),imageWidth*imageHeight,rawMaskFile);
    
    
        //read in output files
        //open the image file containing raw bytes representing floats
        float* outputReferenceImageArray = (float*)malloc(imageWidth*imageHeight*sizeof(float));
        FILE * outputReferenceImageFile;
        char imageOutputLocation [100];
        sprintf (imageOutputLocation, "./convolveRawData/lsstOutputLinCombo%dx%dImg.raw", kernelSize, kernelSize);
        outputReferenceImageFile = fopen ( imageOutputLocation , "rb" );
        if (outputReferenceImageFile==NULL) {fputs ("File error",stderr); exit (1);}
        // copy the image file into the buffer:
        fread (outputReferenceImageArray,sizeof(float),imageWidth*imageHeight,outputReferenceImageFile);
    
    
        //open the variance file containing raw bytes representing floats
        float* outputReferenceVarianceArray = (float*)malloc(imageWidth*imageHeight*sizeof(float));
        FILE * outputReferenceVarianceFile;
        char varianceOutputLocation [100];
        sprintf (varianceOutputLocation, "./convolveRawData/lsstOutputLinCombo%dx%dVar.raw", kernelSize, kernelSize);
        outputReferenceVarianceFile = fopen ( varianceOutputLocation , "rb" );
        if (outputReferenceVarianceFile==NULL) {fputs ("File error",stderr); exit (1);}
        // copy the variance file into the buffer:
        fread (outputReferenceVarianceArray,sizeof(float),imageWidth*imageHeight,outputReferenceVarianceFile);
    
        //open the mask file containing raw bytes representing uint16_t's
        uint16_t* outputReferenceMaskArray = (uint16_t*)malloc(imageWidth*imageHeight*sizeof(uint16_t));
        FILE * outputReferenceMaskFile;
        char maskOutputLocation [100];
        sprintf (maskOutputLocation, "./convolveRawData/lsstOutputLinCombo%dx%dMask.raw", kernelSize, kernelSize);
        outputReferenceMaskFile = fopen ( maskOutputLocation , "rb" );
        if (outputReferenceMaskFile==NULL) {fputs ("File error",stderr); exit (1);}
        // copy the file into the buffer:
        fread (outputReferenceMaskArray,sizeof(uint16_t),imageWidth*imageHeight,outputReferenceMaskFile);
    
    
    #else
        float* funcParams = (float*)malloc(numberOfBasisKernels*numberOfFuncCoef*sizeof(float));
        for(int i = 0; i < numberOfBasisKernels; i++){
            for(int j = 0; j < numberOfFuncCoef; j++){
                funcParams[i*numberOfFuncCoef+j] = 1.0f*i + .001f*(j+1);
            }
        }
    
    
        float* inputImg = (float*)malloc(imageWidth*imageHeight*sizeof(float));
        float* inputVar = (float*)malloc(imageWidth*imageHeight*sizeof(float));
        uint16_t* inputMask = (uint16_t*)malloc(imageWidth*imageHeight*sizeof(uint16_t));
    
        //set input image values
        for (int y = 0; y < imageHeight; y++) {
            for (int x = 0; x < imageWidth; x++){
                inputImg[y*imageWidth + x] = y*x*cos(x/(y+1));
            }
        }
    
        for(int i = 0; i < imageWidth*imageHeight; i++){
    //        inputImg[i] = (float)(i + 1);
            inputVar[i] = (float)(i + 2);
            inputMask[i] = (uint16_t)(i%imageWidth + i/imageWidth);
    
            outputImg[i] = 0.0f;
            outputVar[i] = 0.0f;
            outputMask[i] = 0;
        }
    #endif
    
            clock_t t1 = clock();
            terraFuncNameInC3(inputImg, inputVar, inputMask, outputImg, outputVar, outputMask,
                            imageWidth, imageHeight, kernelArray, funcParams, numberOfBasisKernels,         
                            kernelSize, kernelSize, numberOfFuncCoef);
        
            clock_t t2 = clock();
            printf("computation took: %f ms, kernelSize = %d\n",  (float)(t2-t1)/1000.0, kernelSize);
        

    
//        if(kernelSize == 3){
//            clock_t t1 = clock();
//        
//            terraFuncNameInC3(inputImg, inputVar, inputMask, outputImg, outputVar, outputMask,
//                            imageWidth, imageHeight, kernelArray, funcParams, numberOfBasisKernels,         
//                            kernelSize, kernelSize, numberOfFuncCoef);
//        
//            clock_t t2 = clock();
//            printf("computation took: %f ms, kernelSize = %d\n",  (float)(t2-t1)/1000.0, kernelSize);
//        }
//
//        else if(kernelSize == 5){
//            clock_t t1 = clock();
//        
//            terraFuncNameInC5(inputImg, inputVar, inputMask, outputImg, outputVar, outputMask,
//                            imageWidth, imageHeight, kernelArray, funcParams, numberOfBasisKernels,         
//                            kernelSize, kernelSize, numberOfFuncCoef);
//        
//            clock_t t2 = clock();
//            printf("computation took: %f ms, kernelSize = %d\n",  (float)(t2-t1)/1000.0, kernelSize);
//        }
//
//        else if(kernelSize == 7){
//            clock_t t1 = clock();
//        
//            terraFuncNameInC7(inputImg, inputVar, inputMask, outputImg, outputVar, outputMask,
//                            imageWidth, imageHeight, kernelArray, funcParams, numberOfBasisKernels,         
//                            kernelSize, kernelSize, numberOfFuncCoef);
//        
//            clock_t t2 = clock();
//            printf("computation took: %f ms, kernelSize = %d\n",  (float)(t2-t1)/1000.0, kernelSize);
//        }
//
//        else if(kernelSize == 9){
//            clock_t t1 = clock();
//        
//            terraFuncNameInC9(inputImg, inputVar, inputMask, outputImg, outputVar, outputMask,
//                            imageWidth, imageHeight, kernelArray, funcParams, numberOfBasisKernels,         
//                            kernelSize, kernelSize, numberOfFuncCoef);
//        
//            clock_t t2 = clock();
//            printf("computation took: %f ms, kernelSize = %d\n",  (float)(t2-t1)/1000.0, kernelSize);
//        }
//
//        else if(kernelSize == 11){
//            clock_t t1 = clock();
//        
//            terraFuncNameInC11(inputImg, inputVar, inputMask, outputImg, outputVar, outputMask,
//                            imageWidth, imageHeight, kernelArray, funcParams, numberOfBasisKernels,         
//                            kernelSize, kernelSize, numberOfFuncCoef);
//        
//            clock_t t2 = clock();
//            printf("computation took: %f ms, kernelSize = %d\n",  (float)(t2-t1)/1000.0, kernelSize);
//        }
//
//        else if(kernelSize == 13){
//            clock_t t1 = clock();
//        
//            terraFuncNameInC13(inputImg, inputVar, inputMask, outputImg, outputVar, outputMask,
//                            imageWidth, imageHeight, kernelArray, funcParams, numberOfBasisKernels,         
//                            kernelSize, kernelSize, numberOfFuncCoef);
//        
//            clock_t t2 = clock();
//            printf("computation took: %f ms, kernelSize = %d\n",  (float)(t2-t1)/1000.0, kernelSize);
//        }
//
//        else if(kernelSize == 15){
//            clock_t t1 = clock();
//        
//            terraFuncNameInC15(inputImg, inputVar, inputMask, outputImg, outputVar, outputMask,
//                            imageWidth, imageHeight, kernelArray, funcParams, numberOfBasisKernels,         
//                            kernelSize, kernelSize, numberOfFuncCoef);
//        
//            clock_t t2 = clock();
//            printf("computation took: %f ms, kernelSize = %d\n",  (float)(t2-t1)/1000.0, kernelSize);
//        }
//
//        else if(kernelSize == 17){
//            clock_t t1 = clock();
//        
//            terraFuncNameInC17(inputImg, inputVar, inputMask, outputImg, outputVar, outputMask,
//                            imageWidth, imageHeight, kernelArray, funcParams, numberOfBasisKernels,         
//                            kernelSize, kernelSize, numberOfFuncCoef);
//        
//            clock_t t2 = clock();
//            printf("computation took: %f ms, kernelSize = %d\n",  (float)(t2-t1)/1000.0, kernelSize);
//        }
//
//        else if(kernelSize == 19){
//            clock_t t1 = clock();
//        
//            terraFuncNameInC19(inputImg, inputVar, inputMask, outputImg, outputVar, outputMask,
//                            imageWidth, imageHeight, kernelArray, funcParams, numberOfBasisKernels,         
//                            kernelSize, kernelSize, numberOfFuncCoef);
//        
//            clock_t t2 = clock();
//            printf("computation took: %f ms, kernelSize = %d\n",  (float)(t2-t1)/1000.0, kernelSize);
//        }
//
//        else if(kernelSize == 21){
//            clock_t t1 = clock();
//        
//            terraFuncNameInC21(inputImg, inputVar, inputMask, outputImg, outputVar, outputMask,
//                            imageWidth, imageHeight, kernelArray, funcParams, numberOfBasisKernels,         
//                            kernelSize, kernelSize, numberOfFuncCoef);
//        
//            clock_t t2 = clock();
//            printf("computation took: %f ms, kernelSize = %d\n",  (float)(t2-t1)/1000.0, kernelSize);
//        }
//
//        else if(kernelSize == 23){
//            clock_t t1 = clock();
//        
//            terraFuncNameInC23(inputImg, inputVar, inputMask, outputImg, outputVar, outputMask,
//                            imageWidth, imageHeight, kernelArray, funcParams, numberOfBasisKernels,         
//                            kernelSize, kernelSize, numberOfFuncCoef);
//        
//            clock_t t2 = clock();
//            printf("computation took: %f ms, kernelSize = %d\n",  (float)(t2-t1)/1000.0, kernelSize);
//        }
//
//        else if(kernelSize == 25){
//            clock_t t1 = clock();
//        
//            terraFuncNameInC25(inputImg, inputVar, inputMask, outputImg, outputVar, outputMask,
//                            imageWidth, imageHeight, kernelArray, funcParams, numberOfBasisKernels,         
//                            kernelSize, kernelSize, numberOfFuncCoef);
//        
//            clock_t t2 = clock();
//            printf("computation took: %f ms, kernelSize = %d\n",  (float)(t2-t1)/1000.0, kernelSize);
//        }
//
//        else if(kernelSize == 27){
//            clock_t t1 = clock();
//        
//            terraFuncNameInC27(inputImg, inputVar, inputMask, outputImg, outputVar, outputMask,
//                            imageWidth, imageHeight, kernelArray, funcParams, numberOfBasisKernels,         
//                            kernelSize, kernelSize, numberOfFuncCoef);
//        
//            clock_t t2 = clock();
//            printf("computation took: %f ms, kernelSize = %d\n",  (float)(t2-t1)/1000.0, kernelSize);
//        }
//
//        else
//            printf("error, invalid kernel size reached");
    
    
    #ifdef TEST_VS_LSST
        float maxImgDif = -1.0f;
        float percentErrorAtMaxImgDif = -1.0f;
        float maxImgPercentError = -1.0f;
        float imgAtMaxPercentError = -1.0f;
        float maxImgDifDividedByVar = -1.0f;
        float maxVarDif = -1.0f;
        float maxVarPercentError = -1.0f;
        uint16_t maxMaskDif = 0;
    
        int markX;
        int markY;
    
        for(int y = boundingBox; y < imageHeight - boundingBox; y++){
            for(int x = boundingBox; x < imageWidth - boundingBox; x++){
                    float curImgDif = fabs(outputImg[y*imageWidth + x] - outputReferenceImageArray[y*imageWidth + x]);
                    float curImgPercentError = curImgDif/outputReferenceImageArray[y*imageWidth + x];
                    float curImgDifDividedByVar = curImgDif/outputReferenceVarianceArray[y*imageWidth + x];
                    float curVarDif = fabs(outputVar[y*imageWidth + x] - outputReferenceVarianceArray[y*imageWidth + x]);
                    float curVarPercentError = curVarDif/outputReferenceVarianceArray[y*imageWidth + x];
                    uint16_t curMaskDif = abs(outputMask[y*imageWidth + x] - outputReferenceMaskArray[y*imageWidth + x]);
    
                if(curImgDif > maxImgDif){
                    maxImgDif = curImgDif;
                    percentErrorAtMaxImgDif = 100*curImgDif/outputReferenceImageArray[y*imageWidth + x];
                }
                if(curImgPercentError > maxImgPercentError){
                    maxImgPercentError = curImgPercentError;
                    imgAtMaxPercentError = outputReferenceImageArray[y*imageWidth + x];
                    markX = x;
                    markY = y;
                }
                if(curImgDifDividedByVar > maxImgDifDividedByVar)
                    maxImgDifDividedByVar = curImgDifDividedByVar;
                if(curVarDif > maxVarDif)
                    maxVarDif = curVarDif;
                if(curVarPercentError > maxVarPercentError)
                    maxVarPercentError = curVarPercentError;
                if(curMaskDif > maxMaskDif)
                    maxMaskDif = curMaskDif;
            }
        }
    
    
        printf("maxImgDif = %f (%f%% error), maxImgPercentError = %f (image = %f), maxImgDifDividedByVar = %f, maxVarDif = %f, maxVarPercentError = %f, maxMaskDif = %d, markX = %d, markY = %d\n\n",
                 maxImgDif, percentErrorAtMaxImgDif, 100*maxImgPercentError, imgAtMaxPercentError, maxImgDifDividedByVar, maxVarDif, 100*maxVarPercentError, maxMaskDif, markX, markY);
    
        printf("%f, %f, %f\n\n", outputImg[boundingBox*imageWidth+boundingBox+1],
            outputReferenceImageArray[boundingBox*imageWidth+boundingBox+1], 
            fabs(outputReferenceImageArray[boundingBox*imageWidth+boundingBox+1] - outputImg[boundingBox*imageWidth+boundingBox+1]));
    
        FILE * imgFileP;
        imgFileP = fopen ("./convolveRawData/TerraLinCombo5x5Img.raw", "wb");
        fwrite (outputImg , sizeof(float), imageWidth*imageHeight, imgFileP);
        fclose (imgFileP);
    
    #ifdef PRINT_CORNERS_OF_OUTPUT_AND_REFERENCE
    
        printf("Image output plane, 10x10 box begining at (boundingBox,boundingBox)\n");
        for(int i=boundingBox; i < boundingBox+10; i++){
            for(int j=boundingBox; j < boundingBox+10; j++){
                printf("%f\t", outputImg[i*imageWidth + j]);
            }
            printf("\n");
        }
        printf("\n");
    
        printf("Variance output plane, 10x10 box begining at (boundingBox,boundingBox)\n");
        for(int i=boundingBox; i < boundingBox+10; i++){
            for(int j=boundingBox; j < boundingBox+10; j++){
                printf("%f\t", outputVar[i*imageWidth + j]);
            }
            printf("\n");
        }
        printf("\n");
    
        printf("Mask output plane, 10x10 box begining at (boundingBox,boundingBox)\n");
        for(int i=boundingBox; i < boundingBox+10; i++){
            for(int j=boundingBox; j < boundingBox+10; j++){
                printf("%d\t", outputMask[i*imageWidth + j]);
            }
            printf("\n");
        }
        printf("\n");
    
        printf("Image reference plane, 10x10 box begining at (boundingBox,boundingBox)\n");
        for(int i=boundingBox; i < boundingBox+10; i++){
            for(int j=boundingBox; j < boundingBox+10; j++){
                printf("%f\t", outputReferenceImageArray[i*imageWidth + j]);
            }
            printf("\n");
        }
        printf("\n");
    
        printf("Variance reference plane, 10x10 box begining at (boundingBox,boundingBox)\n");
        for(int i=boundingBox; i < boundingBox+10; i++){
            for(int j=boundingBox; j < boundingBox+10; j++){
                printf("%f\t", outputReferenceVarianceArray[i*imageWidth + j]);
            }
            printf("\n");
        }
        printf("\n");
    
        printf("Mask reference plane, 10x10 box begining at (boundingBox,boundingBox)\n");
        for(int i=boundingBox; i < boundingBox+10; i++){
            for(int j=boundingBox; j < boundingBox+10; j++){
                printf("%d\t", outputReferenceMaskArray[i*imageWidth + j]);
            }
            printf("\n");
        }
        printf("\n");
    
    #endif
        free(outputReferenceImageArray);
        free(outputReferenceVarianceArray);
        free(outputReferenceMaskArray);
    
    #else
    
    #endif
        for(int ker = 0; ker < numberOfBasisKernels; ker++){
            free(kernelArray[ker]);
        }
        free(inputImg);
        free(inputVar);
        free(inputMask);
    
        free(kernelArray);
        free(funcParams);
        free(outputImg);
        free(outputVar);
        free(outputMask);
//    }
}

