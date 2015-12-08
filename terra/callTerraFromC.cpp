//To Run:
//$ terra blurExampleTerraStandalone.t
//$ g++ -c callTerraFromC.cpp
//$ g++ blurExampleTerraStandalone.o callTerraFromC.o
//$ ./a.out
#include "blurExampleTerraStandalone.h"
#include <math.h>       /* cos */
#include <stdlib.h>

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


int main(int argc, char *argv[]) {
    int numberOfBasisKernels = 5;
    int numberOfFuncCoef = 10;
    int kernelSize = 5; //kernel width and height

    float** kernelArray = (float**)malloc(5*sizeof(float*));

//    float (**kernelFuncs)(int, int) = (float (**) (int, int))malloc(5*
//                                            sizeof(float (*)(int, int)));
//    kernelFuncs[0] = &kernel1;
//    kernelFuncs[1] = &kernel2;
//    kernelFuncs[2] = &kernel3;
//    kernelFuncs[3] = &kernel4;
//    kernelFuncs[4] = &kernel5;
//
//    for(int ker = 0; ker < 5; ker++){
//        float * curKernel = (float*)malloc(25*sizeof(float));
//        for(int j = 0; j < 5; j++){
//            for(int i = 0; i < 5; i++){
//                curKernel[j*5+i] = kernelFuncs[ker](i-2, j-2);
//            }
//        }
//        kernelArray[ker] = curKernel;
//    }
    int boundingBox = (kernelSize - 1)/2;
    for(int ker = 0; ker < numberOfBasisKernels; ker++){
        float * curKernel = (float*)malloc(kernelSize*kernelSize*sizeof(float));
        for(int j = 0; j < kernelSize; j++){
            for(int i = 0; i < kernelSize; i++){
                curKernel[j*kernelSize+i] = kernel(i-boundingBox , j-boundingBox,
                                                .12f*(ker+1), .34f*(ker+1), .56f*(ker+1));
            }
        }
        kernelArray[ker] = curKernel;
    }


    float* funcParams = (float*)malloc(numberOfBasisKernels*numberOfFuncCoef*sizeof(float));
    for(int i = 0; i < numberOfBasisKernels; i++){
        for(int j = 0; j < numberOfFuncCoef; j++){
            funcParams[i*numberOfFuncCoef+j] = 1.0f*i + .001f*(j+1);
        }
    }



    int imageWidth = 2048;
    int imageHeight = 1489;
    float* inputImg = (float*)malloc(imageWidth*imageHeight*sizeof(float));
    float* inputVar = (float*)malloc(imageWidth*imageHeight*sizeof(float));
    uint16_t* inputMask = (uint16_t*)malloc(imageWidth*imageHeight*sizeof(uint16_t));
    float* outputImg = (float*)malloc(imageWidth*imageHeight*sizeof(float));
    float* outputVar = (float*)malloc(imageWidth*imageHeight*sizeof(float));
    uint16_t* outputMask = (uint16_t*)malloc(imageWidth*imageHeight*sizeof(uint16_t));

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


    terraFuncNameInC(inputImg, inputVar, inputMask, outputImg, outputVar, outputMask,
                    imageWidth, imageHeight, kernelArray, funcParams, numberOfBasisKernels,         
                    kernelSize, kernelSize, numberOfFuncCoef);

    for(int ker = 0; ker < numberOfBasisKernels; ker++){
        free(kernelArray[ker]);
    }
    free(kernelArray);
    free(funcParams);
    free(inputImg);
    free(inputVar);
    free(inputMask);
    free(outputImg);
    free(outputVar);
    free(outputMask);
}

