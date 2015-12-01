extern "C"

void terraFuncNameInC(float* inputImg, float* inputVar, unsigned short * inputMask,
					float* outputImg, float* outputVar, unsigned short * outputMask,
                    int imageWidth, int imageHeight, float** kernelArray, 
                    float* funcParams, int numberOfBasisKernels, int kernelWidth,
                    int kernelHeight, int numberOfFuncCoef);
