local C = terralib.includecstring [[
    #include<stdio.h>
    #include<stdlib.h>
    #include<time.h>
    #include<string.h>
    #include<math.h>

]]


local terra polynomial1(x:float, y:float)
 	return 0.1f + 0.002f*x + 0.003f*y + 0.4f*x*x + 0.5f*x*y
                     + 0.6f*y*y +  0.0007f*x*x*x + 0.0008f*x*x*y + 0.0009f*x*y*y
                     + 0.00011f*y*y*y;
end

local terra polynomial2(x:float, y:float)
	return 1.1f + 1.002f*x + 1.003f*y + 1.4f*x*x + 1.5f*x*y
                     + 1.6f*y*y +  1.0007f*x*x*x + 1.0008f*x*x*y + 1.0009f*x*y*y
                     + 1.00011f*y*y*y;
end

local terra polynomial3(x:float, y:float)
	return 2.1f + 2.002f*x + 2.003f*y + 2.4f*x*x + 2.5f*x*y
                     + 2.6f*y*y +  2.0007f*x*x*x + 2.0008f*x*x*y + 2.0009f*x*y*y
                     + 2.00011f*y*y*y;
end

local terra polynomial4(x:float, y:float)
	return 3.1f + 3.002f*x + 3.003f*y + 3.4f*x*x + 3.5f*x*y
                     + 3.6f*y*y +  3.0007f*x*x*x + 3.0008f*x*x*y + 3.0009f*x*y*y
                     + 3.00011f*y*y*y;
end

local terra polynomial5(x:float, y:float)
	return 4.1f + 4.002f*x + 4.003f*y + 4.4f*x*x + 4.5f*x*y
                     + 4.6f*y*y +  4.0007f*x*x*x + 4.0008f*x*x*y + 4.0009f*x*y*y
                     + 4.00011f*y*y*y;
end

--return terra function that evaluates f(x, y), where f is the 3rd degree polynomial
--defined by the array of coefficients coefArray
--local function genPolynomial(coefArray)
--	local terms = {1}
--	local terra polynomial5(x:float, y:float)
--		var returnVal : float = [coefArray[1]]
--		escape
--			for i = 2, 10
--		return 4.1f + 4.002f*x + 4.003f*y + 4.4f*x*x + 4.5f*x*y
--	                     + 4.6f*y*y +  4.0007f*x*x*x + 4.0008f*x*x*y + 4.0009f*x*y*y
--	                     + 4.00011f*y*y*y;
--	end
--end

local function kernel1(i, j)
    local sigmaX1 = 2.0;
    local sigmaY1 = 2.0;
    local theta1 = 0.0; --//rotation of sigmaX axis
    return (math.exp(-((i*math.cos(theta1) +j*math.sin(theta1))*(i*math.cos(theta1) +j*math.sin(theta1)))
                    /(2*sigmaX1*sigmaX1)) / (math.sqrt(2*math.pi)*sigmaX1))
                    *(math.exp(-((j*math.cos(theta1) - i*math.sin(theta1))*(j*math.cos(theta1) - i*math.sin(theta1)))
                    /(2*sigmaY1*sigmaY1)) / (math.sqrt(2*math.pi)*sigmaY1));
end



local function kernel2(i, j)
    local sigmaX2 = 0.5;
    local sigmaY2 = 4.0;
    local theta2 = 0.0; --//rotation of sigmaX axis
    return (math.exp(-((i*math.cos(theta2) +j*math.sin(theta2))*(i*math.cos(theta2) +j*math.sin(theta2)))
                    /(2*sigmaX2*sigmaX2)) / (math.sqrt(2*math.pi)*sigmaX2))
                    *(math.exp(-((j*math.cos(theta2) - i*math.sin(theta2))*(j*math.cos(theta2) - i*math.sin(theta2)))
                    /(2*sigmaY2*sigmaY2)) / (math.sqrt(2*math.pi)*sigmaY2));
end

local function kernel3(i, j)
    local sigmaX3 = 0.5;
    local sigmaY3 = 4.0;
    local theta3 = 3.14159/4; --//rotation of sigmaX axis
    return (math.exp(-((i*math.cos(theta3) +j*math.sin(theta3))*(i*math.cos(theta3) +j*math.sin(theta3)))
                    /(2*sigmaX3*sigmaX3)) / (math.sqrt(2*math.pi)*sigmaX3))
                    *(math.exp(-((j*math.cos(theta3) - i*math.sin(theta3))*(j*math.cos(theta3) - i*math.sin(theta3)))
                    /(2*sigmaY3*sigmaY3)) / (math.sqrt(2*math.pi)*sigmaY3));
end

local function kernel4(i, j)
    local sigmaX4 = 0.5;
    local sigmaY4 = 4.0;
    local theta4 = 3.14159/2; --//rotation of sigmaX axis
    return (math.exp(-((i*math.cos(theta4) +j*math.sin(theta4))*(i*math.cos(theta4) +j*math.sin(theta4)))
                    /(2*sigmaX4*sigmaX4)) / (math.sqrt(2*math.pi)*sigmaX4))
                    *(math.exp(-((j*math.cos(theta4) - i*math.sin(theta4))*(j*math.cos(theta4) - i*math.sin(theta4)))
                    /(2*sigmaY4*sigmaY4)) / (math.sqrt(2*math.pi)*sigmaY4));
end


local function kernel5(i, j)
    local sigmaX5 = 4.0;
    local sigmaY5 = 4.0;
    local theta5 = 0.0; --//rotation of sigmaX axis
    return (math.exp(-((i*math.cos(theta5) +j*math.sin(theta5))*(i*math.cos(theta5) +j*math.sin(theta5)))
                    /(2*sigmaX5*sigmaX5)) / (math.sqrt(2*math.pi)*sigmaX5))
                    *(math.exp(-((j*math.cos(theta5) - i*math.sin(theta5))*(j*math.cos(theta5) - i*math.sin(theta5)))
                    /(2*sigmaY5*sigmaY5)) / (math.sqrt(2*math.pi)*sigmaY5));
end

--[[
local function kernel_func(x,y) -- x and y should be quoted terra
	if x:gettype():isintegral() then
		return `x+1
	else
		return `x+y+1
	end
end
--local kernel = macro(kernel_func)

terra usemacro_example()
	var x = 1
	var y = 3+4
	var z = kernel(x+2,y)
	--var z = [ kernel_func( `x+2, `y ) ]
	return x*y*z
end
--]]

--local vectorWidth = 8
local function zeroFloatVec(vectorWidth)
	local zeros = {}
	for i=1,vectorWidth do zeros[i] = `0.f end
	return `vector(zeros)
end

local function zeroUInt16Vec(vectorWidth)
	local zeros = {}
	for i=1,vectorWidth do zeros[i] = `0 end
	return `vector(zeros)
end


local function genGetBitMask(vectorWidth)
	local terra returnFunc(a : vector(float,vectorWidth), b : vector(float,vectorWidth))
	    return terralib.select(a ~= b, [vector(uint16,vectorWidth)](0xFFFFULL),[vector(uint16,vectorWidth)](0) ) 
	end
	return returnFunc
end


local function metaProgramBlurMaskedImageLinComboVectorize(boundingBox, vectorWidth, method, numberOfBasisKernels)
	local getBitMask = genGetBitMask(vectorWidth) --getBitMask is a terra function
	local terra blurMaskedImageLinComboVectorize()
	--	var boundingBox : int = 2
		var width : int = 2048
		var height : int = 1489

		var inputIm : &float
		inputIm = [&float](C.malloc(sizeof(float)*width*height))
		var inputVar : &float
		inputVar = [&float](C.malloc(sizeof(float)*width*height))
		var inputMask : &uint16
		inputMask = [&uint16](C.malloc(sizeof(uint16)*width*height))

		var outputIm : &float
		outputIm = [&float](C.malloc(sizeof(float)*width*height))
		var outputVar : &float
		outputVar = [&float](C.malloc(sizeof(float)*width*height))
		var outputMask : &uint16
		outputMask = [&uint16](C.malloc(sizeof(uint16)*width*height))

		for x = 0, width do
			for y = 0, height do
				inputIm[y*width + x] = [float](y*x*C.cos(x/(y+1)))
				inputVar[y*width + x] = [float](x*y + 2)
				inputMask[y*width + x] = [uint16](x + y)


				outputIm[y*width + x] = 0.0f
				outputVar[y*width + x] = 0.0f
				outputMask[y*width + x] = 0
			end
		end


		var t1 = C.clock()


		var zeroVec : vector(float, vectorWidth) = [zeroFloatVec(vectorWidth)]
		for y = boundingBox, height-boundingBox do
			for x = boundingBox, (width-boundingBox), vectorWidth do
 				escape
 					if(method == "original") then
 						local curKernelVal = symbol()
 						local curKernelValTemp = symbol()
 						local curImOut = symbol()
 						local curVarOut= symbol()
 						local curMaskOut = symbol()
 						local curNorm = symbol()
	 					emit quote
			 				var [curKernelVal] : vector(float,vectorWidth)
							var [curKernelValTemp] : float[vectorWidth]

							var [curImOut] : vector(float,vectorWidth) = [zeroFloatVec(vectorWidth)]
			--				var curImOut : &vector(float,vectorWidth) = [&vector(float,vectorWidth)](&outputIm[y*width + x])
							var [curVarOut] : vector(float,vectorWidth) = [zeroFloatVec(vectorWidth)]
							var [curMaskOut] : vector(uint16, vectorWidth) = [zeroUInt16Vec(vectorWidth)]

							var [curNorm] : vector(float,vectorWidth) = [zeroFloatVec(vectorWidth)]
			--				var curImOutVec = @[&vector(float,4)](&outputIm[y*width + x])


--						 	var curKernelVal : vector(float,vectorWidth)
--							var curKernelValTemp : float[vectorWidth]
--							var curImOut : vector(float,vectorWidth) = [zeroFloatVec(vectorWidth)]
--							var curVarOut : vector(float,vectorWidth) = [zeroFloatVec(vectorWidth)]
--							var curMaskOut : vector(uint16, vectorWidth) = [zeroUInt16Vec(vectorWidth)]
--							var curNorm : vector(float,vectorWidth) = [zeroFloatVec(vectorWidth)]

						end
						for j = -boundingBox, boundingBox do
					        for i = -boundingBox, boundingBox do
					        	for k = 0, vectorWidth - 1 do
					        		if(numberOfBasisKernels==5) then
					        			emit quote
							        		curKernelValTemp[k] = polynomial1(x+k, y)*[kernel1(i, j)] +
								                polynomial2(x+k, y)*[kernel2(i, j)] + polynomial3(x+k, y)*[kernel3(i, j)] + 
								                polynomial4(x+k, y)*[kernel4(i, j)] + polynomial5(x+k, y)*[kernel5(i, j)];
						        		end
						        	elseif(numberOfBasisKernels==4) then
						        		emit quote
							        		curKernelValTemp[k] = polynomial1(x+k, y)*[kernel1(i, j)] +
								                polynomial2(x+k, y)*[kernel2(i, j)] + polynomial3(x+k, y)*[kernel3(i, j)] + 
								                polynomial4(x+k, y)*[kernel4(i, j)];
						        		end
						        	elseif(numberOfBasisKernels==3) then
						        		emit quote
							        		curKernelValTemp[k] = polynomial1(x+k, y)*[kernel1(i, j)] +
								                polynomial2(x+k, y)*[kernel2(i, j)] + polynomial3(x+k, y)*[kernel3(i, j)] 
						        		end
						        	elseif(numberOfBasisKernels==2) then
						        		emit quote
							        		curKernelValTemp[k] = polynomial1(x+k, y)*[kernel1(i, j)] +
								                polynomial2(x+k, y)*[kernel2(i, j)] 
						        		end
						        	else
						        		emit quote
							        		curKernelValTemp[k] = polynomial1(x+k, y)*[kernel1(i, j)]
						        		end
						        	end

					        	end
					        	emit quote
						            curKernelVal = @[&vector(float,vectorWidth)](&curKernelValTemp[0])

						            var curImInVec = @[&vector(float,vectorWidth)](&inputIm[(y+j)*width + x+i])
						            var curVarInVec = @[&vector(float,vectorWidth)](&inputVar[(y+j)*width + x+i])
						            var curMaskInVec = @[&vector(uint16,vectorWidth)](&inputMask[(y+j)*width + x+i])


						            curImOut = curImOut + curImInVec*curKernelVal; 
						            curVarOut = curVarOut + curVarInVec*curKernelVal*curKernelVal;

						            var bitMask : vector(uint16, vectorWidth) = getBitMask(curKernelVal, zeroVec)
						            curMaskOut = curMaskOut or (curMaskInVec and bitMask)

						            curNorm = curNorm + curKernelVal;
						        end
					        end
					    end
					    emit quote
				    		curImOut = curImOut/curNorm
				    		curVarOut = curVarOut/(curNorm*curNorm)
						end
						for i = 0, vectorWidth-1 do
				    		emit quote
					    		outputIm[y*width + x+i] = curImOut[i]	
							    outputVar[y*width + x+i] = curVarOut[i]
				    			outputMask[y*width + x+i] = curMaskOut[i]
				    		end
				    	end
				    end				
				end
			end
		end

		var t2 = C.clock()
--		C.printf("\n\nVectorized %d wide masked image lin combo %dx%d blur Terra:\n", vectorWidth, 1+2*boundingBox, 1+2*boundingBox)
--		C.printf("outputIm[boundingBox*width + boundingBox] = %f, other planes = %f, %f, computation took: %f ms\n",
--			outputIm[boundingBox*width + boundingBox], outputVar[boundingBox*width + boundingBox],
--			outputMask[boundingBox*width + boundingBox],  (t2-t1)/1000.0)
--		C.printf("outputIm[boundingBox*width + boundingBox] = %f, computation took: %f ms\n",
--			outputIm[boundingBox*width + boundingBox], (t2-t1)/1000.0)
		C.printf("outputIm[boundingBox*width + boundingBox] = %f, computation took: %f ms; %f, %d\n",
			outputIm[boundingBox*width + boundingBox], (t2-t1)/1000.0, outputVar[boundingBox*width + boundingBox], outputMask[boundingBox*width + boundingBox])

--		C.printf("C.CLOCKS_PER_SEC = %d\n", C.CLOCKS_PER_SEC)
--
--		C.printf("Image plane, 10x10 box begining at (boundingBox,boundingBox)\n")
--		for i=boundingBox,boundingBox+10 do
--			for j=boundingBox,boundingBox+10 do
--				C.printf("%f\t", outputIm[i*width + j])
--			end
--			C.printf("\n")
--		end
--
--		C.printf("Variance plane, 10x10 box begining at (boundingBox,boundingBox)\n")
--		for i=boundingBox,boundingBox+10 do
--			for j=boundingBox,boundingBox+10 do
--				C.printf("%f\t", outputVar[i*width + j])
--			end
--			C.printf("\n")
--		end
--
--		C.printf("Mask plane, 10x10 box begining at (boundingBox,boundingBox)\n")
--		for i=boundingBox,boundingBox+10 do
--			for j=boundingBox,boundingBox+10 do
--				C.printf("%d\t", outputMask[i*width + j])
--			end
--			C.printf("\n")
--		end

	    C.free(inputIm)
		C.free(outputIm)
	    C.free(inputVar)
		C.free(outputVar)
	    C.free(inputMask)
		C.free(outputMask)
	end
	return blurMaskedImageLinComboVectorize
end



local function metaProgramBlurMaskedImageLinComboVectorizeUnrollFlexiblePolyCoef(boundingBox, vectorWidth, method)
	local getBitMask = genGetBitMask(vectorWidth) --getBitMask is a terra function
	local terra blurMaskedImageLinComboVectorize()
	--	var boundingBox : int = 2
		var width : int = 2048
		var height : int = 1489

		var inputIm : &float
		inputIm = [&float](C.malloc(sizeof(float)*width*height))
		var inputVar : &float
		inputVar = [&float](C.malloc(sizeof(float)*width*height))
		var inputMask : &uint16
		inputMask = [&uint16](C.malloc(sizeof(uint16)*width*height))

		var outputIm : &float
		outputIm = [&float](C.malloc(sizeof(float)*width*height))
		var outputVar : &float
		outputVar = [&float](C.malloc(sizeof(float)*width*height))
		var outputMask : &uint16
		outputMask = [&uint16](C.malloc(sizeof(uint16)*width*height))

		for x = 0, width do
			for y = 0, height do
				inputIm[y*width + x] = [float](y*x*C.cos(x/(y+1)))
				inputVar[y*width + x] = [float](x*y + 2)
				inputMask[y*width + x] = [uint16](x + y)


				outputIm[y*width + x] = 0.0f
				outputVar[y*width + x] = 0.0f
				outputMask[y*width + x] = 0
			end
		end


		var t1 = C.clock()


		var zeroVec : vector(float, vectorWidth) = [zeroFloatVec(vectorWidth)]
		for y = boundingBox, height-boundingBox do
			for x = boundingBox, (width-boundingBox), vectorWidth do
 				escape
 					if(method == "original") then
 						local curKernelVal = symbol()
 						local curKernelValTemp = symbol()
 						local curImOut = symbol()
 						local curVarOut= symbol()
 						local curMaskOut = symbol()
 						local curNorm = symbol()
 						local curPolyVals = symbol()
	 					emit quote
			 				var [curKernelVal] : vector(float,vectorWidth)
							var [curKernelValTemp] : float[vectorWidth]

							var [curImOut] : vector(float,vectorWidth) = [zeroFloatVec(vectorWidth)]
			--				var curImOut : &vector(float,vectorWidth) = [&vector(float,vectorWidth)](&outputIm[y*width + x])
							var [curVarOut] : vector(float,vectorWidth) = [zeroFloatVec(vectorWidth)]
							var [curMaskOut] : vector(uint16, vectorWidth) = [zeroUInt16Vec(vectorWidth)]

							var [curNorm] : vector(float,vectorWidth) = [zeroFloatVec(vectorWidth)]
			--				var curImOutVec = @[&vector(float,4)](&outputIm[y*width + x])


--						 	var curKernelVal : vector(float,vectorWidth)
--							var curKernelValTemp : float[vectorWidth]
--							var curImOut : vector(float,vectorWidth) = [zeroFloatVec(vectorWidth)]
--							var curVarOut : vector(float,vectorWidth) = [zeroFloatVec(vectorWidth)]
--							var curMaskOut : vector(uint16, vectorWidth) = [zeroUInt16Vec(vectorWidth)]
--							var curNorm : vector(float,vectorWidth) = [zeroFloatVec(vectorWidth)]

							var [curPolyVals] : float[vectorWidth*5] --curPolyVals[5*i + j] contains
																	--the jth polynomial value for the ith vector position
							for i = 0, vectorWidth do
								curPolyVals[5*i] = polynomial1(x+i, y)
								curPolyVals[5*i+1] = polynomial2(x+i, y)
								curPolyVals[5*i+2] = polynomial3(x+i, y)
								curPolyVals[5*i+3] = polynomial4(x+i, y)
								curPolyVals[5*i+4] = polynomial5(x+i, y)
							end

						end
						for j = -boundingBox, boundingBox do
					        for i = -boundingBox, boundingBox do
					        	for k = 0, vectorWidth - 1 do
					        		emit quote
						        		curKernelValTemp[k] = curPolyVals[5*k]*[kernel1(i, j)] +
							                curPolyVals[5*k+1]*[kernel2(i, j)] + curPolyVals[5*k+2]*[kernel3(i, j)] + 
							                curPolyVals[5*k+3]*[kernel4(i, j)] + curPolyVals[5*k+4]*[kernel5(i, j)];
					        		end
					        	end
					        	emit quote
						            curKernelVal = @[&vector(float,vectorWidth)](&curKernelValTemp[0])

						            var curImInVec = @[&vector(float,vectorWidth)](&inputIm[(y+j)*width + x+i])
						            var curVarInVec = @[&vector(float,vectorWidth)](&inputVar[(y+j)*width + x+i])
						            var curMaskInVec = @[&vector(uint16,vectorWidth)](&inputMask[(y+j)*width + x+i])


						            curImOut = curImOut + curImInVec*curKernelVal; 
						            curVarOut = curVarOut + curVarInVec*curKernelVal*curKernelVal;

						            var bitMask : vector(uint16, vectorWidth) = getBitMask(curKernelVal, zeroVec)
						            curMaskOut = curMaskOut or (curMaskInVec and bitMask)

						            curNorm = curNorm + curKernelVal;
						        end
					        end
					    end
					    emit quote
				    		curImOut = curImOut/curNorm
				    		curVarOut = curVarOut/(curNorm*curNorm)
						end
						for i = 0, vectorWidth-1 do
				    		emit quote
					    		outputIm[y*width + x+i] = curImOut[i]	
							    outputVar[y*width + x+i] = curVarOut[i]
				    			outputMask[y*width + x+i] = curMaskOut[i]
				    		end
				    	end
				    end				
				end
			end
		end

		var t2 = C.clock()
--		C.printf("\n\nVectorized %d wide masked image lin combo %dx%d blur Terra:\n", vectorWidth, 1+2*boundingBox, 1+2*boundingBox)
--		C.printf("outputIm[boundingBox*width + boundingBox] = %f, other planes = %f, %f, computation took: %f ms\n",
--			outputIm[boundingBox*width + boundingBox], outputVar[boundingBox*width + boundingBox],
--			outputMask[boundingBox*width + boundingBox],  (t2-t1)/1000.0)
--		C.printf("outputIm[boundingBox*width + boundingBox] = %f, computation took: %f ms\n",
--			outputIm[boundingBox*width + boundingBox], (t2-t1)/1000.0)
		C.printf("outputIm[boundingBox*width + boundingBox] = %f, computation took: %f ms; %f, %d\n",
			outputIm[boundingBox*width + boundingBox], (t2-t1)/1000.0, outputVar[boundingBox*width + boundingBox], outputMask[boundingBox*width + boundingBox])

--		C.printf("C.CLOCKS_PER_SEC = %d\n", C.CLOCKS_PER_SEC)
--
--		C.printf("Image plane, 10x10 box begining at (boundingBox,boundingBox)\n")
--		for i=boundingBox,boundingBox+10 do
--			for j=boundingBox,boundingBox+10 do
--				C.printf("%f\t", outputIm[i*width + j])
--			end
--			C.printf("\n")
--		end
--
--		C.printf("Variance plane, 10x10 box begining at (boundingBox,boundingBox)\n")
--		for i=boundingBox,boundingBox+10 do
--			for j=boundingBox,boundingBox+10 do
--				C.printf("%f\t", outputVar[i*width + j])
--			end
--			C.printf("\n")
--		end
--
--		C.printf("Mask plane, 10x10 box begining at (boundingBox,boundingBox)\n")
--		for i=boundingBox,boundingBox+10 do
--			for j=boundingBox,boundingBox+10 do
--				C.printf("%d\t", outputMask[i*width + j])
--			end
--			C.printf("\n")
--		end

	    C.free(inputIm)
		C.free(outputIm)
	    C.free(inputVar)
		C.free(outputVar)
	    C.free(inputMask)
		C.free(outputMask)
	end
	return blurMaskedImageLinComboVectorize
end


local kernelVals0 = symbol()
local kernelArea0 = symbol()
local kernelWidth0 = symbol()

local function metaProgramBlurMaskedImageLinComboVectorizeUnrollFlexiblePolyAndKernel(boundingBox, vectorWidth, method)
	local getBitMask = genGetBitMask(vectorWidth) --getBitMask is a terra function
	local terra blurMaskedImageLinComboVectorize()
	--	var boundingBox : int = 2
		var width : int = 2048
		var height : int = 1489

		var inputIm : &float
		inputIm = [&float](C.malloc(sizeof(float)*width*height))
		var inputVar : &float
		inputVar = [&float](C.malloc(sizeof(float)*width*height))
		var inputMask : &uint16
		inputMask = [&uint16](C.malloc(sizeof(uint16)*width*height))

		var outputIm : &float
		outputIm = [&float](C.malloc(sizeof(float)*width*height))
		var outputVar : &float
		outputVar = [&float](C.malloc(sizeof(float)*width*height))
		var outputMask : &uint16
		outputMask = [&uint16](C.malloc(sizeof(uint16)*width*height))

		for x = 0, width do
			for y = 0, height do
				inputIm[y*width + x] = [float](y*x*C.cos(x/(y+1)))
				inputVar[y*width + x] = [float](x*y + 2)
				inputMask[y*width + x] = [uint16](x + y)


				outputIm[y*width + x] = 0.0f
				outputVar[y*width + x] = 0.0f
				outputMask[y*width + x] = 0
			end
		end


		var t1 = C.clock()


		var [kernelArea0] = (boundingBox*2+1)*(boundingBox*2+1)
		var [kernelWidth0] = boundingBox*2+1

		var [kernelVals0] : &float = [&float](C.malloc(kernelArea0*5*sizeof(float)))--float[kernelArea0*5]
	--							var [kernelVals0] : float[121*5]
		var tA = C.clock()
		escape
			for j = -boundingBox, boundingBox do
				for i = -boundingBox, boundingBox do
					emit quote
						kernelVals0[kernelArea0*0 + kernelWidth0*(j+boundingBox) + (i+boundingBox)] = [kernel1(i, j)]
						kernelVals0[kernelArea0*1 + kernelWidth0*(j+boundingBox) + (i+boundingBox)] = [kernel2(i, j)]
						kernelVals0[kernelArea0*2 + kernelWidth0*(j+boundingBox) + (i+boundingBox)] = [kernel3(i, j)]
						kernelVals0[kernelArea0*3 + kernelWidth0*(j+boundingBox) + (i+boundingBox)] = [kernel4(i, j)]
						kernelVals0[kernelArea0*4 + kernelWidth0*(j+boundingBox) + (i+boundingBox)] = [kernel5(i, j)]
					end
				end
			end
		end
		var tB = C.clock()
		C.printf("calculating kernel values took: %fms", (tB-tA)/C.CLOCKS_PER_SEC*1000)

		var zeroVec : vector(float, vectorWidth) = [zeroFloatVec(vectorWidth)]
		for y = boundingBox, height-boundingBox do
			for x = boundingBox, (width-boundingBox), vectorWidth do
 				escape
 					if(method == "original") then
 						local curKernelVal = symbol()
 						local curKernelValTemp = symbol()
 						local curImOut = symbol()
 						local curVarOut= symbol()
 						local curMaskOut = symbol()
 						local curNorm = symbol()

 						local curPolyVals = symbol()

	 					emit quote
			 				var [curKernelVal] : vector(float,vectorWidth)
							var [curKernelValTemp] : float[vectorWidth]

							var [curImOut] : vector(float,vectorWidth) = [zeroFloatVec(vectorWidth)]
			--				var curImOut : &vector(float,vectorWidth) = [&vector(float,vectorWidth)](&outputIm[y*width + x])
							var [curVarOut] : vector(float,vectorWidth) = [zeroFloatVec(vectorWidth)]
							var [curMaskOut] : vector(uint16, vectorWidth) = [zeroUInt16Vec(vectorWidth)]

							var [curNorm] : vector(float,vectorWidth) = [zeroFloatVec(vectorWidth)]
			--				var curImOutVec = @[&vector(float,4)](&outputIm[y*width + x])
						end
						emit quote
							--HARD CODED TO 5 BASIS KERNELS!!!!!!!!!!!!!!!

							var [curPolyVals] : float[vectorWidth*5] --curPolyVals[5*i + j] contains
																	--the jth polynomial value for the ith vector position
							for i = 0, vectorWidth do
								curPolyVals[5*i] = polynomial1(x+i, y)
								curPolyVals[5*i+1] = polynomial2(x+i, y)
								curPolyVals[5*i+2] = polynomial3(x+i, y)
								curPolyVals[5*i+3] = polynomial4(x+i, y)
								curPolyVals[5*i+4] = polynomial5(x+i, y)
							end


							escape
								for j = -boundingBox, boundingBox do
						        	for i = -boundingBox, boundingBox do
							--for j = -boundingBox, boundingBox + 1 do
						    --    for i = -boundingBox, boundingBox + 1 do
							---        escape
							        	for k = 0, vectorWidth - 1 do
							        		emit quote
								        		curKernelValTemp[k] = curPolyVals[5*k]*kernelVals0[kernelArea0*0 + kernelWidth0*(j+boundingBox) + (i+boundingBox)] +
									                curPolyVals[5*k+1]*kernelVals0[kernelArea0*1 + kernelWidth0*(j+boundingBox) + (i+boundingBox)] + curPolyVals[5*k+2]*kernelVals0[kernelArea0*2 + kernelWidth0*(j+boundingBox) + (i+boundingBox)] + 
									                curPolyVals[5*k+3]*kernelVals0[kernelArea0*3 + kernelWidth0*(j+boundingBox) + (i+boundingBox)] + curPolyVals[5*k+4]*kernelVals0[kernelArea0*4 + kernelWidth0*(j+boundingBox) + (i+boundingBox)];


								        		--curKernelValTemp[k] = polynomial1(x+k, y)*[kernel1(i, j)] +
									            --    polynomial2(x+k, y)*[kernel2(i, j)] + polynomial3(x+k, y)*[kernel3(i, j)] + 
									            --    polynomial4(x+k, y)*[kernel4(i, j)] + polynomial5(x+k, y)*[kernel5(i, j)];
							        		end
							        	end
							        	emit quote
								            curKernelVal = @[&vector(float,vectorWidth)](&curKernelValTemp[0])

								            var curImInVec = @[&vector(float,vectorWidth)](&inputIm[(y+j)*width + x+i])
								            var curVarInVec = @[&vector(float,vectorWidth)](&inputVar[(y+j)*width + x+i])
								            var curMaskInVec = @[&vector(uint16,vectorWidth)](&inputMask[(y+j)*width + x+i])


								            curImOut = curImOut + curImInVec*curKernelVal; 
								            curVarOut = curVarOut + curVarInVec*curKernelVal*curKernelVal;

								            var bitMask : vector(uint16, vectorWidth) = getBitMask(curKernelVal, zeroVec)
								            curMaskOut = curMaskOut or (curMaskInVec and bitMask)

								            curNorm = curNorm + curKernelVal;
								        end
							        end
						    	end
						    end
							--end
						end
					    emit quote
				    		curImOut = curImOut/curNorm
				    		curVarOut = curVarOut/(curNorm*curNorm)
						end
						for i = 0, vectorWidth-1 do
				    		emit quote
					    		outputIm[y*width + x+i] = curImOut[i]	
							    outputVar[y*width + x+i] = curVarOut[i]
				    			outputMask[y*width + x+i] = curMaskOut[i]
				    		end
				    	end
				    end				
				end
			end
		end

		var t2 = C.clock()

		C.printf("outputIm[boundingBox*width + boundingBox] = %f, computation took: %f ms; %f, %d\n",
			outputIm[boundingBox*width + boundingBox], (t2-t1)/1000.0, outputVar[boundingBox*width + boundingBox], outputMask[boundingBox*width + boundingBox])


--		C.printf("\n\nVectorized %d wide masked image lin combo %dx%d blur Terra NO UNROLL:\n", vectorWidth, 1+2*boundingBox, 1+2*boundingBox)
--		C.printf("outputIm[boundingBox*width + boundingBox] = %f, computation took: %f ms\n", outputIm[boundingBox*width + boundingBox],  (t2-t1)/1000.0)
--		C.printf("C.CLOCKS_PER_SEC = %d\n", C.CLOCKS_PER_SEC)
--
--		C.printf("Image plane, 10x10 box begining at (boundingBox,boundingBox)\n")
--		for i=boundingBox,boundingBox+10 do
--			for j=boundingBox,boundingBox+10 do
--				C.printf("%f\t", outputIm[i*width + j])
--			end
--			C.printf("\n")
--		end
--
--		C.printf("Variance plane, 10x10 box begining at (boundingBox,boundingBox)\n")
--		for i=boundingBox,boundingBox+10 do
--			for j=boundingBox,boundingBox+10 do
--				C.printf("%f\t", outputVar[i*width + j])
--			end
--			C.printf("\n")
--		end
--
--		C.printf("Mask plane, 10x10 box begining at (boundingBox,boundingBox)\n")
--		for i=boundingBox,boundingBox+10 do
--			for j=boundingBox,boundingBox+10 do
--				C.printf("%d\t", outputMask[i*width + j])
--			end
--			C.printf("\n")
--		end

	    C.free(inputIm)
		C.free(outputIm)
	    C.free(inputVar)
		C.free(outputVar)
	    C.free(inputMask)
		C.free(outputMask)
		C.free([kernelVals0])

	end
	return blurMaskedImageLinComboVectorize
end


local kernelVals = symbol()
local kernelArea = symbol()
local kernelWidth = symbol()

local function metaProgramBlurMaskedImageLinComboVectorizeNoUnroll(boundingBox, vectorWidth, method)
	local getBitMask = genGetBitMask(vectorWidth) --getBitMask is a terra function
	local terra blurMaskedImageLinComboVectorize()
	--	var boundingBox : int = 2
		var width : int = 2048
		var height : int = 1489

		var inputIm : &float
		inputIm = [&float](C.malloc(sizeof(float)*width*height))
		var inputVar : &float
		inputVar = [&float](C.malloc(sizeof(float)*width*height))
		var inputMask : &uint16
		inputMask = [&uint16](C.malloc(sizeof(uint16)*width*height))

		var outputIm : &float
		outputIm = [&float](C.malloc(sizeof(float)*width*height))
		var outputVar : &float
		outputVar = [&float](C.malloc(sizeof(float)*width*height))
		var outputMask : &uint16
		outputMask = [&uint16](C.malloc(sizeof(uint16)*width*height))

		for x = 0, width do
			for y = 0, height do
				inputIm[y*width + x] = [float](y*x*C.cos(x/(y+1)))
				inputVar[y*width + x] = [float](x*y + 2)
				inputMask[y*width + x] = [uint16](x + y)


				outputIm[y*width + x] = 0.0f
				outputVar[y*width + x] = 0.0f
				outputMask[y*width + x] = 0
			end
		end


		var t1 = C.clock()

		var [kernelArea] = (boundingBox*2+1)*(boundingBox*2+1)
		var [kernelWidth] = boundingBox*2+1

		var [kernelVals] : &float = [&float](C.malloc(kernelArea*5*sizeof(float)))--float[kernelArea*5]
	--							var [kernelVals] : float[121*5]

		escape
			for j = -boundingBox, boundingBox do
				for i = -boundingBox, boundingBox do
					emit quote
						kernelVals[kernelArea*0 + kernelWidth*(j+boundingBox) + (i+boundingBox)] = [kernel1(i, j)]
						kernelVals[kernelArea*1 + kernelWidth*(j+boundingBox) + (i+boundingBox)] = [kernel2(i, j)]
						kernelVals[kernelArea*2 + kernelWidth*(j+boundingBox) + (i+boundingBox)] = [kernel3(i, j)]
						kernelVals[kernelArea*3 + kernelWidth*(j+boundingBox) + (i+boundingBox)] = [kernel4(i, j)]
						kernelVals[kernelArea*4 + kernelWidth*(j+boundingBox) + (i+boundingBox)] = [kernel5(i, j)]
					end
				end
			end
		end

		var zeroVec : vector(float, vectorWidth) = [zeroFloatVec(vectorWidth)]
		for y = boundingBox, height-boundingBox do
			for x = boundingBox, (width-boundingBox), vectorWidth do
 				escape
 					if(method == "original") then
 						local curKernelVal = symbol()
 						local curKernelValTemp = symbol()
 						local curImOut = symbol()
 						local curVarOut= symbol()
 						local curMaskOut = symbol()
 						local curNorm = symbol()

 						local curPolyVals = symbol()

	 					emit quote
			 				var [curKernelVal] : vector(float,vectorWidth)
							var [curKernelValTemp] : float[vectorWidth]

							var [curImOut] : vector(float,vectorWidth) = [zeroFloatVec(vectorWidth)]
			--				var curImOut : &vector(float,vectorWidth) = [&vector(float,vectorWidth)](&outputIm[y*width + x])
							var [curVarOut] : vector(float,vectorWidth) = [zeroFloatVec(vectorWidth)]
							var [curMaskOut] : vector(uint16, vectorWidth) = [zeroUInt16Vec(vectorWidth)]

							var [curNorm] : vector(float,vectorWidth) = [zeroFloatVec(vectorWidth)]
			--				var curImOutVec = @[&vector(float,4)](&outputIm[y*width + x])
						end
						emit quote
							--HARD CODED TO 5 BASIS KERNELS!!!!!!!!!!!!!!!

							var [curPolyVals] : float[vectorWidth*5] --curPolyVals[5*i + j] contains
																	--the jth polynomial value for the ith vector position
							for i = 0, vectorWidth do
								curPolyVals[5*i] = polynomial1(x+i, y)
								curPolyVals[5*i+1] = polynomial2(x+i, y)
								curPolyVals[5*i+2] = polynomial3(x+i, y)
								curPolyVals[5*i+3] = polynomial4(x+i, y)
								curPolyVals[5*i+4] = polynomial5(x+i, y)
							end


							--escape
							--for j = -boundingBox, boundingBox do
						    --    for i = -boundingBox, boundingBox do
							for j = -boundingBox, boundingBox + 1 do
						        for i = -boundingBox, boundingBox + 1 do
							        escape
							        	for k = 0, vectorWidth - 1 do
							        		emit quote
								        		curKernelValTemp[k] = curPolyVals[5*k]*kernelVals[kernelArea*0 + kernelWidth*(j+boundingBox) + (i+boundingBox)] +
									                curPolyVals[5*k+1]*kernelVals[kernelArea*1 + kernelWidth*(j+boundingBox) + (i+boundingBox)] + curPolyVals[5*k+2]*kernelVals[kernelArea*2 + kernelWidth*(j+boundingBox) + (i+boundingBox)] + 
									                curPolyVals[5*k+3]*kernelVals[kernelArea*3 + kernelWidth*(j+boundingBox) + (i+boundingBox)] + curPolyVals[5*k+4]*kernelVals[kernelArea*4 + kernelWidth*(j+boundingBox) + (i+boundingBox)];


								        		--curKernelValTemp[k] = polynomial1(x+k, y)*[kernel1(i, j)] +
									            --    polynomial2(x+k, y)*[kernel2(i, j)] + polynomial3(x+k, y)*[kernel3(i, j)] + 
									            --    polynomial4(x+k, y)*[kernel4(i, j)] + polynomial5(x+k, y)*[kernel5(i, j)];
							        		end
							        	end
							        	emit quote
								            curKernelVal = @[&vector(float,vectorWidth)](&curKernelValTemp[0])

								            var curImInVec = @[&vector(float,vectorWidth)](&inputIm[(y+j)*width + x+i])
								            var curVarInVec = @[&vector(float,vectorWidth)](&inputVar[(y+j)*width + x+i])
								            var curMaskInVec = @[&vector(uint16,vectorWidth)](&inputMask[(y+j)*width + x+i])


								            curImOut = curImOut + curImInVec*curKernelVal; 
								            curVarOut = curVarOut + curVarInVec*curKernelVal*curKernelVal;

								            var bitMask : vector(uint16, vectorWidth) = getBitMask(curKernelVal, zeroVec)
								            curMaskOut = curMaskOut or (curMaskInVec and bitMask)

								            curNorm = curNorm + curKernelVal;
								        end
							        end
						    	end
						    end
							--end
						end
					    emit quote
				    		curImOut = curImOut/curNorm
				    		curVarOut = curVarOut/(curNorm*curNorm)
						end
						for i = 0, vectorWidth-1 do
				    		emit quote
					    		outputIm[y*width + x+i] = curImOut[i]	
							    outputVar[y*width + x+i] = curVarOut[i]
				    			outputMask[y*width + x+i] = curMaskOut[i]
				    		end
				    	end
				    end				
				end
			end
		end

		var t2 = C.clock()

		C.printf("outputIm[boundingBox*width + boundingBox] = %f, computation took: %f ms; %f, %d\n",
			outputIm[boundingBox*width + boundingBox], (t2-t1)/1000.0, outputVar[boundingBox*width + boundingBox], outputMask[boundingBox*width + boundingBox])


--		C.printf("\n\nVectorized %d wide masked image lin combo %dx%d blur Terra NO UNROLL:\n", vectorWidth, 1+2*boundingBox, 1+2*boundingBox)
--		C.printf("outputIm[boundingBox*width + boundingBox] = %f, computation took: %f ms\n", outputIm[boundingBox*width + boundingBox],  (t2-t1)/1000.0)
--		C.printf("C.CLOCKS_PER_SEC = %d\n", C.CLOCKS_PER_SEC)
--
--		C.printf("Image plane, 10x10 box begining at (boundingBox,boundingBox)\n")
--		for i=boundingBox,boundingBox+10 do
--			for j=boundingBox,boundingBox+10 do
--				C.printf("%f\t", outputIm[i*width + j])
--			end
--			C.printf("\n")
--		end
--
--		C.printf("Variance plane, 10x10 box begining at (boundingBox,boundingBox)\n")
--		for i=boundingBox,boundingBox+10 do
--			for j=boundingBox,boundingBox+10 do
--				C.printf("%f\t", outputVar[i*width + j])
--			end
--			C.printf("\n")
--		end
--
--		C.printf("Mask plane, 10x10 box begining at (boundingBox,boundingBox)\n")
--		for i=boundingBox,boundingBox+10 do
--			for j=boundingBox,boundingBox+10 do
--				C.printf("%d\t", outputMask[i*width + j])
--			end
--			C.printf("\n")
--		end

	    C.free(inputIm)
		C.free(outputIm)
	    C.free(inputVar)
		C.free(outputVar)
	    C.free(inputMask)
		C.free(outputMask)
		C.free([kernelVals])

	end
	return blurMaskedImageLinComboVectorize
end


local function metaProgramBlurMaskedImageInterpolate(boundingBox, vectorWidth, method, interpDist)
	local getBitMask = genGetBitMask(vectorWidth) --getBitMask is a terra function
	local terra blurMaskedImageInterpolate()
	--	var boundingBox : int = 2
		var width : int = 2048
		var height : int = 1489

		var inputIm : &float
		inputIm = [&float](C.malloc(sizeof(float)*width*height))
		var inputVar : &float
		inputVar = [&float](C.malloc(sizeof(float)*width*height))
		var inputMask : &uint16
		inputMask = [&uint16](C.malloc(sizeof(uint16)*width*height))

		var outputIm : &float
		outputIm = [&float](C.malloc(sizeof(float)*width*height))
		var outputVar : &float
		outputVar = [&float](C.malloc(sizeof(float)*width*height))
		var outputMask : &uint16
		outputMask = [&uint16](C.malloc(sizeof(uint16)*width*height))

		for x = 0, width do
			for y = 0, height do
				inputIm[y*width + x] = [float](y*x*C.cos(x/(y+1)))
				inputVar[y*width + x] = [float](x*y + 2)
				inputMask[y*width + x] = [uint16](x + y)


				outputIm[y*width + x] = 0.0f
				outputVar[y*width + x] = 0.0f
				outputMask[y*width + x] = 0
			end
		end


		var t1 = C.clock()


		var zeroVec : vector(float, vectorWidth) = [zeroFloatVec(vectorWidth)]
		for y = boundingBox, height-boundingBox, interpDist do
			for x = boundingBox, (width-boundingBox), interpDist do
 				escape
 					if(method == "original") then
 						local curKernelVal = symbol()
 						local curKernelValTemp = symbol()
 						local curImOut = symbol()
 						local curVarOut= symbol()
 						local curMaskOut = symbol()
 						local curNorm = symbol()
	 					emit quote
			 				var [curKernelVal] : vector(float,vectorWidth)
							var [curKernelValTemp] : float[vectorWidth]

							var [curImOut] : vector(float,vectorWidth) = [zeroFloatVec(vectorWidth)]
			--				var curImOut : &vector(float,vectorWidth) = [&vector(float,vectorWidth)](&outputIm[y*width + x])
							var [curVarOut] : vector(float,vectorWidth) = [zeroFloatVec(vectorWidth)]
							var [curMaskOut] : vector(uint16, vectorWidth) = [zeroUInt16Vec(vectorWidth)]

							var [curNorm] : vector(float,vectorWidth) = [zeroFloatVec(vectorWidth)]
			--				var curImOutVec = @[&vector(float,4)](&outputIm[y*width + x])
						end
						for j = -boundingBox, boundingBox do
					        for i = -boundingBox, boundingBox do
					        	for k = 0, vectorWidth - 1 do
					        		emit quote
						        		curKernelValTemp[k] = polynomial1(x+k, y)*[kernel1(i, j)] +
							                polynomial2(x+k, y)*[kernel2(i, j)] + polynomial3(x+k, y)*[kernel3(i, j)] + 
							                polynomial4(x+k, y)*[kernel4(i, j)] + polynomial5(x+k, y)*[kernel5(i, j)];
					        		end
					        	end
					        	emit quote
						            curKernelVal = @[&vector(float,vectorWidth)](&curKernelValTemp[0])

						            var curImInVec = @[&vector(float,vectorWidth)](&inputIm[(y+j)*width + x+i])
						            var curVarInVec = @[&vector(float,vectorWidth)](&inputVar[(y+j)*width + x+i])
						            var curMaskInVec = @[&vector(uint16,vectorWidth)](&inputMask[(y+j)*width + x+i])


						            curImOut = curImOut + curImInVec*curKernelVal; 
						            curVarOut = curVarOut + curVarInVec*curKernelVal*curKernelVal;

						            var bitMask : vector(uint16, vectorWidth) = getBitMask(curKernelVal, zeroVec)
						            curMaskOut = curMaskOut or (curMaskInVec and bitMask)

						            curNorm = curNorm + curKernelVal;
						        end
					        end
					    end
					    emit quote
				    		curImOut = curImOut/curNorm
				    		curVarOut = curVarOut/(curNorm*curNorm)
						end
						for i = 0, vectorWidth-1 do
				    		emit quote
					    		outputIm[y*width + x+i] = curImOut[i]	
							    outputVar[y*width + x+i] = curVarOut[i]
				    			outputMask[y*width + x+i] = curMaskOut[i]
				    		end
				    	end
				    end				
				end
			end
		end

		var t2 = C.clock()
		C.printf("\n\nVectorized %d wide masked image lin combo %dx%d blur Terra:\n", vectorWidth, 1+2*boundingBox, 1+2*boundingBox)
		C.printf("outputIm[boundingBox*width + boundingBox] = %f, computation took: %f ms\n", outputIm[boundingBox*width + boundingBox],  (t2-t1)/1000.0)
		C.printf("C.CLOCKS_PER_SEC = %d\n", C.CLOCKS_PER_SEC)

		C.printf("Image plane, 10x10 box begining at (boundingBox,boundingBox)\n")
		for i=boundingBox,boundingBox+10 do
			for j=boundingBox,boundingBox+10 do
				C.printf("%f\t", outputIm[i*width + j])
			end
			C.printf("\n")
		end

		C.printf("Variance plane, 10x10 box begining at (boundingBox,boundingBox)\n")
		for i=boundingBox,boundingBox+10 do
			for j=boundingBox,boundingBox+10 do
				C.printf("%f\t", outputVar[i*width + j])
			end
			C.printf("\n")
		end

		C.printf("Mask plane, 10x10 box begining at (boundingBox,boundingBox)\n")
		for i=boundingBox,boundingBox+10 do
			for j=boundingBox,boundingBox+10 do
				C.printf("%d\t", outputMask[i*width + j])
			end
			C.printf("\n")
		end

	    C.free(inputIm)
		C.free(outputIm)
	    C.free(inputVar)
		C.free(outputVar)
	    C.free(inputMask)
		C.free(outputMask)
	end
	return blurMaskedImageInterpolate
end



--local blurMaskedImageLinComboVectorize = metaProgramBlurMaskedImageLinComboVectorize(1, 8, "original")
--blurMaskedImageLinComboVectorize()

--local blurMaskedImageLinComboVectorize = metaProgramBlurMaskedImageLinComboVectorize(5, 8, "original")
--blurMaskedImageLinComboVectorize()

--local blurMaskedImageLinComboVectorize = metaProgramBlurMaskedImageLinComboVectorize(3, 8, "original")
--blurMaskedImageLinComboVectorize()
--
--local blurMaskedImageLinComboVectorize = metaProgramBlurMaskedImageLinComboVectorize(4, 8, "original")
--blurMaskedImageLinComboVectorize()
--
--local blurMaskedImageLinComboVectorize = metaProgramBlurMaskedImageLinComboVectorize(5, 8, "original")
--blurMaskedImageLinComboVectorize()
--
--local blurMaskedImageLinComboVectorize = metaProgramBlurMaskedImageLinComboVectorize(6, 8, "original")
--blurMaskedImageLinComboVectorize()
--
--local blurMaskedImageLinComboVectorize = metaProgramBlurMaskedImageLinComboVectorize(7, 8, "original")
--blurMaskedImageLinComboVectorize()
--
--local blurMaskedImageLinComboVectorize = metaProgramBlurMaskedImageLinComboVectorize(8, 8, "original")
--blurMaskedImageLinComboVectorize()
--
--local blurMaskedImageLinComboVectorize = metaProgramBlurMaskedImageLinComboVectorize(9, 8, "original")
--blurMaskedImageLinComboVectorize()
--
--local blurMaskedImageLinComboVectorize = metaProgramBlurMaskedImageLinComboVectorize(13, 8, "original")
--blurMaskedImageLinComboVectorize()



--local blurMaskedImageLinComboVectorize = metaProgramBlurMaskedImageLinComboVectorize(5, 8, "original")
--blurMaskedImageLinComboVectorize()
--
--
--local blurMaskedImageLinComboVectorizeNoUnroll = metaProgramBlurMaskedImageLinComboVectorizeNoUnroll(5, 8, "original")
--blurMaskedImageLinComboVectorizeNoUnroll()

local blurMaskedImageLinComboVectorize = metaProgramBlurMaskedImageLinComboVectorize(9, 8, "original", 5)
blurMaskedImageLinComboVectorize()
local blurMaskedImageLinComboVectorize = metaProgramBlurMaskedImageLinComboVectorize(9, 8, "original", 4)
blurMaskedImageLinComboVectorize()
local blurMaskedImageLinComboVectorize = metaProgramBlurMaskedImageLinComboVectorize(9, 8, "original", 3)
blurMaskedImageLinComboVectorize()
local blurMaskedImageLinComboVectorize = metaProgramBlurMaskedImageLinComboVectorize(9, 8, "original", 2)
blurMaskedImageLinComboVectorize()
local blurMaskedImageLinComboVectorize = metaProgramBlurMaskedImageLinComboVectorize(9, 8, "original", 1)
blurMaskedImageLinComboVectorize()

local NUMBER_OF_RUNS = 1




--blurMaskedImageLinComboVectorizeNoUnroll()
print("NoUnroll")
local totalTimeNoUnroll = 0
for i = 1, NUMBER_OF_RUNS do
	local blurMaskedImageLinComboVectorizeNoUnroll = metaProgramBlurMaskedImageLinComboVectorizeNoUnroll(9, 8, "original")
	local t3 = terralib.currenttimeinseconds()
	blurMaskedImageLinComboVectorizeNoUnroll:compile()
	blurMaskedImageLinComboVectorizeNoUnroll()
	local t4 = terralib.currenttimeinseconds()
	totalTimeNoUnroll = totalTimeNoUnroll + t4 - t3
end




print("Unroll, polynomial is calculated at (x, y) then array accessed")
local totalTimeUnrollFlexiblePoly = 0
for i = 1, NUMBER_OF_RUNS do
	local blurMaskedImageLinComboVectorizeFlexPoly = metaProgramBlurMaskedImageLinComboVectorizeUnrollFlexiblePolyCoef(9, 8, "original")
	local t1 = terralib.currenttimeinseconds()
	blurMaskedImageLinComboVectorizeFlexPoly:compile()
	blurMaskedImageLinComboVectorizeFlexPoly()
	local t2 = terralib.currenttimeinseconds()
	totalTimeUnrollFlexiblePoly = totalTimeUnrollFlexiblePoly + t2-t1
end

print("Unroll, polynomial AND kernel are calculated at (x, y) then array accessed")
local totalTimeUnrollFlexiblePolyAndKernel = 0
for i = 1, NUMBER_OF_RUNS do
	local blurMaskedImageLinComboVectorizeFlexPolyAndKernel = metaProgramBlurMaskedImageLinComboVectorizeUnrollFlexiblePolyAndKernel(9, 8, "original")
	local t1 = terralib.currenttimeinseconds()
	blurMaskedImageLinComboVectorizeFlexPolyAndKernel:compile()
	blurMaskedImageLinComboVectorizeFlexPolyAndKernel()
	local t2 = terralib.currenttimeinseconds()
	totalTimeUnrollFlexiblePolyAndKernel = totalTimeUnrollFlexiblePolyAndKernel + t2-t1
end

--blurMaskedImageLinComboVectorize:compile()
print("Unroll")
local totalTimeUnroll = 0
for i = 1, NUMBER_OF_RUNS do
	local blurMaskedImageLinComboVectorize = metaProgramBlurMaskedImageLinComboVectorize(9, 8, "original")
	local t1 = terralib.currenttimeinseconds()
	blurMaskedImageLinComboVectorize:compile()
	blurMaskedImageLinComboVectorize()
	local t2 = terralib.currenttimeinseconds()
	totalTimeUnroll = totalTimeUnroll + t2-t1
end


print("unrolling time = "..(totalTimeUnroll/NUMBER_OF_RUNS*1000).."ms")
print("unrolling, polynomial array access; time = "..(totalTimeUnrollFlexiblePoly/NUMBER_OF_RUNS*1000).."ms")
print("unrolling, polynomial AND kernel array access; time = "..(totalTimeUnrollFlexiblePolyAndKernel/NUMBER_OF_RUNS*1000).."ms")
print("No unrolling time = "..(totalTimeNoUnroll/NUMBER_OF_RUNS*1000).."ms")
--print("no unrolling time = "..((t4-t3)/NUMBER_OF_RUNS*1000).."ms")


--				var polynomial1plusx0 : float = polynomial1(x, y);
--				var polynomial2plusx0 : float = polynomial2(x, y);
--				var polynomial3plusx0 : float = polynomial3(x, y);
--				var polynomial4plusx0 : float = polynomial4(x, y);
--				var polynomial5plusx0 : float = polynomial5(x, y);
--				 
--				var polynomial1plusx1 : float = polynomial1(x+1, y);
--				var polynomial2plusx1 : float = polynomial2(x+1, y);
--				var polynomial3plusx1 : float = polynomial3(x+1, y);
--				var polynomial4plusx1 : float = polynomial4(x+1, y);
--				var polynomial5plusx1 : float = polynomial5(x+1, y);
--
--				var polynomial1plusx2 : float = polynomial1(x+2, y);
--				var polynomial2plusx2 : float = polynomial2(x+2, y);
--				var polynomial3plusx2 : float = polynomial3(x+2, y);
--				var polynomial4plusx2 : float = polynomial4(x+2, y);
--				var polynomial5plusx2 : float = polynomial5(x+2, y);
--
--				var polynomial1plusx3 : float = polynomial1(x+3, y);
--				var polynomial2plusx3 : float = polynomial2(x+3, y);
--				var polynomial3plusx3 : float = polynomial3(x+3, y);
--				var polynomial4plusx3 : float = polynomial4(x+3, y);
--				var polynomial5plusx3 : float = polynomial5(x+3, y);
--
--				var polynomial1plusx4 : float = polynomial1(x+4, y);
--				var polynomial2plusx4 : float = polynomial2(x+4, y);
--				var polynomial3plusx4 : float = polynomial3(x+4, y);
--				var polynomial4plusx4 : float = polynomial4(x+4, y);
--				var polynomial5plusx4 : float = polynomial5(x+4, y);
--
--				var polynomial1plusx5 : float = polynomial1(x+5, y);
--				var polynomial2plusx5 : float = polynomial2(x+5, y);
--				var polynomial3plusx5 : float = polynomial3(x+5, y);
--				var polynomial4plusx5 : float = polynomial4(x+5, y);
--				var polynomial5plusx5 : float = polynomial5(x+5, y);
--
--				var polynomial1plusx6 : float = polynomial1(x+6, y);
--				var polynomial2plusx6 : float = polynomial2(x+6, y);
--				var polynomial3plusx6 : float = polynomial3(x+6, y);
--				var polynomial4plusx6 : float = polynomial4(x+6, y);
--				var polynomial5plusx6 : float = polynomial5(x+6, y);
--
--				var polynomial1plusx7 : float = polynomial1(x+7, y);
--				var polynomial2plusx7 : float = polynomial2(x+7, y);
--				var polynomial3plusx7 : float = polynomial3(x+7, y);
--				var polynomial4plusx7 : float = polynomial4(x+7, y);
--				var polynomial5plusx7 : float = polynomial5(x+7, y);


