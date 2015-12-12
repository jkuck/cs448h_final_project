local C = terralib.includecstring [[
    #include<stdio.h>
    #include<stdlib.h>
    #include<time.h>
    #include<string.h>
    #include<math.h>

]]

--local boundingBox = 2
local vectorWidth = 8
--local luaNumberOfBasisKernels = 27
--local luaBoundingBox = 11
local luaNumberOfBasisKernels = 5
local luaBoundingBox = 2
local method = "original"

local terra getInfoFromLSST(inputImg:&float, inputVar:&float, inputMask:&uint16,
							outputImg:&float, outputVar:&float, outputMask:&uint16,
                    		imageWidth:int, imageHeight:int, kernelArray:&&float , 
                    		funcParams:&float, numberOfBasisKernels:int,
							kernelWidth:int, kernelHeight:int, numberOfFuncCoef:int)
	C.printf("received function paramters in Terra:\n")
	for i = 0, numberOfBasisKernels do
		C.printf("function %d:", i)
		for j = 0, numberOfFuncCoef do
			C.printf("%f\t", funcParams[i*numberOfFuncCoef + j])
		end
		C.printf("\n")
	end
end


local function zerofloatVec(vectorWidth)
	local zeros = {}
	for i=1,vectorWidth do zeros[i] = `0.f end
	return `vector(zeros)
end

local function zeroUInt16Vec(vectorWidth)
	local zeros = {}
	for i=1,vectorWidth do zeros[i] = `0 end
	return `vector(zeros)
end

local function gen_convolveOneVec_original(x, y, inputImg, inputVar, inputMask,
							outputImg, outputVar, outputMask,
                    		imageWidth, imageHeight, kernelArray, 
                    		funcParams, numberOfBasisKernels,
							kernelWidth, kernelHeight, numberOfFuncCoef, zeroVec, getBitMask)
	local returnQuote = {}
    local curQuote

	local curKernelVal = symbol()
	local curKernelValTemp = symbol()
	local curImOut = symbol()
	local curVarOut= symbol()
	local curMaskOut = symbol()
	local curNorm = symbol()
		
	curQuote = quote
		var [curKernelVal] : vector(float,vectorWidth)
		var [curKernelValTemp] : float[vectorWidth]

		var [curImOut] : vector(float,vectorWidth) = [zerofloatVec(vectorWidth)]
--				var curImOut : &vector(float,vectorWidth) = [&vector(float,vectorWidth)](&outputImg[y*imageWidth + x])
		var [curVarOut] : vector(float,vectorWidth) = [zerofloatVec(vectorWidth)]
		var [curMaskOut] : vector(uint16, vectorWidth) = [zeroUInt16Vec(vectorWidth)]

		var [curNorm] : vector(float,vectorWidth) = [zerofloatVec(vectorWidth)]
--				var curImOutVec = @[&vector(float,4)](&outputImg[y*imageWidth + x])
	end
	table.insert(returnQuote, curQuote)


	for j = -luaBoundingBox, luaBoundingBox do
		for i = -luaBoundingBox, luaBoundingBox do
        	for k = 0, vectorWidth - 1 do
        		curQuote = quote
						curKernelValTemp[k] = kernelArray[0][(j+luaBoundingBox)*kernelWidth+(i+luaBoundingBox)]*(funcParams[0*numberOfFuncCoef + 0]
        				+ funcParams[0*numberOfFuncCoef + 1]*x + funcParams[0*numberOfFuncCoef + 2]*y
        				+ funcParams[0*numberOfFuncCoef + 3]*x*x + funcParams[0*numberOfFuncCoef + 4]*x*y
        				+ funcParams[0*numberOfFuncCoef + 5]*y*y + funcParams[0*numberOfFuncCoef + 6]*x*x*x
        				+ funcParams[0*numberOfFuncCoef + 7]*x*x*y + funcParams[0*numberOfFuncCoef + 8]*x*y*y
        				+ funcParams[0*numberOfFuncCoef + 9]*y*y*y)	
        		end
				table.insert(returnQuote, curQuote)
        		for l = 1, luaNumberOfBasisKernels-1 do 
        			curQuote = quote
        			curKernelValTemp[k] = curKernelValTemp[k] + kernelArray[l][(j+luaBoundingBox)*kernelWidth+(i+luaBoundingBox)]*(funcParams[l*numberOfFuncCoef + 0]
        				+ funcParams[l*numberOfFuncCoef + 1]*x + funcParams[l*numberOfFuncCoef + 2]*y
        				+ funcParams[l*numberOfFuncCoef + 3]*x*x + funcParams[l*numberOfFuncCoef + 4]*x*y
        				+ funcParams[l*numberOfFuncCoef + 5]*y*y + funcParams[l*numberOfFuncCoef + 6]*x*x*x
        				+ funcParams[l*numberOfFuncCoef + 7]*x*x*y + funcParams[l*numberOfFuncCoef + 8]*x*y*y
        				+ funcParams[l*numberOfFuncCoef + 9]*y*y*y)
        			end
    				table.insert(returnQuote, curQuote)
        		end
        	end
        	curQuote = quote
	            curKernelVal = @[&vector(float,vectorWidth)](&curKernelValTemp[0])

	            var curImInVec = @[&vector(float,vectorWidth)](&inputImg[(y+j)*imageWidth + x+i])
	            var curVarInVec = @[&vector(float,vectorWidth)](&inputVar[(y+j)*imageWidth + x+i])
	            var curMaskInVec = @[&vector(uint16,vectorWidth)](&inputMask[(y+j)*imageWidth + x+i])


	            curImOut = curImOut + curImInVec*curKernelVal; 
	            curVarOut = curVarOut + curVarInVec*curKernelVal*curKernelVal;

	            var bitMask : vector(uint16, vectorWidth) = getBitMask(curKernelVal, zeroVec)
	            curMaskOut = curMaskOut or (curMaskInVec and bitMask)

	            curNorm = curNorm + curKernelVal;
	        end
	        table.insert(returnQuote, curQuote)

    	end
    end
    curQuote = quote
		curImOut = curImOut/curNorm
		curVarOut = curVarOut/(curNorm*curNorm)
	end
	table.insert(returnQuote, curQuote)
	for i = 0, vectorWidth-1 do
		curQuote = quote
    		outputImg[y*imageWidth + x+i] = curImOut[i]	
		    outputVar[y*imageWidth + x+i] = curVarOut[i]
			outputMask[y*imageWidth + x+i] = curMaskOut[i]
		end
		table.insert(returnQuote, curQuote)
	end

end



local function genGetBitMask(vectorWidth)
	local terra returnFunc(a : vector(float,vectorWidth), b : vector(float,vectorWidth))
	    return terralib.select(a ~= b, [vector(uint16,vectorWidth)](0xFFFFULL),[vector(uint16,vectorWidth)](0) ) 
	end
	return returnFunc
end

    --store the jth coefficient of the ith function in
    --funcParams[i*numberOfFuncCoef + j]
local terra getInfoFromLSST(inputImg:&float, inputVar:&float, inputMask:&uint16,
							outputImg:&float, outputVar:&float, outputMask:&uint16,
                    		imageWidth:int, imageHeight:int, kernelArray:&&float, 
                    		funcParams:&float, numberOfBasisKernels:int,
							kernelWidth:int, kernelHeight:int, numberOfFuncCoef:int)
	C.printf("Begin Terra function\n")
	C.printf("imageWidth = %d, imageHeight = %d, numberOfBasisKernels = %d, kernelWidth = %d, kernelHeight = %d, numberOfFuncCoef = %d\n", imageWidth,
		imageHeight, numberOfBasisKernels, kernelWidth, kernelHeight, numberOfFuncCoef)
	escape
		local getBitMask = genGetBitMask(vectorWidth) --getBitMask is a terra function
--		assert(kernelWidth == kernelHeight, "kernel is not square")
--		assert(kernelWidth%2 == 1, "kernel width is even")
		local boundingBox = symbol()
		emit quote
			var [boundingBox] = (kernelWidth - 1)/2

			var t1 = C.clock()


			var zeroVec : vector(float, vectorWidth) = [zerofloatVec(vectorWidth)]
			for y = boundingBox, imageHeight-boundingBox do
				for x = boundingBox, (imageWidth-boundingBox), vectorWidth do
	 				[gen_convolveOneVec_original(x, y, inputImg, inputVar, inputMask,
							outputImg, outputVar, outputMask,
                    		imageWidth, imageHeight, kernelArray, 
                    		funcParams, numberOfBasisKernels,
							kernelWidth, kernelHeight, numberOfFuncCoef, zeroVec, getBitMask)]	 						
				end
			end

			var t2 = C.clock()
			C.printf("\n\nVectorized %d wide masked image lin combo %dx%d blur Terra:\n", vectorWidth, 1+2*boundingBox, 1+2*boundingBox)
			C.printf("outputImg[boundingBox*imageWidth + boundingBox] = %f, computation took: %f ms\n", outputImg[boundingBox*imageWidth + boundingBox],  (t2-t1)/1000.0)
			C.printf("C.CLOCKS_PER_SEC = %d\n", C.CLOCKS_PER_SEC)
	
			C.printf("Input image plane, 10x10 box begining at (boundingBox,boundingBox)\n")
			for i=boundingBox,boundingBox+10 do
				for j=boundingBox,boundingBox+10 do
					C.printf("%f\t", inputImg[i*imageWidth + j])
				end
				C.printf("\n")
			end

			C.printf("Output image plane, 10x10 box begining at (boundingBox,boundingBox)\n")
			for i=boundingBox,boundingBox+10 do
				for j=boundingBox,boundingBox+10 do
					C.printf("%f\t", outputImg[i*imageWidth + j])
				end
				C.printf("\n")
			end
	
			C.printf("Variance plane, 10x10 box begining at (boundingBox,boundingBox)\n")
			for i=boundingBox,boundingBox+10 do
				for j=boundingBox,boundingBox+10 do
					C.printf("%f\t", outputVar[i*imageWidth + j])
				end
				C.printf("\n")
			end
	
			C.printf("Mask plane, 10x10 box begining at (boundingBox,boundingBox)\n")
			for i=boundingBox,boundingBox+10 do
				for j=boundingBox,boundingBox+10 do
					C.printf("%d\t", outputMask[i*imageWidth + j])
				end
				C.printf("\n")
			end

			C.printf("Function coefficients:\n")
			for i = 0, numberOfBasisKernels do
				for j = 0, numberOfFuncCoef do
					C.printf("%f\t", funcParams[i*numberOfFuncCoef + j])
				end
				C.printf("\n")
			end

			C.printf("Basis kernels:\n")
			for h = 0, numberOfBasisKernels do
				for j = 0, (luaBoundingBox*2+1) do
					for i = 0, (luaBoundingBox*2+1) do
						C.printf("%f\t", kernelArray[h][j*(luaBoundingBox*2+1) + i])
					end
					C.printf("\n")
				end
				C.printf("\n")
			end



		end
	end
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

--getInfoFromLSST:printpretty()

--the first term is the name in the .o and .h files, the second name is the name in this file
terralib.saveobj("blurExampleTerraStandalone.o",{ terraFuncNameInC = getInfoFromLSST })




