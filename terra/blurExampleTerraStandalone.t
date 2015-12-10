local C = terralib.includecstring [[
    #include<stdio.h>
    #include<stdlib.h>
    #include<time.h>
    #include<string.h>
    #include<math.h>

]]

--local boundingBox = 2
local vectorWidth = 8
local kernelSize = 19
local luaNumberOfBasisKernels = 5
local luaNumberOfFuncParams = 10
local luaBoundingBox = (kernelSize-1)/2
local boundingBox = (kernelSize-1)/2
local kernelWidth = kernelSize
local luaKernelArea = kernelSize*kernelSize
--local luaNumberOfBasisKernels = 5
--local luaNumberOfFuncParams = 10
--local luaBoundingBox = 2
local interpDist = 10
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


local function genGetBitMask(vectorWidth)
	local terra returnFunc(a : vector(float,vectorWidth), b : vector(float,vectorWidth))
	    return terralib.select(a ~= b, [vector(uint16,vectorWidth)](0xFFFFULL),[vector(uint16,vectorWidth)](0) ) 
	end
	return returnFunc
end


local function genGetBitMaskDouble(vectorWidth)
	local terra returnFunc(a : vector(double,vectorWidth), b : vector(double,vectorWidth))
	    return terralib.select(a ~= b, [vector(uint16,vectorWidth)](0xFFFFULL),[vector(uint16,vectorWidth)](0) ) 
	end
	return returnFunc
end

    --store the jth coefficient of the ith function in
    --funcParams[i*numberOfFuncCoef + j]
local function genConvolve(kenelSize)
	local luaBoundingBox = (kenelSize-1)/2
	local returnFunc = terra(inputImg:&float, inputVar:&float, inputMask:&uint16,
								outputImg:&float, outputVar:&float, outputMask:&uint16,
	                    		imageWidth:int, imageHeight:int, kernelArray:&&float, 
	                    		funcParams:&float, numberOfBasisKernels:int,
								kernelWidth:int, kernelHeight:int, numberOfFuncCoef:int)
--		C.printf("Begin Terra function\n")
--		C.printf("imageWidth = %d, imageHeight = %d, numberOfBasisKernels = %d, kernelWidth = %d, kernelHeight = %d, numberOfFuncCoef = %d\n", imageWidth,
--			imageHeight, numberOfBasisKernels, kernelWidth, kernelHeight, numberOfFuncCoef)
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
		 				escape
	--	 					if(method == "original") then
		 					if(true) then
		 						local curKernelVal = symbol()
		 						local curKernelValTemp = symbol()
		 						local curImOut = symbol()
		 						local curVarOut= symbol()
		 						local curMaskOut = symbol()
		 						local curNorm = symbol()
			 					emit quote
					 				var [curKernelVal] : vector(float,vectorWidth)
									var [curKernelValTemp] : float[vectorWidth]

									var [curImOut] : vector(float,vectorWidth) = [zerofloatVec(vectorWidth)]
					--				var curImOut : &vector(float,vectorWidth) = [&vector(float,vectorWidth)](&outputImg[y*imageWidth + x])
									var [curVarOut] : vector(float,vectorWidth) = [zerofloatVec(vectorWidth)]
									var [curMaskOut] : vector(uint16, vectorWidth) = [zeroUInt16Vec(vectorWidth)]

									var [curNorm] : vector(float,vectorWidth) = [zerofloatVec(vectorWidth)]
					--				var curImOutVec = @[&vector(float,4)](&outputImg[y*imageWidth + x])
								end
	--							for j = -boundingBox, boundingBox do
	--						        for i = -boundingBox, boundingBox do
								emit quote
								for j = -luaBoundingBox, luaBoundingBox +1 do
								escape
									emit quote
							        for i = -luaBoundingBox, luaBoundingBox + 1 do
							        escape
							        	for k = 0, vectorWidth - 1 do
							        		emit quote
	       										curKernelValTemp[k] = kernelArray[0][(j+boundingBox)*kernelWidth+(i+boundingBox)]*(funcParams[0*numberOfFuncCoef + 0]
							        				+ funcParams[0*numberOfFuncCoef + 1]*x + funcParams[0*numberOfFuncCoef + 2]*y
							        				+ funcParams[0*numberOfFuncCoef + 3]*x*x + funcParams[0*numberOfFuncCoef + 4]*x*y
							        				+ funcParams[0*numberOfFuncCoef + 5]*y*y + funcParams[0*numberOfFuncCoef + 6]*x*x*x
							        				+ funcParams[0*numberOfFuncCoef + 7]*x*x*y + funcParams[0*numberOfFuncCoef + 8]*x*y*y
							        				+ funcParams[0*numberOfFuncCoef + 9]*y*y*y)	
							        		end
	--						        		for l = 1, numberOfBasisKernels-1 do 
							        		for l = 1, luaNumberOfBasisKernels-1 do 
							        			emit quote
							        			curKernelValTemp[k] = curKernelValTemp[k] + kernelArray[l][(j+boundingBox)*kernelWidth+(i+boundingBox)]*(funcParams[l*numberOfFuncCoef + 0]
							        				+ funcParams[l*numberOfFuncCoef + 1]*x + funcParams[l*numberOfFuncCoef + 2]*y
							        				+ funcParams[l*numberOfFuncCoef + 3]*x*x + funcParams[l*numberOfFuncCoef + 4]*x*y
							        				+ funcParams[l*numberOfFuncCoef + 5]*y*y + funcParams[l*numberOfFuncCoef + 6]*x*x*x
							        				+ funcParams[l*numberOfFuncCoef + 7]*x*x*y + funcParams[l*numberOfFuncCoef + 8]*x*y*y
							        				+ funcParams[l*numberOfFuncCoef + 9]*y*y*y)
							        			end
							        		end
							        	end
							        	emit quote
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
							        end
							    	end
							        end
							    end
								end
								end
							    emit quote
						    		curImOut = curImOut/curNorm
						    		curVarOut = curVarOut/(curNorm*curNorm)
								end
								for i = 0, vectorWidth-1 do
						    		emit quote
							    		outputImg[y*imageWidth + x+i] = curImOut[i]	
									    outputVar[y*imageWidth + x+i] = curVarOut[i]
						    			outputMask[y*imageWidth + x+i] = curMaskOut[i]
						    		end
						    	end
						    end				
						end
					end
				end

				var t2 = C.clock()
--				C.printf("\n\nVectorized %d wide masked image lin combo %dx%d blur Terra:\n", vectorWidth, 1+2*boundingBox, 1+2*boundingBox)
--				C.printf("outputImg[boundingBox*imageWidth + boundingBox] = %f, computation took: %f ms\n", outputImg[boundingBox*imageWidth + boundingBox],  (t2-t1)/1000.0)
--				C.printf("C.CLOCKS_PER_SEC = %d\n", C.CLOCKS_PER_SEC)
--		
--				C.printf("Input image plane, 10x10 box begining at (boundingBox,boundingBox)\n")
--				for i=boundingBox,boundingBox+10 do
--					for j=boundingBox,boundingBox+10 do
--						C.printf("%f\t", inputImg[i*imageWidth + j])
--					end
--					C.printf("\n")
--				end
--
--				C.printf("Output image plane, 10x10 box begining at (boundingBox,boundingBox)\n")
--				for i=boundingBox,boundingBox+10 do
--					for j=boundingBox,boundingBox+10 do
--						C.printf("%f\t", outputImg[i*imageWidth + j])
--					end
--					C.printf("\n")
--				end
--		
--				C.printf("Output Variance plane, 10x10 box begining at (boundingBox,boundingBox)\n")
--				for i=boundingBox,boundingBox+10 do
--					for j=boundingBox,boundingBox+10 do
--						C.printf("%f\t", outputVar[i*imageWidth + j])
--					end
--					C.printf("\n")
--				end
--		
--				C.printf("Output Mask plane, 10x10 box begining at (boundingBox,boundingBox)\n")
--				for i=boundingBox,boundingBox+10 do
--					for j=boundingBox,boundingBox+10 do
--						C.printf("%d\t", outputMask[i*imageWidth + j])
--					end
--					C.printf("\n")
--				end
--
--				C.printf("Function coefficients:\n")
--				for i = 0, numberOfBasisKernels do
--					for j = 0, numberOfFuncCoef do
--						C.printf("%f\t", funcParams[i*numberOfFuncCoef + j])
--					end
--					C.printf("\n")
--				end
--
--				C.printf("Basis kernels:\n")
--				for h = 0, numberOfBasisKernels do
--					for j = 0, (luaBoundingBox*2+1) do
--						for i = 0, (luaBoundingBox*2+1) do
--							C.printf("%f\t", kernelArray[h][j*(luaBoundingBox*2+1) + i])
--						end
--						C.printf("\n")
--					end
--					C.printf("\n")
--				end



			end
		end
	end

	return returnFunc
end


local terra convolveDoublePrecisionKernel(inputImg:&float, inputVar:&float, inputMask:&uint16,
							outputImg:&float, outputVar:&float, outputMask:&uint16,
                    		imageWidth:int, imageHeight:int, kernelArray:&&double, 
                    		funcParams:&float, numberOfBasisKernels:int,
							kernelWidth:int, kernelHeight:int, numberOfFuncCoef:int)
	C.printf("Begin Terra function\n")
	C.printf("imageWidth = %d, imageHeight = %d, numberOfBasisKernels = %d, kernelWidth = %d, kernelHeight = %d, numberOfFuncCoef = %d\n", imageWidth,
		imageHeight, numberOfBasisKernels, kernelWidth, kernelHeight, numberOfFuncCoef)
	escape
		local getBitMask = genGetBitMaskDouble(vectorWidth) --getBitMask is a terra function
--		assert(kernelWidth == kernelHeight, "kernel is not square")
--		assert(kernelWidth%2 == 1, "kernel width is even")
		local boundingBox = symbol()
		emit quote
			var [boundingBox] = (kernelWidth - 1)/2

			var t1 = C.clock()


			var zeroVec : vector(double, vectorWidth) = [zeroDoubleVec(vectorWidth)]
			for y = boundingBox, imageHeight-boundingBox do
				for x = boundingBox, (imageWidth-boundingBox), vectorWidth do
	 				escape
--	 					if(method == "original") then
	 					if(true) then
	 						local curKernelVal = symbol()
	 						local curKernelValTemp = symbol()
	 						local curImOut = symbol()
	 						local curVarOut= symbol()
	 						local curMaskOut = symbol()
	 						local curNorm = symbol()
		 					emit quote
				 				var [curKernelVal] : vector(double,vectorWidth)
								var [curKernelValTemp] : double[vectorWidth]

								var [curImOut] : vector(float,vectorWidth) = [zerofloatVec(vectorWidth)]
				--				var curImOut : &vector(float,vectorWidth) = [&vector(float,vectorWidth)](&outputImg[y*imageWidth + x])
								var [curVarOut] : vector(float,vectorWidth) = [zerofloatVec(vectorWidth)]
								var [curMaskOut] : vector(uint16, vectorWidth) = [zeroUInt16Vec(vectorWidth)]

								var [curNorm] : vector(double,vectorWidth) = [zeroDoubleVec(vectorWidth)]
							end


--							for j = -boundingBox, boundingBox do
--						        for i = -boundingBox, boundingBox do
							emit quote
								for j = -boundingBox, boundingBox +1 do
							        for i = -boundingBox, boundingBox + 1 do
								        escape
								        	for k = 0, vectorWidth - 1 do
								        		emit quote
		       										curKernelValTemp[k] = kernelArray[0][(j+boundingBox)*kernelWidth+(i+boundingBox)]*(funcParams[0*numberOfFuncCoef + 0]
								        				+ funcParams[0*numberOfFuncCoef + 1]*x + funcParams[0*numberOfFuncCoef + 2]*y
								        				+ funcParams[0*numberOfFuncCoef + 3]*x*x + funcParams[0*numberOfFuncCoef + 4]*x*y
								        				+ funcParams[0*numberOfFuncCoef + 5]*y*y + funcParams[0*numberOfFuncCoef + 6]*x*x*x
								        				+ funcParams[0*numberOfFuncCoef + 7]*x*x*y + funcParams[0*numberOfFuncCoef + 8]*x*y*y
								        				+ funcParams[0*numberOfFuncCoef + 9]*y*y*y)	
								        		end
		--						        		for l = 1, numberOfBasisKernels-1 do 
								        		for l = 1, luaNumberOfBasisKernels-1 do 
								        			emit quote
								        			curKernelValTemp[k] = curKernelValTemp[k] + kernelArray[l][(j+boundingBox)*kernelWidth+(i+boundingBox)]*(funcParams[l*numberOfFuncCoef + 0]
								        				+ funcParams[l*numberOfFuncCoef + 1]*x + funcParams[l*numberOfFuncCoef + 2]*y
								        				+ funcParams[l*numberOfFuncCoef + 3]*x*x + funcParams[l*numberOfFuncCoef + 4]*x*y
								        				+ funcParams[l*numberOfFuncCoef + 5]*y*y + funcParams[l*numberOfFuncCoef + 6]*x*x*x
								        				+ funcParams[l*numberOfFuncCoef + 7]*x*x*y + funcParams[l*numberOfFuncCoef + 8]*x*y*y
								        				+ funcParams[l*numberOfFuncCoef + 9]*y*y*y)
								        			end
								        		end
								        	end
								        	emit quote
									            curKernelVal = @[&vector(double,vectorWidth)](&curKernelValTemp[0])

									            var curImInVec = @[&vector(float,vectorWidth)](&inputImg[(y+j)*imageWidth + x+i])
									            var curVarInVec = @[&vector(float,vectorWidth)](&inputVar[(y+j)*imageWidth + x+i])
									            var curMaskInVec = @[&vector(uint16,vectorWidth)](&inputMask[(y+j)*imageWidth + x+i])


									            curImOut = curImOut + curImInVec*curKernelVal; 
									            curVarOut = curVarOut + curVarInVec*curKernelVal*curKernelVal;

									            var bitMask : vector(uint16, vectorWidth) = getBitMask(curKernelVal, zeroVec)
									            curMaskOut = curMaskOut or (curMaskInVec and bitMask)

									            curNorm = curNorm + curKernelVal;
									        end
								        end
							    	end
								end
							end
						    emit quote
					    		curImOut = curImOut/curNorm
					    		curVarOut = curVarOut/(curNorm*curNorm)
							end
							for i = 0, vectorWidth-1 do
					    		emit quote
						    		outputImg[y*imageWidth + x+i] = curImOut[i]	
								    outputVar[y*imageWidth + x+i] = curVarOut[i]
					    			outputMask[y*imageWidth + x+i] = curMaskOut[i]
					    		end
					    	end
					    end				
					end
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
	
			C.printf("Output Variance plane, 10x10 box begining at (boundingBox,boundingBox)\n")
			for i=boundingBox,boundingBox+10 do
				for j=boundingBox,boundingBox+10 do
					C.printf("%f\t", outputVar[i*imageWidth + j])
				end
				C.printf("\n")
			end
	
			C.printf("Output Mask plane, 10x10 box begining at (boundingBox,boundingBox)\n")
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
				for j = 0, (boundingBox*2+1) do
					for i = 0, (boundingBox*2+1) do
						C.printf("%f\t", kernelArray[h][j*(boundingBox*2+1) + i])
					end
					C.printf("\n")
				end
				C.printf("\n")
			end



		end
	end
end

    --store the jth coefficient of the ith function in
    --funcParams[i*numberOfFuncCoef + j]
local terra convolveRefactor(inputImg:&float, inputVar:&float, inputMask:&uint16,
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
		local newKernelArray = symbol()
		emit quote
			var [boundingBox] = (kernelWidth - 1)/2

			var t1 = C.clock()

			--refactor kernels
--			var [newKernelArray] : &&float = [&&float](C.malloc(numberOfFuncCoef*C.sizeof(&float)))
			var [newKernelArray] : &&float = [&&float](C.malloc(numberOfFuncCoef*8))
			for h = 0, numberOfFuncCoef do
--				newKernelArray[h] = [&float](C.malloc(kernelWidth*kernelHeight*C.sizeof(float)))
				newKernelArray[h] = [&float](C.malloc(kernelWidth*kernelHeight*4))
				for j = 0, kernelHeight do
					for i = 0, kernelWidth do
						newKernelArray[h][j*kernelWidth+i] = kernelArray[0][j*kernelWidth+i] * funcParams[h]
						for k = 1, numberOfBasisKernels do
							newKernelArray[h][j*kernelWidth+i] = newKernelArray[h][j*kernelWidth+i] 
								+ kernelArray[k][j*kernelWidth+i] * funcParams[k*numberOfFuncCoef + h]
						end
					end
				end
			end


			var zeroVec : vector(float, vectorWidth) = [zerofloatVec(vectorWidth)]
			for y = boundingBox, imageHeight-boundingBox do
				for x = boundingBox, (imageWidth-boundingBox), vectorWidth do
	 				escape
--	 					if(method == "original") then
	 					if(true) then
	 						local curKernelVal = symbol()
	 						local curKernelValTemp = symbol()
	 						local curImOut = symbol()
	 						local curVarOut= symbol()
	 						local curMaskOut = symbol()
	 						local curNorm = symbol()
		 					emit quote
				 				var [curKernelVal] : vector(float,vectorWidth)
								var [curKernelValTemp] : float[vectorWidth]

								var [curImOut] : vector(float,vectorWidth) = [zerofloatVec(vectorWidth)]
				--				var curImOut : &vector(float,vectorWidth) = [&vector(float,vectorWidth)](&outputImg[y*imageWidth + x])
								var [curVarOut] : vector(float,vectorWidth) = [zerofloatVec(vectorWidth)]
								var [curMaskOut] : vector(uint16, vectorWidth) = [zeroUInt16Vec(vectorWidth)]

								var [curNorm] : vector(float,vectorWidth) = [zerofloatVec(vectorWidth)]
				--				var curImOutVec = @[&vector(float,4)](&outputImg[y*imageWidth + x])
							end
--							for j = -boundingBox, boundingBox do
--						        for i = -boundingBox, boundingBox do
							for j = -luaBoundingBox, luaBoundingBox do
								--emit quote
						        for i = -luaBoundingBox, luaBoundingBox do--+ 1 do
						        --escape
						        	for k = 0, vectorWidth - 1 do
						        		emit quote
       										curKernelValTemp[k] = newKernelArray[0][(j+luaBoundingBox)*kernelWidth+(i+luaBoundingBox)]	
						        				+ newKernelArray[1][(j+luaBoundingBox)*kernelWidth+(i+luaBoundingBox)]*([float](x) - [float](imageWidth/2))/[float](imageWidth/2) 
						        				+ newKernelArray[2][(j+luaBoundingBox)*kernelWidth+(i+luaBoundingBox)]*([float](y) - [float](imageHeight/2))/[float](imageHeight/2)
						        				+ newKernelArray[3][(j+luaBoundingBox)*kernelWidth+(i+luaBoundingBox)]*x*x 
						        				+ newKernelArray[4][(j+luaBoundingBox)*kernelWidth+(i+luaBoundingBox)]*x*y
						        				+ newKernelArray[5][(j+luaBoundingBox)*kernelWidth+(i+luaBoundingBox)]*y*y 
						        				+ newKernelArray[6][(j+luaBoundingBox)*kernelWidth+(i+luaBoundingBox)]*x*x*x
						        				+ newKernelArray[7][(j+luaBoundingBox)*kernelWidth+(i+luaBoundingBox)]*x*x*y 
						        				+ newKernelArray[8][(j+luaBoundingBox)*kernelWidth+(i+luaBoundingBox)]*x*y*y
						        				+ newKernelArray[9][(j+luaBoundingBox)*kernelWidth+(i+luaBoundingBox)]*y*y*y
						        		end
						        	end
						        	emit quote
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
						        --end
						    	end
						        --end
						    end
						    emit quote
					    		curImOut = curImOut/curNorm
					    		curVarOut = curVarOut/(curNorm*curNorm)
							end
							for i = 0, vectorWidth-1 do
					    		emit quote
						    		outputImg[y*imageWidth + x+i] = curImOut[i]	
								    outputVar[y*imageWidth + x+i] = curVarOut[i]
					    			outputMask[y*imageWidth + x+i] = curMaskOut[i]
					    		end
					    	end
					    end				
					end
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

			C.printf("Difference Observed!!:\n")


		end
	end
end



local terra convolveWithInterp(inputImg:&float, inputVar:&float, inputMask:&uint16,
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
	 				escape
--	 					if(method == "original") then
	 					if(true) then
	 						local curKernelVal = symbol()
	 						local curKernelValTemp = symbol()
	 						local curImOut = symbol()
	 						local curVarOut= symbol()
	 						local curMaskOut = symbol()
	 						local curNorm = symbol()
		 					emit quote
				 				var [curKernelVal] : vector(float,vectorWidth)
								var [curKernelValTemp] : float[vectorWidth]

								var [curImOut] : vector(float,vectorWidth) = [zerofloatVec(vectorWidth)]
				--				var curImOut : &vector(float,vectorWidth) = [&vector(float,vectorWidth)](&outputImg[y*imageWidth + x])
								var [curVarOut] : vector(float,vectorWidth) = [zerofloatVec(vectorWidth)]
								var [curMaskOut] : vector(uint16, vectorWidth) = [zeroUInt16Vec(vectorWidth)]

								var [curNorm] : vector(float,vectorWidth) = [zerofloatVec(vectorWidth)]
				--				var curImOutVec = @[&vector(float,4)](&outputImg[y*imageWidth + x])
							end
--							for j = -boundingBox, boundingBox do
--						        for i = -boundingBox, boundingBox do
							for j = -luaBoundingBox, luaBoundingBox do
								emit quote
						        for i = -luaBoundingBox, luaBoundingBox + 1 do
						        escape
						        	for k = 0, vectorWidth - 1 do
						        		emit quote
       										curKernelValTemp[k] = kernelArray[0][(j+luaBoundingBox)*kernelWidth+(i+luaBoundingBox)]*(funcParams[0*numberOfFuncCoef + 0]
						        				+ funcParams[0*numberOfFuncCoef + 1]*x + funcParams[0*numberOfFuncCoef + 2]*y
						        				+ funcParams[0*numberOfFuncCoef + 3]*x*x + funcParams[0*numberOfFuncCoef + 4]*x*y
						        				+ funcParams[0*numberOfFuncCoef + 5]*y*y + funcParams[0*numberOfFuncCoef + 6]*x*x*x
						        				+ funcParams[0*numberOfFuncCoef + 7]*x*x*y + funcParams[0*numberOfFuncCoef + 8]*x*y*y
						        				+ funcParams[0*numberOfFuncCoef + 9]*y*y*y)	
						        		end
--						        		for l = 1, numberOfBasisKernels-1 do 
						        		for l = 1, luaNumberOfBasisKernels-1 do 
						        			emit quote
						        			curKernelValTemp[k] = curKernelValTemp[k] + kernelArray[l][(j+luaBoundingBox)*kernelWidth+(i+luaBoundingBox)]*(funcParams[l*numberOfFuncCoef + 0]
						        				+ funcParams[l*numberOfFuncCoef + 1]*x + funcParams[l*numberOfFuncCoef + 2]*y
						        				+ funcParams[l*numberOfFuncCoef + 3]*x*x + funcParams[l*numberOfFuncCoef + 4]*x*y
						        				+ funcParams[l*numberOfFuncCoef + 5]*y*y + funcParams[l*numberOfFuncCoef + 6]*x*x*x
						        				+ funcParams[l*numberOfFuncCoef + 7]*x*x*y + funcParams[l*numberOfFuncCoef + 8]*x*y*y
						        				+ funcParams[l*numberOfFuncCoef + 9]*y*y*y)
						        			end
						        		end
						        	end
						        	emit quote
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
						        end
						    	end
						        end
						    end
						    emit quote
					    		curImOut = curImOut/curNorm
					    		curVarOut = curVarOut/(curNorm*curNorm)
							end
							for i = 0, vectorWidth-1 do
					    		emit quote
						    		outputImg[y*imageWidth + x+i] = curImOut[i]	
								    outputVar[y*imageWidth + x+i] = curVarOut[i]
					    			outputMask[y*imageWidth + x+i] = curMaskOut[i]
					    		end
					    	end
					    end				
					end
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


--    --store the jth coefficient of the ith function in
--    --funcParams[i*numberOfFuncCoef + j]
----refactor the kernel if there are more basis kernels than polynomial coefficients
----unroll convolution at one point by metaprogramming kernel values
--local kernelValSymbols = {}
--
--local terra convolveFinal(inputImg:&float, inputVar:&float, inputMask:&uint16,
--							outputImg:&float, outputVar:&float, outputMask:&uint16,
--                    		imageWidth:int, imageHeight:int, kernelArray:&&float, 
--                    		funcParams:&float, numberOfBasisKernels:int,
--							kernelWidth:int, kernelHeight:int, numberOfFuncCoef:int)
--	C.printf("Begin Terra function\n")
--	C.printf("imageWidth = %d, imageHeight = %d, numberOfBasisKernels = %d, kernelWidth = %d, kernelHeight = %d, numberOfFuncCoef = %d\n", imageWidth,
--		imageHeight, numberOfBasisKernels, kernelWidth, kernelHeight, numberOfFuncCoef)
--	escape
--		local newKernelArray = symbol()
--		emit quote
--
--			var t1 = C.clock()
--
--			--refactor kernels
----			var [newKernelArray] : &&float = [&&float](C.malloc(numberOfFuncCoef*C.sizeof(&float)))
--			var [newKernelArray] : &&float = [&&float](C.malloc(numberOfFuncCoef*8))
--			for h = 0, numberOfFuncCoef do
----				newKernelArray[h] = [&float](C.malloc(kernelWidth*kernelHeight*C.sizeof(float)))
--				newKernelArray[h] = [&float](C.malloc(kernelWidth*kernelHeight*4))
--				for j = 0, kernelHeight do
--					for i = 0, kernelWidth do
--						newKernelArray[h][j*kernelWidth+i] = kernelArray[0][j*kernelWidth+i] * funcParams[h]
--						for k = 1, numberOfBasisKernels do
--							newKernelArray[h][j*kernelWidth+i] = newKernelArray[h][j*kernelWidth+i] 
--								+ kernelArray[k][j*kernelWidth+i] * funcParams[k*numberOfFuncCoef + h]
--						end
--					end
--				end
--			end
--		
--			escape
--		
--		        for i = 1, luaKernelArea*5 do
--		            local cur_kernelValue_symbol = symbol(float, "kVal"..i)
--		            table.insert(kernelValSymbols, cur_kernelValue_symbol)
--		        end
--		
--		
--				for j = -boundingBox, boundingBox do
--					for i = -boundingBox, boundingBox do
--						emit quote
--							var [kernelValSymbols[(kernelWidth*(j+boundingBox) + (i+boundingBox))*5+1]] = [kernel1(i, j)]
--							var [kernelValSymbols[(kernelWidth*(j+boundingBox) + (i+boundingBox))*5+2]] = [kernel2(i, j)]
--							var [kernelValSymbols[(kernelWidth*(j+boundingBox) + (i+boundingBox))*5+3]] = [kernel3(i, j)]
--							var [kernelValSymbols[(kernelWidth*(j+boundingBox) + (i+boundingBox))*5+4]] = [kernel4(i, j)]
--							var [kernelValSymbols[(kernelWidth*(j+boundingBox) + (i+boundingBox))*5+5]] = [kernel5(i, j)]
--						end
--					end
--				end
--			end
--		
--			var tB = C.clock()
--			C.printf("time spent calculating kernel = %f ms!!!!!!!!!!!!!!@@@@@@@@@@@@", (float)(tB-t1)/C.CLOCKS_PER_SEC*1000.0f)
--		
--		
--			for y = boundingBox, imageHeight-boundingBox do
--				for x = boundingBox, (imageWidth-boundingBox) do
--					var curKernelVal : float
--	
--					var curImOut : float = 0.0f
--					var curVarOut : float = 0.0f
--					var curMaskOut : uint16 = 0
--	
--					var curNorm : float = 0.0f
--	
--					var kernelArea = (boundingBox*2+1)*(boundingBox*2+1)
--					var kernelWidth = boundingBox*2+1		
--	
--					escape 	
--					    for j = -boundingBox, boundingBox do
--					        for i = -boundingBox, boundingBox do
--					        	emit quote
--					        		curKernelVal = polynomial1(x, y)*[kernelValSymbols[(kernelWidth*(j+boundingBox) + (i+boundingBox))*5+1]] +
--						                polynomial2(x, y)*[kernelValSymbols[(kernelWidth*(j+boundingBox) + (i+boundingBox))*5+2]] + polynomial3(x, y)*[kernelValSymbols[(kernelWidth*(j+boundingBox) + (i+boundingBox))*5+3]] + 
--						                polynomial4(x, y)*[kernelValSymbols[(kernelWidth*(j+boundingBox) + (i+boundingBox))*5+4]] + polynomial5(x, y)*[kernelValSymbols[(kernelWidth*(j+boundingBox) + (i+boundingBox))*5+5]];
--	
--						            var curImIn = inputImg[(y+j)*imageWidth + x+i]
--						            var curVarIn = inputVar[(y+j)*imageWidth + x+i]
--						            var curMaskIn = inputMask[(y+j)*imageWidth + x+i]
--	
--	
--						            curImOut = curImOut + curImIn*curKernelVal; 
--						            curVarOut = curVarOut + curVarIn*curKernelVal*curKernelVal;
--	
--						            if curKernelVal ~= 0.0f then
--						            	curMaskOut = curMaskOut or curMaskIn 
--						            end
--	
--						            curNorm = curNorm + curKernelVal;
--						        end
--					        end
--					    end
--					end
--				    curImOut = curImOut/curNorm
--				    outputImg[y*imageWidth + x] = curImOut
--	
--				    curVarOut = curVarOut/(curNorm*curNorm)
--				    outputVar[y*imageWidth + x] = curVarOut
--	
--				    outputMask[y*imageWidth + x] = curMaskOut
--	
--				end
--			end
--
--			var t2 = C.clock()
--			C.printf("\n\nVectorized %d wide masked image lin combo %dx%d blur Terra:\n", vectorWidth, 1+2*boundingBox, 1+2*boundingBox)
--			C.printf("outputImg[boundingBox*imageWidth + boundingBox] = %f, computation took: %f ms\n", outputImg[boundingBox*imageWidth + boundingBox],  (t2-t1)/1000.0)
--			C.printf("C.CLOCKS_PER_SEC = %d\n", C.CLOCKS_PER_SEC)
--	
--			C.printf("Input image plane, 10x10 box begining at (boundingBox,boundingBox)\n")
--			for i=boundingBox,boundingBox+10 do
--				for j=boundingBox,boundingBox+10 do
--					C.printf("%f\t", inputImg[i*imageWidth + j])
--				end
--				C.printf("\n")
--			end
--
--			C.printf("Output image plane, 10x10 box begining at (boundingBox,boundingBox)\n")
--			for i=boundingBox,boundingBox+10 do
--				for j=boundingBox,boundingBox+10 do
--					C.printf("%f\t", outputImg[i*imageWidth + j])
--				end
--				C.printf("\n")
--			end
--	
--			C.printf("Variance plane, 10x10 box begining at (boundingBox,boundingBox)\n")
--			for i=boundingBox,boundingBox+10 do
--				for j=boundingBox,boundingBox+10 do
--					C.printf("%f\t", outputVar[i*imageWidth + j])
--				end
--				C.printf("\n")
--			end
--	
--			C.printf("Mask plane, 10x10 box begining at (boundingBox,boundingBox)\n")
--			for i=boundingBox,boundingBox+10 do
--				for j=boundingBox,boundingBox+10 do
--					C.printf("%d\t", outputMask[i*imageWidth + j])
--				end
--				C.printf("\n")
--			end
--
--			C.printf("Function coefficients:\n")
--			for i = 0, numberOfBasisKernels do
--				for j = 0, numberOfFuncCoef do
--					C.printf("%f\t", funcParams[i*numberOfFuncCoef + j])
--				end
--				C.printf("\n")
--			end
--
--			C.printf("Basis kernels:\n")
--			for h = 0, numberOfBasisKernels do
--				for j = 0, (luaBoundingBox*2+1) do
--					for i = 0, (luaBoundingBox*2+1) do
--						C.printf("%f\t", kernelArray[h][j*(luaBoundingBox*2+1) + i])
--					end
--					C.printf("\n")
--				end
--				C.printf("\n")
--			end
--		end
--	end
--end



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

--convolve:printpretty()

--the first term is the name in the .o and .h files, the second name is the name in this file

local kerneSize3 = genConvolve(27)


--local kerneSize5 = genConvolve(5)
--local kerneSize7 = genConvolve(7)
--local kerneSize9 = genConvolve(9)
--local kerneSize11 = genConvolve(11)
--local kerneSize13 = genConvolve(13)
--local kerneSize15 = genConvolve(15)
--local kerneSize17 = genConvolve(17)
--local kerneSize19 = genConvolve(19)
--local kerneSize21 = genConvolve(21)
--local kerneSize23 = genConvolve(23)
--local kerneSize25 = genConvolve(25)
--local kerneSize27 = genConvolve(27)


terralib.saveobj("blurExampleTerraStandalone.o",{ terraFuncNameInC3  = convolveDoublePrecisionKernel})
--terralib.saveobj("blurExampleTerraStandalone5.o",{ terraFuncNameInC5  = kerneSize5})
--terralib.saveobj("blurExampleTerraStandalone.o",{ terraFuncNameInC3 = kerneSize3,
--												  terraFuncNameInC5 = kerneSize5})


--terralib.saveobj("blurExampleTerraStandalone.o",{ terraFuncNameInC3  = kerneSize3})
--terralib.saveobj("blurExampleTerraStandalone.o",{ terraFuncNameInC5  = kerneSize5})
--terralib.saveobj("blurExampleTerraStandalone.o",{ terraFuncNameInC7  = kerneSize7})
--terralib.saveobj("blurExampleTerraStandalone.o",{ terraFuncNameInC9  = kerneSize9})
--terralib.saveobj("blurExampleTerraStandalone.o",{ terraFuncNameInC11 = kerneSize11})
--terralib.saveobj("blurExampleTerraStandalone.o",{ terraFuncNameInC13 = kerneSize13})
--terralib.saveobj("blurExampleTerraStandalone.o",{ terraFuncNameInC15 = kerneSize15})
--terralib.saveobj("blurExampleTerraStandalone.o",{ terraFuncNameInC17 = kerneSize17})
--terralib.saveobj("blurExampleTerraStandalone.o",{ terraFuncNameInC19 = kerneSize19})
--terralib.saveobj("blurExampleTerraStandalone.o",{ terraFuncNameInC21 = kerneSize21})
--terralib.saveobj("blurExampleTerraStandalone.o",{ terraFuncNameInC23 = kerneSize23})
--terralib.saveobj("blurExampleTerraStandalone.o",{ terraFuncNameInC25 = kerneSize25})
--terralib.saveobj("blurExampleTerraStandalone.o",{ terraFuncNameInC27 = kerneSize27})




