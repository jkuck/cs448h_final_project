local C = terralib.includecstring [[
    #include<stdio.h>
    #include<stdlib.h>
    #include<time.h>
    #include<string.h>
    #include<math.h>

]]

local boundingBox = 9

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


local function kernel1(x, y)
    local sigmaX1 = 2.0;
    local sigmaY1 = 2.0;
    local theta1 = 0.0; --//rotation of sigmaX axis
    return (math.exp(-((x*math.cos(theta1) +y*math.sin(theta1))*(x*math.cos(theta1) +y*math.sin(theta1)))
                    /(2*sigmaX1*sigmaX1)) / (math.sqrt(2*math.pi)*sigmaX1))
                    *(math.exp(-((y*math.cos(theta1) - x*math.sin(theta1))*(y*math.cos(theta1) - x*math.sin(theta1)))
                    /(2*sigmaY1*sigmaY1)) / (math.sqrt(2*math.pi)*sigmaY1));
end



local function kernel2(x, y)
    local sigmaX2 = 0.5;
    local sigmaY2 = 4.0;
    local theta2 = 0.0; --//rotation of sigmaX axis
    return (math.exp(-((x*math.cos(theta2) +y*math.sin(theta2))*(x*math.cos(theta2) +y*math.sin(theta2)))
                    /(2*sigmaX2*sigmaX2)) / (math.sqrt(2*math.pi)*sigmaX2))
                    *(math.exp(-((y*math.cos(theta2) - x*math.sin(theta2))*(y*math.cos(theta2) - x*math.sin(theta2)))
                    /(2*sigmaY2*sigmaY2)) / (math.sqrt(2*math.pi)*sigmaY2));
end

local function kernel3(x, y)
    local sigmaX3 = 0.5;
    local sigmaY3 = 4.0;
    local theta3 = 3.14159/4; --//rotation of sigmaX axis
    return (math.exp(-((x*math.cos(theta3) +y*math.sin(theta3))*(x*math.cos(theta3) +y*math.sin(theta3)))
                    /(2*sigmaX3*sigmaX3)) / (math.sqrt(2*math.pi)*sigmaX3))
                    *(math.exp(-((y*math.cos(theta3) - x*math.sin(theta3))*(y*math.cos(theta3) - x*math.sin(theta3)))
                    /(2*sigmaY3*sigmaY3)) / (math.sqrt(2*math.pi)*sigmaY3));
end

local function kernel4(x, y)
    local sigmaX4 = 0.5;
    local sigmaY4 = 4.0;
    local theta4 = 3.14159/2; --//rotation of sigmaX axis
    return (math.exp(-((x*math.cos(theta4) +y*math.sin(theta4))*(x*math.cos(theta4) +y*math.sin(theta4)))
                    /(2*sigmaX4*sigmaX4)) / (math.sqrt(2*math.pi)*sigmaX4))
                    *(math.exp(-((y*math.cos(theta4) - x*math.sin(theta4))*(y*math.cos(theta4) - x*math.sin(theta4)))
                    /(2*sigmaY4*sigmaY4)) / (math.sqrt(2*math.pi)*sigmaY4));
end


local function kernel5(x, y)
    local sigmaX5 = 4.0;
    local sigmaY5 = 4.0;
    local theta5 = 0.0; --//rotation of sigmaX axis
    return (math.exp(-((x*math.cos(theta5) +y*math.sin(theta5))*(x*math.cos(theta5) +y*math.sin(theta5)))
                    /(2*sigmaX5*sigmaX5)) / (math.sqrt(2*math.pi)*sigmaX5))
                    *(math.exp(-((y*math.cos(theta5) - x*math.sin(theta5))*(y*math.cos(theta5) - x*math.sin(theta5)))
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


local function vecnf(n)
	local zeros = {}
	for i=1,n do zeros[i] = `0.f end
	return `vector(zeros)
end



terra vectormask8(a : vector(float,8), b : vector(float,8))
    return terralib.select(a ~= b, [vector(uint16,8)](0xFFFFULL),[vector(uint16,8)](0) ) 
end

local function metaProgramBlurMaskedImageLinComboVectorize8()
	return quote
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

		for i = 0, width*height do
			inputIm[i] = [float](i + 1)
			inputVar[i] = [float](i + 2)
			inputMask[i] = [uint16](i%width + i/width)


			outputIm[i] = 0.0f
			outputVar[i] = 0.0f
			outputMask[i] = 0
		end

		var t1 = C.clock()

		for y = boundingBox, height-boundingBox do
			for x = boundingBox, (width-boundingBox), 8 do
				--var curOut : float = 0
				--var curNorm : float = 0
				var curKernelVal : vector(float,8)
				var curKernelValTemp : float[8]
				var curImOut : vector(float,8) = vector(0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f)
				var curVarOut : vector(float,8) = vector(0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f)

	--			var curMaskOut : vector(uint16,4) = vector(0,0,0,0)

				var curMaskOut : vector(uint16, 8) = vector(0, 0, 0, 0,0, 0, 0, 0)
				var zeroVec : vector(float, 8) = vector(0.f, 0.f, 0.f, 0.f,0.f,0.f,0.f,0.f)

				var curNorm : vector(float,8) = vector(0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f)
	--			var curImOutVec = @[&vector(float,4)](&outputIm[y*width + x])

				var p1 : float = polynomial1(x, y);
				var p2 : float = polynomial2(x, y);
				var p3 : float = polynomial3(x, y);
				var p4 : float = polynomial4(x, y);
				var p5 : float = polynomial5(x, y);
				 
				escape
	--				local vectorWidth = 8
	--				local pol
	--		        for i = 1, vectorWidth do
	--		            local cur_func_symbol = symbol(&float, "function"..i)
	--		            table.insert(function_symbols, cur_func_symbol)
	--		        end



					for j = -boundingBox, boundingBox do
				        for i = -boundingBox, boundingBox do
				        	emit quote
				        		curKernelValTemp[0] = polynomial1(x, y)*[kernel1(i, j)] +
					                polynomial2(x, y)*[kernel2(i, j)] + polynomial3(x, y)*[kernel3(i, j)] + 
					                polynomial4(x, y)*[kernel4(i, j)] + polynomial5(x, y)*[kernel5(i, j)];
				        		curKernelValTemp[1] = polynomial1(x+1, y)*[kernel1(i, j)] +
					                polynomial2(x+1, y)*[kernel2(i, j)] + polynomial3(x+1, y)*[kernel3(i, j)] + 
					                polynomial4(x+1, y)*[kernel4(i, j)] + polynomial5(x+1, y)*[kernel5(i, j)];
				        		curKernelValTemp[2] = polynomial1(x+2, y)*[kernel1(i, j)] +
					                polynomial2(x+2, y)*[kernel2(i, j)] + polynomial3(x+2, y)*[kernel3(i, j)] + 
					                polynomial4(x+2, y)*[kernel4(i, j)] + polynomial5(x+2, y)*[kernel5(i, j)];
				        		curKernelValTemp[3] = polynomial1(x+3, y)*[kernel1(i, j)] +
					                polynomial2(x+3, y)*[kernel2(i, j)] + polynomial3(x+3, y)*[kernel3(i, j)] + 
					                polynomial4(x+3, y)*[kernel4(i, j)] + polynomial5(x+3, y)*[kernel5(i, j)];
				        		curKernelValTemp[4] = polynomial1(x+4, y)*[kernel1(i, j)] +
					                polynomial2(x+4, y)*[kernel2(i, j)] + polynomial3(x+4, y)*[kernel3(i, j)] + 
					                polynomial4(x+4, y)*[kernel4(i, j)] + polynomial5(x+4, y)*[kernel5(i, j)];		
				        		curKernelValTemp[5] = polynomial1(x+5, y)*[kernel1(i, j)] +
					                polynomial2(x+5, y)*[kernel2(i, j)] + polynomial3(x+5, y)*[kernel3(i, j)] + 
					                polynomial4(x+5, y)*[kernel4(i, j)] + polynomial5(x+5, y)*[kernel5(i, j)];
				        		curKernelValTemp[6] = polynomial1(x+6, y)*[kernel1(i, j)] +
					                polynomial2(x+6, y)*[kernel2(i, j)] + polynomial3(x+6, y)*[kernel3(i, j)] + 
					                polynomial4(x+6, y)*[kernel4(i, j)] + polynomial5(x+6, y)*[kernel5(i, j)];
				        		curKernelValTemp[7] = polynomial1(x+7, y)*[kernel1(i, j)] +
					                polynomial2(x+7, y)*[kernel2(i, j)] + polynomial3(x+7, y)*[kernel3(i, j)] + 
					                polynomial4(x+7, y)*[kernel4(i, j)] + polynomial5(x+7, y)*[kernel5(i, j)];
					            
					            curKernelVal = @[&vector(float,8)](&curKernelValTemp[0])

					            var curImInVec = @[&vector(float,8)](&inputIm[(y+j)*width + x+i])
					            var curVarInVec = @[&vector(float,8)](&inputVar[(y+j)*width + x+i])
					            var curMaskInVec = @[&vector(uint16,8)](&inputMask[(y+j)*width + x+i])


					            curImOut = curImOut + curImInVec*curKernelVal; 
					            curVarOut = curVarOut + curVarInVec*curKernelVal*curKernelVal;

					            var bitMask : vector(uint16, 8) = vectormask8(curKernelVal, zeroVec)
					            curMaskOut = curMaskOut or (curMaskInVec and bitMask)

					            curNorm = curNorm + curKernelVal;
					        end
				        end
				    end
				end
			    curImOut = curImOut/curNorm
			    outputIm[y*width + x] = curImOut[0]
			    outputIm[y*width + x+1] = curImOut[1]
			    outputIm[y*width + x+2] = curImOut[2]
			    outputIm[y*width + x+3] = curImOut[3]
			    outputIm[y*width + x+4] = curImOut[4]
			    outputIm[y*width + x+5] = curImOut[5]
			    outputIm[y*width + x+6] = curImOut[6]
			    outputIm[y*width + x+7] = curImOut[7]

--				@[&vector(float,8)](&outputIm[y*width + x]) = curImOut
--				outputIm[y*width + x] = curImOut
	--		    outputIm[x + width*y] = curImOut/curNorm

			    curVarOut = curVarOut/(curNorm*curNorm)
			    outputVar[y*width + x] = curVarOut[0]
			    outputVar[y*width + x+1] = curVarOut[1]
			    outputVar[y*width + x+2] = curVarOut[2]
			    outputVar[y*width + x+3] = curVarOut[3]
			    outputVar[y*width + x+4] = curVarOut[4]
			    outputVar[y*width + x+5] = curVarOut[5]
			    outputVar[y*width + x+6] = curVarOut[6]
			    outputVar[y*width + x+7] = curVarOut[7]

			    outputMask[y*width + x] = curMaskOut[0]
			    outputMask[y*width + x+1] = curMaskOut[1]
			    outputMask[y*width + x+2] = curMaskOut[2]
			    outputMask[y*width + x+3] = curMaskOut[3]
			    outputMask[y*width + x+4] = curMaskOut[4]
			    outputMask[y*width + x+5] = curMaskOut[5]
			    outputMask[y*width + x+6] = curMaskOut[6]
			    outputMask[y*width + x+7] = curMaskOut[7]

			end
		end

		var t2 = C.clock()
		C.printf("\n\nVectorized 8wide masked image lin combo %dx%d blur Terra:\n", boundingBox, boundingBox)
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
	end
end

local terra blurMaskedImageLinComboVectorize8()
	[metaProgramBlurMaskedImageLinComboVectorize8()]
end


blurMaskedImageLinComboVectorize8()