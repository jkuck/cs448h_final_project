local C = terralib.includecstring [[
    #include<stdio.h>
    #include<stdlib.h>
    #include<time.h>
    #include<string.h>
    #include<math.h>

]]

local boundingBox = 2

--local terra loadFitsImage()
--    var im = lsst::afw::image::MaskedImage<float>("./images/calexp-004207-g3-0123.fits");
--    var width = im.getWidth(), height = im.getHeight();
--end

local terra blurArrayTerra()
--	var sum : float = 0.0f
	var width : int = 1280
	var height : int = 874

	var input : &float
	input = [&float](C.malloc(sizeof(float)*width*height))

	var intermediate : &float
	intermediate = [&float](C.malloc(sizeof(float)*width*height))

	var output : &float
	output = [&float](C.malloc(sizeof(float)*width*height))

	for i = 0, width*height do
		input[i] = [float](i + 1)
		output[i] = 0.0f
		intermediate[i] = 0.0f
	end

	var t1 = C.clock()

	for y = 0, height do
		for x=2, width-2 do

			intermediate[y*width + x] = (input[y*width + (x-2)] + input[y*width + (x-1)] + input[y*width + x]
									  + input[y*width + (x+1)] + input[y*width + (x+2)])/5.0f

		end
	end

	for y = 2, height-2 do
		for x=2, width-2 do

			output[y*width + x] = (intermediate[(y-2)*width + (x)] + intermediate[(y-1)*width + (x)] + intermediate[y*width + x]
									  + intermediate[(y+1)*width + x] + intermediate[(y+2)*width + x])/5.0f

		end
	end

	var t2 = C.clock()
	C.printf("\n\nseparable %dx%d blur Terra:\n", 2, 2)
	C.printf("output[2*width + 2] = %f, computation took: %d micro-seconds\n", output[2*width + 2],  (t2-t1))

	C.printf("C.CLOCKS_PER_SEC = %d\n", C.CLOCKS_PER_SEC)

	C.free(input)
	C.free(intermediate)
	C.free(output)

end



local terra blurArrayTerraVectorizedAcrossX()
	var sum : float = 0.0f
	var input : &float
	input = [&float](C.malloc(10000000*sizeof(float)))
	var output : &float
	output = [&float](C.malloc(10000000*sizeof(float)))

	for i = 0, 10000000 do
		input[i] = [float](i + 1)
		output[i] = 0.0f
	end

	var t1 = C.clock()

--	var curInVec2 : vector(float,4)

	for i = 0, (1249992)do
		var curOutVec = @[&vector(float,4)](&output[i*8])
		var curOutVec2 = @[&vector(float,4)](&output[i*8+4])

		for j=0,64 do
			var curAddVec = @[&vector(float,4)](&input[i*8+j])
			var curAddVec2 = @[&vector(float,4)](&input[i*8+4+j])

			curOutVec = curOutVec + curAddVec
			curOutVec2 = curOutVec2 + curAddVec2
		end
			output[i*8] = curOutVec[0]
			output[i*8+1] = curOutVec[1]
			output[i*8+2] = curOutVec[2]
			output[i*8+3] = curOutVec[3]

			output[i*8+4] = curOutVec2[0]
			output[i*8+5] = curOutVec2[1]
			output[i*8+6] = curOutVec2[2]
			output[i*8+7] = curOutVec2[3]

	end

	var t2 = C.clock()
	C.printf("\n\nblurArrayTerraVectorizedAcrossX:\n")
	C.printf("output[10] = %f, computation took: %d micro-seconds\n", output[10],  (t2-t1))

	C.printf("C.CLOCKS_PER_SEC = %d\n", C.CLOCKS_PER_SEC)

	C.free(input)
	C.free(output)

end




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



local terra vecPolynomial1(xVec:vector(float,8), y:float)
	var yVec = vector(y, y, y, y, y, y, y, y)
 	return [vector(float, 8)](0.1f) + [vector(float, 8)](0.002f)*xVec + [vector(float, 8)](0.003f)*yVec + [vector(float, 8)](0.4f)*xVec*xVec + [vector(float, 8)](0.5f)*xVec*yVec
                     + [vector(float, 8)](0.6f)*yVec*yVec + [vector(float, 8)]( 0.0007f)*xVec*xVec*xVec + [vector(float, 8)](0.0008f)*xVec*xVec*yVec + [vector(float, 8)](0.0009f)*xVec*yVec*yVec
                     + [vector(float, 8)](0.00011f)*yVec*yVec*yVec;
end

local terra vecPolynomial2(xVec:vector(float,8), y:float)
	var yVec = vector(y, y, y, y, y, y, y, y)
	return [vector(float, 8)](1.1f) + [vector(float, 8)](1.002f)*xVec + [vector(float, 8)](1.003f)*yVec + [vector(float, 8)](1.4f)*xVec*xVec + [vector(float, 8)](1.5f)*xVec*yVec
                     + [vector(float, 8)](1.6f)*yVec*yVec + [vector(float, 8)]( 1.0007f)*xVec*xVec*xVec + [vector(float, 8)](1.0008f)*xVec*xVec*yVec + [vector(float, 8)](1.0009f)*xVec*yVec*yVec
                     + [vector(float, 8)](1.00011f)*yVec*yVec*yVec;
end

local terra vecPolynomial3(xVec:vector(float,8), y:float)
	var yVec = vector(y, y, y, y, y, y, y, y)
	return [vector(float, 8)](2.1f) + [vector(float, 8)](2.002f)*xVec + [vector(float, 8)](2.003f)*yVec + [vector(float, 8)](2.4f)*xVec*xVec + [vector(float, 8)](2.5f)*xVec*yVec
                     + [vector(float, 8)](2.6f)*yVec*yVec + [vector(float, 8)]( 2.0007f)*xVec*xVec*xVec + [vector(float, 8)](2.0008f)*xVec*xVec*yVec + [vector(float, 8)](2.0009f)*xVec*yVec*yVec
                     + [vector(float, 8)](2.00011f)*yVec*yVec*yVec;
end

local terra vecPolynomial4(xVec:vector(float,8), y:float)
	var yVec = vector(y, y, y, y, y, y, y, y)
	return [vector(float, 8)](3.1f) + [vector(float, 8)](3.002f)*xVec + [vector(float, 8)](3.003f)*yVec + [vector(float, 8)](3.4f)*xVec*xVec + [vector(float, 8)](3.5f)*xVec*yVec
                     + [vector(float, 8)](3.6f)*yVec*yVec + [vector(float, 8)]( 3.0007f)*xVec*xVec*xVec + [vector(float, 8)](3.0008f)*xVec*xVec*yVec + [vector(float, 8)](3.0009f)*xVec*yVec*yVec
                     + [vector(float, 8)](3.00011f)*yVec*yVec*yVec;
end

local terra vecPolynomial5(xVec:vector(float,8), y:float)
	var yVec = vector(y, y, y, y, y, y, y, y)
	return [vector(float, 8)](4.1f) + [vector(float, 8)](4.002f)*xVec + [vector(float, 8)](4.003f)*yVec + [vector(float, 8)](4.4f)*xVec*xVec + [vector(float, 8)](4.5f)*xVec*yVec
                     + [vector(float, 8)](4.6f)*yVec*yVec + [vector(float, 8)]( 4.0007f)*xVec*xVec*xVec + [vector(float, 8)](4.0008f)*xVec*xVec*yVec + [vector(float, 8)](4.0009f)*xVec*yVec*yVec
                     + [vector(float, 8)](4.00011f)*yVec*yVec*yVec;
end


local terra kernel1(x:float, y:float)
    var sigmaX1 : float = 2.0f;
    var sigmaY1 : float = 2.0f;
    var theta1 : float = 0.0f; --//rotation of sigmaX axis
    return (C.exp(-((x*C.cos(theta1) +y*C.sin(theta1))*(x*C.cos(theta1) +y*C.sin(theta1)))
                    /(2*sigmaX1*sigmaX1)) / (C.sqrtf(2*C.M_PI)*sigmaX1))
                    *(C.exp(-((y*C.cos(theta1) - x*C.sin(theta1))*(y*C.cos(theta1) - x*C.sin(theta1)))
                    /(2*sigmaY1*sigmaY1)) / (C.sqrtf(2*C.M_PI)*sigmaY1));
end



local terra kernel2(x:float, y:float)
    var sigmaX2 : float = 0.5f;
    var sigmaY2 : float = 4.0f;
    var theta2 : float = 0.0f; --//rotation of sigmaX axis
    return (C.exp(-((x*C.cos(theta2) +y*C.sin(theta2))*(x*C.cos(theta2) +y*C.sin(theta2)))
                    /(2*sigmaX2*sigmaX2)) / (C.sqrtf(2*C.M_PI)*sigmaX2))
                    *(C.exp(-((y*C.cos(theta2) - x*C.sin(theta2))*(y*C.cos(theta2) - x*C.sin(theta2)))
                    /(2*sigmaY2*sigmaY2)) / (C.sqrtf(2*C.M_PI)*sigmaY2));
end

local terra kernel3(x:float, y:float)
    var sigmaX3 : float = 0.5f;
    var sigmaY3 : float = 4.0f;
    var theta3 : float = 3.14159f/4; --//rotation of sigmaX axis
    return (C.exp(-((x*C.cos(theta3) +y*C.sin(theta3))*(x*C.cos(theta3) +y*C.sin(theta3)))
                    /(2*sigmaX3*sigmaX3)) / (C.sqrtf(2*C.M_PI)*sigmaX3))
                    *(C.exp(-((y*C.cos(theta3) - x*C.sin(theta3))*(y*C.cos(theta3) - x*C.sin(theta3)))
                    /(2*sigmaY3*sigmaY3)) / (C.sqrtf(2*C.M_PI)*sigmaY3));
end

local terra kernel4(x:float, y:float)
    var sigmaX4 : float = 0.5f;
    var sigmaY4 : float = 4.0f;
    var theta4 : float = 3.14159f/2; --//rotation of sigmaX axis
    return (C.exp(-((x*C.cos(theta4) +y*C.sin(theta4))*(x*C.cos(theta4) +y*C.sin(theta4)))
                    /(2*sigmaX4*sigmaX4)) / (C.sqrtf(2*C.M_PI)*sigmaX4))
                    *(C.exp(-((y*C.cos(theta4) - x*C.sin(theta4))*(y*C.cos(theta4) - x*C.sin(theta4)))
                    /(2*sigmaY4*sigmaY4)) / (C.sqrtf(2*C.M_PI)*sigmaY4));
end


local terra kernel5(x:float, y:float)
    var sigmaX5 : float = 4.0f;
    var sigmaY5 : float = 4.0f;
    var theta5 : float = 0.0; --//rotation of sigmaX axis
    return (C.exp(-((x*C.cos(theta5) +y*C.sin(theta5))*(x*C.cos(theta5) +y*C.sin(theta5)))
                    /(2*sigmaX5*sigmaX5)) / (C.sqrtf(2*C.M_PI)*sigmaX5))
                    *(C.exp(-((y*C.cos(theta5) - x*C.sin(theta5))*(y*C.cos(theta5) - x*C.sin(theta5)))
                    /(2*sigmaY5*sigmaY5)) / (C.sqrtf(2*C.M_PI)*sigmaY5));
end

local function luaKernel1(x, y)
    local sigmaX1 = 2.0;
    local sigmaY1 = 2.0;
    local theta1 = 0.0; --//rotation of sigmaX axis
    return (math.exp(-((x*math.cos(theta1) +y*math.sin(theta1))*(x*math.cos(theta1) +y*math.sin(theta1)))
                    /(2*sigmaX1*sigmaX1)) / (math.sqrt(2*math.pi)*sigmaX1))
                    *(math.exp(-((y*math.cos(theta1) - x*math.sin(theta1))*(y*math.cos(theta1) - x*math.sin(theta1)))
                    /(2*sigmaY1*sigmaY1)) / (math.sqrt(2*math.pi)*sigmaY1));
end



local function luaKernel2(x, y)
    local sigmaX2 = 0.5;
    local sigmaY2 = 4.0;
    local theta2 = 0.0; --//rotation of sigmaX axis
    return (math.exp(-((x*math.cos(theta2) +y*math.sin(theta2))*(x*math.cos(theta2) +y*math.sin(theta2)))
                    /(2*sigmaX2*sigmaX2)) / (math.sqrt(2*math.pi)*sigmaX2))
                    *(math.exp(-((y*math.cos(theta2) - x*math.sin(theta2))*(y*math.cos(theta2) - x*math.sin(theta2)))
                    /(2*sigmaY2*sigmaY2)) / (math.sqrt(2*math.pi)*sigmaY2));
end

local function luaKernel3(x, y)
    local sigmaX3 = 0.5;
    local sigmaY3 = 4.0;
    local theta3 = 3.14159/4; --//rotation of sigmaX axis
    return (math.exp(-((x*math.cos(theta3) +y*math.sin(theta3))*(x*math.cos(theta3) +y*math.sin(theta3)))
                    /(2*sigmaX3*sigmaX3)) / (math.sqrt(2*math.pi)*sigmaX3))
                    *(math.exp(-((y*math.cos(theta3) - x*math.sin(theta3))*(y*math.cos(theta3) - x*math.sin(theta3)))
                    /(2*sigmaY3*sigmaY3)) / (math.sqrt(2*math.pi)*sigmaY3));
end

local function luaKernel4(x, y)
    local sigmaX4 = 0.5;
    local sigmaY4 = 4.0;
    local theta4 = 3.14159/2; --//rotation of sigmaX axis
    return (math.exp(-((x*math.cos(theta4) +y*math.sin(theta4))*(x*math.cos(theta4) +y*math.sin(theta4)))
                    /(2*sigmaX4*sigmaX4)) / (math.sqrt(2*math.pi)*sigmaX4))
                    *(math.exp(-((y*math.cos(theta4) - x*math.sin(theta4))*(y*math.cos(theta4) - x*math.sin(theta4)))
                    /(2*sigmaY4*sigmaY4)) / (math.sqrt(2*math.pi)*sigmaY4));
end


local function luaKernel5(x, y)
    local sigmaX5 = 4.0;
    local sigmaY5 = 4.0;
    local theta5 = 0.0; --//rotation of sigmaX axis
    return (math.exp(-((x*math.cos(theta5) +y*math.sin(theta5))*(x*math.cos(theta5) +y*math.sin(theta5)))
                    /(2*sigmaX5*sigmaX5)) / (math.sqrt(2*math.pi)*sigmaX5))
                    *(math.exp(-((y*math.cos(theta5) - x*math.sin(theta5))*(y*math.cos(theta5) - x*math.sin(theta5)))
                    /(2*sigmaY5*sigmaY5)) / (math.sqrt(2*math.pi)*sigmaY5));
end



--terra vectormask8(a : vector(float,8), b : vector(float,8))
--    return terralib.select(a ~= b, [vector(uint16,8)](0xFFFFULL),[vector(uint16,8)](0) ) 
--end

--return a vector of kernel values at locations:
--[(x[0], y, i, j), (x[1], y, i, j), ...]
terra calcKernelLocation(xVec:vector(int, 8), y:int, i:int, j:int)
	var curKernelVec1 : vector(float, 8) = [vector(float, 8)](kernel1(i, j))
	var curKernelVec2 : vector(float, 8) = [vector(float, 8)](kernel2(i, j))
	var curKernelVec3 : vector(float, 8) = [vector(float, 8)](kernel3(i, j))
	var curKernelVec4 : vector(float, 8) = [vector(float, 8)](kernel4(i, j))
	var curKernelVec5 : vector(float, 8) = [vector(float, 8)](kernel5(i, j))

	return (vecPolynomial1(xVec, y)*curKernelVec1 +
			vecPolynomial2(xVec, y)*curKernelVec2 + vecPolynomial3(xVec, y)*curKernelVec3 + 
			vecPolynomial4(xVec, y)*curKernelVec4 + vecPolynomial5(xVec, y)*curKernelVec5)

end


--emit incorrect code if the bool testingVectorizationLR == true to test
--automatic vectorization on lightroast
local function blurImageLinCombo(method)
	local terra terraBlur()
		var width : int = 2048
		var height : int = 1489

		var input : &float
		input = [&float](C.malloc(sizeof(float)*width*height))

		var output : &float
		output = [&float](C.malloc(sizeof(float)*width*height))

		for i = 0, width*height do
			input[i] = [float](i + 1)
			output[i] = 0.0f
		end

		var t1 = C.clock()

		for y = boundingBox, height-boundingBox do
			for x = boundingBox, width-boundingBox do
				var curOut : float = 0
				var curNorm : float = 0
				var curKernelVal : float

				escape 
					if method == "testingVectorizationLR" then
						emit quote
							var prevKernelVal : float = 0.0f
							escape
								for j = -boundingBox, boundingBox do
				        			for i = -boundingBox, boundingBox do
				        				emit quote
							        		curKernelVal = prevKernelVal + polynomial1(x, y)*kernel1(i, j) +
								                polynomial2(x, y)*kernel2(i, j) + polynomial3(x, y)*kernel3(i, j) + 
								                polynomial4(x, y)*kernel4(i, j) + polynomial5(x, y)*kernel5(i, j);
								            curOut = curOut + input[(x + i) + width*(y + j)]*curKernelVal ; 
								            curNorm = curNorm + curKernelVal;
								            prevKernelVal = curKernelVal
							        	end
				        			end
				        		end
							end
						end
					elseif method == "testingRefactor" then
						local k1Out = symbol()
						local k2Out = symbol()
						local k3Out = symbol()
						local k4Out = symbol()
						local k5Out = symbol()
						local k1Sum = 0
						local k2Sum = 0
						local k3Sum = 0
						local k4Sum = 0
						local k5Sum = 0
						for j = -boundingBox, boundingBox do
					        for i = -boundingBox, boundingBox do
								k1Sum = k1Sum + luaKernel1(i, j)
								k2Sum = k2Sum + luaKernel2(i, j)
								k3Sum = k3Sum + luaKernel3(i, j)
								k4Sum = k4Sum + luaKernel4(i, j)
								k5Sum = k5Sum + luaKernel5(i, j)
					        end
					    end
						emit quote
							var [k1Out] = 0.0f
							var [k2Out] = 0.0f
							var [k3Out] = 0.0f
							var [k4Out] = 0.0f
							var [k5Out] = 0.0f
						end	
					    for j = -boundingBox, boundingBox do
					        for i = -boundingBox, boundingBox do
					        	emit quote
					        		k1Out = k1Out + input[(x + i) + width*(y + j)]*[luaKernel1(i, j)]
					        		k2Out = k2Out + input[(x + i) + width*(y + j)]*[luaKernel2(i, j)]
					        		k3Out = k3Out + input[(x + i) + width*(y + j)]*[luaKernel3(i, j)]
					        		k4Out = k4Out + input[(x + i) + width*(y + j)]*[luaKernel4(i, j)]
					        		k5Out = k5Out + input[(x + i) + width*(y + j)]*[luaKernel5(i, j)]
						        end
					        end
					    end
					    emit quote
						 	curOut = k1Out * polynomial1(x, y) +
						 				k2Out * polynomial2(x, y) + 
						 				k3Out * polynomial3(x, y) + 
						 				k4Out * polynomial4(x, y) + 
						 				k5Out * polynomial5(x, y); 
							curNorm = k1Sum * polynomial1(x, y) +
						 				k2Sum * polynomial2(x, y) + 
						 				k3Sum * polynomial3(x, y) + 
						 				k4Sum * polynomial4(x, y) + 
						 				k5Sum * polynomial5(x, y);
						end
					else
					    for j = -boundingBox, boundingBox do
					        for i = -boundingBox, boundingBox do
					        	emit quote
					        		curKernelVal = polynomial1(x, y)*kernel1(i, j) +
						                polynomial2(x, y)*kernel2(i, j) + polynomial3(x, y)*kernel3(i, j) + 
						                polynomial4(x, y)*kernel4(i, j) + polynomial5(x, y)*kernel5(i, j);
						            curOut = curOut + input[(x + i) + width*(y + j)]*curKernelVal ; 
						            curNorm = curNorm + curKernelVal;
						        end
					        end
					    end
					end
				end
			    output[x + width*y] = curOut/curNorm
			end
		end

		var fp : &C.FILE
		fp=C.fopen("./terraBlur.txt", "w");
		for i=boundingBox,height - boundingBox do
			for j=boundingBox,width - boundingBox do
				C.fprintf(fp, "%f\t", output[i*width + j]);		
			end
			C.fprintf(fp, "\n");		
		end
		C.fclose(fp);

		var t2 = C.clock()
		C.printf("\n\nlin combo %dx%d blur Terra, method = %s:\n", boundingBox, boundingBox, method)
		C.printf("output[boundingBox*width + boundingBox] = %f, computation took: %f ms\n", output[boundingBox*width + boundingBox],  (t2-t1)/1000.0f)
		C.printf("C.CLOCKS_PER_SEC = %d\n", C.CLOCKS_PER_SEC)

	    C.free(input)
		C.free(output)
	end
--	terraBlur:disas()
	terraBlur()
end
--emit incorrect code if testingVectorizationLR == true to test
--automatic vectorization on lightroast
--local terra blurImageLinCombo(testingVectorizationLR:bool)
--	[genBlurImageLinCombo(testingVectorizationLR)]
--end


local terra blurImageLinComboVectorize2()
--	var boundingBox : int = 2
	var width : int = 2048
	var height : int = 1489

	var input : &float
	input = [&float](C.malloc(sizeof(float)*width*height))

	var output : &float
	output = [&float](C.malloc(sizeof(float)*width*height))

	for i = 0, width*height do
		input[i] = [float](i + 1)
		output[i] = 0.0f
	end

	var t1 = C.clock()

	for y = boundingBox, height-boundingBox do
		for x = boundingBox, (width-boundingBox), 2 do
			--var curOut : float = 0
			--var curNorm : float = 0
			var curKernelVal : vector(float,2)
			var curKernelValTemp : float[2]
			var curOut : vector(float,2) = vector(0.f,0.f)
			var curNorm : vector(float,2) = vector(0.f,0.f)
--			var curOutVec = @[&vector(float,4)](&output[y*width + x])
--			var curOutVec : vector(float,4)
			escape 
			    for j = -boundingBox, boundingBox do
			        for i = -boundingBox, boundingBox do
			        	emit quote

--correct original
			        		curKernelValTemp[0] = polynomial1(x, y)*kernel1(i, j) +
				                polynomial2(x, y)*kernel2(i, j) + polynomial3(x, y)*kernel3(i, j) + 
				                polynomial4(x, y)*kernel4(i, j) + polynomial5(x, y)*kernel5(i, j);
			        		curKernelValTemp[1] = polynomial1(x+1, y)*kernel1(i, j) +
				                polynomial2(x+1, y)*kernel2(i, j) + polynomial3(x+1, y)*kernel3(i, j) + 
				                polynomial4(x+1, y)*kernel4(i, j) + polynomial5(x+1, y)*kernel5(i, j);

--testing destroying vectorization with incorrect code below ************************
--(it works, performance is 4x worse than above)
--				        	curKernelValTemp[0] = polynomial1(x, y)*kernel1(i, j) +
--				                polynomial2(x, y)*kernel2(i, j) + polynomial3(x, y)*kernel3(i, j) + 
--				                polynomial4(x, y)*kernel4(i, j) + polynomial5(x, y)*kernel5(i, j);
--			        		curKernelValTemp[1] = curKernelValTemp[0] + polynomial1(x+1, y)*kernel1(i, j) +
--				                polynomial2(x+1, y)*kernel2(i, j) + polynomial3(x+1, y)*kernel3(i, j) + 
--				                polynomial4(x+1, y)*kernel4(i, j) + polynomial5(x+1, y)*kernel5(i, j);
--			        		curKernelValTemp[2] = curKernelValTemp[1] + polynomial1(x+2, y)*kernel1(i, j) +
--				                polynomial2(x+2, y)*kernel2(i, j) + polynomial3(x+2, y)*kernel3(i, j) + 
--				                polynomial4(x+2, y)*kernel4(i, j) + polynomial5(x+2, y)*kernel5(i, j);
--			        		curKernelValTemp[3] = curKernelValTemp[2] + polynomial1(x+3, y)*kernel1(i, j) +
--				                polynomial2(x+3, y)*kernel2(i, j) + polynomial3(x+3, y)*kernel3(i, j) + 
--				                polynomial4(x+3, y)*kernel4(i, j) + polynomial5(x+3, y)*kernel5(i, j);
--testing destroying vectorization with incorrect code above ************************
	            

				            curKernelVal = @[&vector(float,2)](&curKernelValTemp[0])

--							curKernelVal = [vector(float, 4)](polynomial1(x, y)*kernel1(i, j) +
--				                polynomial2(x, y)*kernel2(i, j) + polynomial3(x, y)*kernel3(i, j) + 
--				                polynomial4(x, y)*kernel4(i, j) + polynomial5(x, y)*kernel5(i, j));

				            var curInVec = @[&vector(float,2)](&input[(y+j)*width + x+i])

				            curOut = curOut + curInVec*curKernelVal ; 
				            curNorm = curNorm + curKernelVal;
				        end
			        end
			    end
			end
		    curOut = curOut/curNorm
		    output[y*width + x] = curOut[0]
		    output[y*width + x+1] = curOut[1]
--		    output[x + width*y] = curOut/curNorm
		end
	end

	var t2 = C.clock()
	C.printf("\n\nVectorized2 lin combo %dx%d blur Terra:\n", boundingBox, boundingBox)
	C.printf("output[boundingBox*width + boundingBox] = %f, computation took: %f ms\n", output[boundingBox*width + boundingBox],  (t2-t1)/1000.0f)
	C.printf("C.CLOCKS_PER_SEC = %d\n", C.CLOCKS_PER_SEC)

	for i=boundingBox,boundingBox+10 do
		for j=boundingBox,boundingBox+10 do
			C.printf("%f\t", output[i*width + j])
		end
		C.printf("\n")
	end

	C.printf(":)\n")
    C.free(input)
	C.free(output)
end


local terra blurImageLinComboVectorize4()
--	var boundingBox : int = 2
	var width : int = 2048
	var height : int = 1489

	var input : &float
	input = [&float](C.malloc(sizeof(float)*width*height))

	var output : &float
	output = [&float](C.malloc(sizeof(float)*width*height))

	for i = 0, width*height do
		input[i] = [float](i + 1)
		output[i] = 0.0f
	end

	var t1 = C.clock()

	for y = boundingBox, height-boundingBox do
		for x = boundingBox, (width-boundingBox), 4 do
			--var curOut : float = 0
			--var curNorm : float = 0
			var curKernelVal : vector(float,4)
			var curKernelValTemp : float[4]
			var curOut : vector(float,4) = vector(0.f,0.f,0.f,0.f)
			var curNorm : vector(float,4) = vector(0.f,0.f,0.f,0.f)
--			var curOutVec = @[&vector(float,4)](&output[y*width + x])
--			var curOutVec : vector(float,4)
			escape 
			    for j = -boundingBox, boundingBox do
			        for i = -boundingBox, boundingBox do
			        	emit quote

--correct original
			        		curKernelValTemp[0] = polynomial1(x, y)*kernel1(i, j) +
				                polynomial2(x, y)*kernel2(i, j) + polynomial3(x, y)*kernel3(i, j) + 
				                polynomial4(x, y)*kernel4(i, j) + polynomial5(x, y)*kernel5(i, j);
			        		curKernelValTemp[1] = polynomial1(x+1, y)*kernel1(i, j) +
				                polynomial2(x+1, y)*kernel2(i, j) + polynomial3(x+1, y)*kernel3(i, j) + 
				                polynomial4(x+1, y)*kernel4(i, j) + polynomial5(x+1, y)*kernel5(i, j);
			        		curKernelValTemp[2] = polynomial1(x+2, y)*kernel1(i, j) +
				                polynomial2(x+2, y)*kernel2(i, j) + polynomial3(x+2, y)*kernel3(i, j) + 
				                polynomial4(x+2, y)*kernel4(i, j) + polynomial5(x+2, y)*kernel5(i, j);
			        		curKernelValTemp[3] = polynomial1(x+3, y)*kernel1(i, j) +
				                polynomial2(x+3, y)*kernel2(i, j) + polynomial3(x+3, y)*kernel3(i, j) + 
				                polynomial4(x+3, y)*kernel4(i, j) + polynomial5(x+3, y)*kernel5(i, j);

--testing destroying vectorization with incorrect code below ************************
--(it works, performance is 4x worse than above)
--				        	curKernelValTemp[0] = polynomial1(x, y)*kernel1(i, j) +
--				                polynomial2(x, y)*kernel2(i, j) + polynomial3(x, y)*kernel3(i, j) + 
--				                polynomial4(x, y)*kernel4(i, j) + polynomial5(x, y)*kernel5(i, j);
--			        		curKernelValTemp[1] = curKernelValTemp[0] + polynomial1(x+1, y)*kernel1(i, j) +
--				                polynomial2(x+1, y)*kernel2(i, j) + polynomial3(x+1, y)*kernel3(i, j) + 
--				                polynomial4(x+1, y)*kernel4(i, j) + polynomial5(x+1, y)*kernel5(i, j);
--			        		curKernelValTemp[2] = curKernelValTemp[1] + polynomial1(x+2, y)*kernel1(i, j) +
--				                polynomial2(x+2, y)*kernel2(i, j) + polynomial3(x+2, y)*kernel3(i, j) + 
--				                polynomial4(x+2, y)*kernel4(i, j) + polynomial5(x+2, y)*kernel5(i, j);
--			        		curKernelValTemp[3] = curKernelValTemp[2] + polynomial1(x+3, y)*kernel1(i, j) +
--				                polynomial2(x+3, y)*kernel2(i, j) + polynomial3(x+3, y)*kernel3(i, j) + 
--				                polynomial4(x+3, y)*kernel4(i, j) + polynomial5(x+3, y)*kernel5(i, j);
--testing destroying vectorization with incorrect code above ************************
	            

				            curKernelVal = @[&vector(float,4)](&curKernelValTemp[0])

--							curKernelVal = [vector(float, 4)](polynomial1(x, y)*kernel1(i, j) +
--				                polynomial2(x, y)*kernel2(i, j) + polynomial3(x, y)*kernel3(i, j) + 
--				                polynomial4(x, y)*kernel4(i, j) + polynomial5(x, y)*kernel5(i, j));

				            var curInVec = @[&vector(float,4)](&input[(y+j)*width + x+i])

				            curOut = curOut + curInVec*curKernelVal ; 
				            curNorm = curNorm + curKernelVal;
				        end
			        end
			    end
			end
		    curOut = curOut/curNorm
		    output[y*width + x] = curOut[0]
		    output[y*width + x+1] = curOut[1]
		    output[y*width + x+2] = curOut[2]
		    output[y*width + x+3] = curOut[3]
--		    output[x + width*y] = curOut/curNorm
		end
	end

	var t2 = C.clock()
	C.printf("\n\nVectorized84 lin combo %dx%d blur Terra:\n", boundingBox, boundingBox)
	C.printf("output[boundingBox*width + boundingBox] = %f, computation took: %f ms\n", output[boundingBox*width + boundingBox],  (t2-t1)/1000.0f)
	C.printf("C.CLOCKS_PER_SEC = %d\n", C.CLOCKS_PER_SEC)

	for i=boundingBox,boundingBox+10 do
		for j=boundingBox,boundingBox+10 do
			C.printf("%f\t", output[i*width + j])
		end
		C.printf("\n")
	end

	C.printf(":)\n")
    C.free(input)
	C.free(output)
end


local terra blurImageLinComboVectorize8()
--	var boundingBox : int = 2
	var width : int = 2048
	var height : int = 1489

	var input : &float
	input = [&float](C.malloc(sizeof(float)*width*height))

	var output : &float
	output = [&float](C.malloc(sizeof(float)*width*height))

	for i = 0, width*height do
		input[i] = [float](i + 1)
		output[i] = 0.0f
	end

	var t1 = C.clock()

	for y = boundingBox, height-boundingBox do
		for x = boundingBox, (width-boundingBox), 8 do
			--var curOut : float = 0
			--var curNorm : float = 0
			var curKernelVal : vector(float,8)
			var curKernelValTemp : float[8]
			var curOut : vector(float,8) = vector(0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f)
			var curNorm : vector(float,8) = vector(0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f)
--			var curOutVec = @[&vector(float,4)](&output[y*width + x])
--			var curOutVec : vector(float,4)
			escape 
			    for j = -boundingBox, boundingBox do
			        for i = -boundingBox, boundingBox do
			        	emit quote

--correct original
			        		curKernelValTemp[0] = polynomial1(x, y)*kernel1(i, j) +
				                polynomial2(x, y)*kernel2(i, j) + polynomial3(x, y)*kernel3(i, j) + 
				                polynomial4(x, y)*kernel4(i, j) + polynomial5(x, y)*kernel5(i, j);
			        		curKernelValTemp[1] = polynomial1(x+1, y)*kernel1(i, j) +
				                polynomial2(x+1, y)*kernel2(i, j) + polynomial3(x+1, y)*kernel3(i, j) + 
				                polynomial4(x+1, y)*kernel4(i, j) + polynomial5(x+1, y)*kernel5(i, j);
			        		curKernelValTemp[2] = polynomial1(x+2, y)*kernel1(i, j) +
				                polynomial2(x+2, y)*kernel2(i, j) + polynomial3(x+2, y)*kernel3(i, j) + 
				                polynomial4(x+2, y)*kernel4(i, j) + polynomial5(x+2, y)*kernel5(i, j);
			        		curKernelValTemp[3] = polynomial1(x+3, y)*kernel1(i, j) +
				                polynomial2(x+3, y)*kernel2(i, j) + polynomial3(x+3, y)*kernel3(i, j) + 
				                polynomial4(x+3, y)*kernel4(i, j) + polynomial5(x+3, y)*kernel5(i, j);
			        		curKernelValTemp[4] = polynomial1(x+4, y)*kernel1(i, j) +
				                polynomial2(x+4, y)*kernel2(i, j) + polynomial3(x+4, y)*kernel3(i, j) + 
				                polynomial4(x+4, y)*kernel4(i, j) + polynomial5(x+4, y)*kernel5(i, j);
			        		curKernelValTemp[5] = polynomial1(x+5, y)*kernel1(i, j) +
				                polynomial2(x+5, y)*kernel2(i, j) + polynomial3(x+5, y)*kernel3(i, j) + 
				                polynomial4(x+5, y)*kernel4(i, j) + polynomial5(x+5, y)*kernel5(i, j);
			        		curKernelValTemp[6] = polynomial1(x+6, y)*kernel1(i, j) +
				                polynomial2(x+6, y)*kernel2(i, j) + polynomial3(x+6, y)*kernel3(i, j) + 
				                polynomial4(x+6, y)*kernel4(i, j) + polynomial5(x+6, y)*kernel5(i, j);
			        		curKernelValTemp[7] = polynomial1(x+7, y)*kernel1(i, j) +
				                polynomial2(x+7, y)*kernel2(i, j) + polynomial3(x+7, y)*kernel3(i, j) + 
				                polynomial4(x+7, y)*kernel4(i, j) + polynomial5(x+7, y)*kernel5(i, j);

--testing destroying vectorization with incorrect code below ************************
--(it works, performance is 4x worse than above)
--				        	curKernelValTemp[0] = polynomial1(x, y)*kernel1(i, j) +
--				                polynomial2(x, y)*kernel2(i, j) + polynomial3(x, y)*kernel3(i, j) + 
--				                polynomial4(x, y)*kernel4(i, j) + polynomial5(x, y)*kernel5(i, j);
--			        		curKernelValTemp[1] = curKernelValTemp[0] + polynomial1(x+1, y)*kernel1(i, j) +
--				                polynomial2(x+1, y)*kernel2(i, j) + polynomial3(x+1, y)*kernel3(i, j) + 
--				                polynomial4(x+1, y)*kernel4(i, j) + polynomial5(x+1, y)*kernel5(i, j);
--			        		curKernelValTemp[2] = curKernelValTemp[1] + polynomial1(x+2, y)*kernel1(i, j) +
--				                polynomial2(x+2, y)*kernel2(i, j) + polynomial3(x+2, y)*kernel3(i, j) + 
--				                polynomial4(x+2, y)*kernel4(i, j) + polynomial5(x+2, y)*kernel5(i, j);
--			        		curKernelValTemp[3] = curKernelValTemp[2] + polynomial1(x+3, y)*kernel1(i, j) +
--				                polynomial2(x+3, y)*kernel2(i, j) + polynomial3(x+3, y)*kernel3(i, j) + 
--				                polynomial4(x+3, y)*kernel4(i, j) + polynomial5(x+3, y)*kernel5(i, j);
--testing destroying vectorization with incorrect code above ************************
	            

				            curKernelVal = @[&vector(float,8)](&curKernelValTemp[0])

--							curKernelVal = [vector(float, 4)](polynomial1(x, y)*kernel1(i, j) +
--				                polynomial2(x, y)*kernel2(i, j) + polynomial3(x, y)*kernel3(i, j) + 
--				                polynomial4(x, y)*kernel4(i, j) + polynomial5(x, y)*kernel5(i, j));

				            var curInVec = @[&vector(float,8)](&input[(y+j)*width + x+i])

				            curOut = curOut + curInVec*curKernelVal ; 
				            curNorm = curNorm + curKernelVal;
				        end
			        end
			    end
			end
		    curOut = curOut/curNorm
		    output[y*width + x] = curOut[0]
		    output[y*width + x+1] = curOut[1]
		    output[y*width + x+2] = curOut[2]
		    output[y*width + x+3] = curOut[3]
		    output[y*width + x+1] = curOut[4]
		    output[y*width + x+1] = curOut[5]
		    output[y*width + x+2] = curOut[6]
		    output[y*width + x+3] = curOut[7]
--		    output[x + width*y] = curOut/curNorm
		end
	end

	var t2 = C.clock()
	C.printf("\n\nVectorized8 lin combo %dx%d blur Terra:\n", boundingBox, boundingBox)
	C.printf("output[boundingBox*width + boundingBox] = %f, computation took: %f ms\n", output[boundingBox*width + boundingBox],  (t2-t1)/1000.0f)
	C.printf("C.CLOCKS_PER_SEC = %d\n", C.CLOCKS_PER_SEC)

	for i=boundingBox,boundingBox+10 do
		for j=boundingBox,boundingBox+10 do
			C.printf("%f\t", output[i*width + j])
		end
		C.printf("\n")
	end

	C.printf(":)\n")
    C.free(input)
	C.free(output)
end


local terra blurImageLinComboVectorize8Wide1()
--	var boundingBox : int = 2
	var width : int = 2048
	var height : int = 1489

	var input : &float
	input = [&float](C.malloc(sizeof(float)*width*height))

	var output : &float
	output = [&float](C.malloc(sizeof(float)*width*height))

	for i = 0, width*height do
		input[i] = [float](i + 1)
		output[i] = 0.0f
	end

	var t1 = C.clock()

	for y = boundingBox, height-boundingBox do
		for x = boundingBox, (width-boundingBox), 8 do
			--var curOut : float = 0
			--var curNorm : float = 0
			var curKernelVal : vector(float,8)
			var curOut : vector(float,8) = vector(0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f)
			var curNorm : vector(float,8) = vector(0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f)
--			var curOutVec = @[&vector(float,4)](&output[y*width + x])
--			var curOutVec : vector(float,4)
			escape 
			    for j = -boundingBox, boundingBox do
			        for i = -boundingBox, boundingBox do
			        	emit quote
			        		var xVec : vector(int, 8) = vector(x, x+1, x+2, x+3,
			        											x+4, x+5, x+6, x+7)
				            curKernelVal = calcKernelLocation(xVec, y, i, j)

				            var curInVec = @[&vector(float,8)](&input[(y+j)*width + x+i])

				            curOut = curOut + curInVec*curKernelVal; 
				            curNorm = curNorm + curKernelVal;
				        end
			        end
			    end
			end
		    curOut = curOut/curNorm
		    output[y*width + x] = curOut[0]
		    output[y*width + x+1] = curOut[1]
		    output[y*width + x+2] = curOut[2]
		    output[y*width + x+3] = curOut[3]
		    output[y*width + x+4] = curOut[4]
		    output[y*width + x+5] = curOut[5]
		    output[y*width + x+6] = curOut[6]
		    output[y*width + x+7] = curOut[7]
--		    output[x + width*y] = curOut/curNorm
		end
	end

	var t2 = C.clock()
	C.printf("\n\nVectorized 8wide lin combo %dx%d blur Terra:\n", boundingBox, boundingBox)
	C.printf("output[boundingBox*width + boundingBox] = %f, computation took: %f ms\n", output[boundingBox*width + boundingBox],  (t2-t1)/1000.0f)
	C.printf("C.CLOCKS_PER_SEC = %d\n", C.CLOCKS_PER_SEC)

	for i=boundingBox,boundingBox+10 do
		for j=boundingBox,boundingBox+10 do
			C.printf("%f\t", output[i*width + j])
		end
		C.printf("\n")
	end

	C.printf(":)\n")
    C.free(input)
	C.free(output)
end



local terra blurImageLinComboVectorize16Wide()
--	var boundingBox : int = 2
	var width : int = 2048
	var height : int = 1489

	var input : &float
	input = [&float](C.malloc(sizeof(float)*width*height))

	var output : &float
	output = [&float](C.malloc(sizeof(float)*width*height))

	for i = 0, width*height do
		input[i] = [float](i + 1)
		output[i] = 0.0f
	end

	var t1 = C.clock()

	for y = boundingBox, height-boundingBox do
		for x = boundingBox, (width-boundingBox), 16 do
			--var curOut : float = 0
			--var curNorm : float = 0
			var curKernelVal : vector(float,16)
			var curKernelValTemp : float[16]
			var curOut : vector(float,16) = vector(0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f)
			var curNorm : vector(float,16) = vector(0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f)
--			var curOutVec = @[&vector(float,4)](&output[y*width + x])
--			var curOutVec : vector(float,4)
			escape 
			    for j = -boundingBox, boundingBox do
			        for i = -boundingBox, boundingBox do
			        	emit quote
			        		curKernelValTemp[0] = polynomial1(x, y)*kernel1(i, j) +
				                polynomial2(x, y)*kernel2(i, j) + polynomial3(x, y)*kernel3(i, j) + 
				                polynomial4(x, y)*kernel4(i, j) + polynomial5(x, y)*kernel5(i, j);
			        		curKernelValTemp[1] = polynomial1(x+1, y)*kernel1(i, j) +
				                polynomial2(x+1, y)*kernel2(i, j) + polynomial3(x+1, y)*kernel3(i, j) + 
				                polynomial4(x+1, y)*kernel4(i, j) + polynomial5(x+1, y)*kernel5(i, j);
			        		curKernelValTemp[2] = polynomial1(x+2, y)*kernel1(i, j) +
				                polynomial2(x+2, y)*kernel2(i, j) + polynomial3(x+2, y)*kernel3(i, j) + 
				                polynomial4(x+2, y)*kernel4(i, j) + polynomial5(x+2, y)*kernel5(i, j);
			        		curKernelValTemp[3] = polynomial1(x+3, y)*kernel1(i, j) +
				                polynomial2(x+3, y)*kernel2(i, j) + polynomial3(x+3, y)*kernel3(i, j) + 
				                polynomial4(x+3, y)*kernel4(i, j) + polynomial5(x+3, y)*kernel5(i, j);
			        		curKernelValTemp[4] = polynomial1(x+4, y)*kernel1(i, j) +
				                polynomial2(x+4, y)*kernel2(i, j) + polynomial3(x+4, y)*kernel3(i, j) + 
				                polynomial4(x+4, y)*kernel4(i, j) + polynomial5(x+4, y)*kernel5(i, j);		
			        		curKernelValTemp[5] = polynomial1(x+5, y)*kernel1(i, j) +
				                polynomial2(x+5, y)*kernel2(i, j) + polynomial3(x+5, y)*kernel3(i, j) + 
				                polynomial4(x+5, y)*kernel4(i, j) + polynomial5(x+5, y)*kernel5(i, j);
			        		curKernelValTemp[6] = polynomial1(x+6, y)*kernel1(i, j) +
				                polynomial2(x+6, y)*kernel2(i, j) + polynomial3(x+6, y)*kernel3(i, j) + 
				                polynomial4(x+6, y)*kernel4(i, j) + polynomial5(x+6, y)*kernel5(i, j);
			        		curKernelValTemp[7] = polynomial1(x+7, y)*kernel1(i, j) +
				                polynomial2(x+7, y)*kernel2(i, j) + polynomial3(x+7, y)*kernel3(i, j) + 
				                polynomial4(x+7, y)*kernel4(i, j) + polynomial5(x+7, y)*kernel5(i, j);				            

			        		curKernelValTemp[8] = polynomial1(x+8, y)*kernel1(i, j) +
				                polynomial2(x+8, y)*kernel2(i, j) + polynomial3(x+8, y)*kernel3(i, j) + 
				                polynomial4(x+8, y)*kernel4(i, j) + polynomial5(x+8, y)*kernel5(i, j);	
			        		curKernelValTemp[9] = polynomial1(x+9, y)*kernel1(i, j) +
				                polynomial2(x+9, y)*kernel2(i, j) + polynomial3(x+9, y)*kernel3(i, j) + 
				                polynomial4(x+9, y)*kernel4(i, j) + polynomial5(x+9, y)*kernel5(i, j);
			        		curKernelValTemp[10] = polynomial1(x+10, y)*kernel1(i, j) +
				                polynomial2(x+10, y)*kernel2(i, j) + polynomial3(x+10, y)*kernel3(i, j) + 
				                polynomial4(x+10, y)*kernel4(i, j) + polynomial5(x+10, y)*kernel5(i, j);
			        		curKernelValTemp[11] = polynomial1(x+11, y)*kernel1(i, j) +
				                polynomial2(x+11, y)*kernel2(i, j) + polynomial3(x+11, y)*kernel3(i, j) + 
				                polynomial4(x+11, y)*kernel4(i, j) + polynomial5(x+11, y)*kernel5(i, j);
			        		curKernelValTemp[12] = polynomial1(x+12, y)*kernel1(i, j) +
				                polynomial2(x+12, y)*kernel2(i, j) + polynomial3(x+12, y)*kernel3(i, j) + 
				                polynomial4(x+12, y)*kernel4(i, j) + polynomial5(x+12, y)*kernel5(i, j);		
			        		curKernelValTemp[13] = polynomial1(x+13, y)*kernel1(i, j) +
				                polynomial2(x+13, y)*kernel2(i, j) + polynomial3(x+13, y)*kernel3(i, j) + 
				                polynomial4(x+13, y)*kernel4(i, j) + polynomial5(x+13, y)*kernel5(i, j);
			        		curKernelValTemp[14] = polynomial1(x+14, y)*kernel1(i, j) +
				                polynomial2(x+14, y)*kernel2(i, j) + polynomial3(x+14, y)*kernel3(i, j) + 
				                polynomial4(x+14, y)*kernel4(i, j) + polynomial5(x+14, y)*kernel5(i, j);
			        		curKernelValTemp[15] = polynomial1(x+15, y)*kernel1(i, j) +
				                polynomial2(x+15, y)*kernel2(i, j) + polynomial3(x+15, y)*kernel3(i, j) + 
				                polynomial4(x+15, y)*kernel4(i, j) + polynomial5(x+15, y)*kernel5(i, j);

				            curKernelVal = @[&vector(float,16)](&curKernelValTemp[0])

				            var curInVec = @[&vector(float,16)](&input[(y+j)*width + x+i])

				            curOut = curOut + curInVec*curKernelVal; 
				            curNorm = curNorm + curKernelVal;
				        end
			        end
			    end
			end
		    curOut = curOut/curNorm
		    output[y*width + x] = curOut[0]
		    output[y*width + x+1] = curOut[1]
		    output[y*width + x+2] = curOut[2]
		    output[y*width + x+3] = curOut[3]
		    output[y*width + x+4] = curOut[4]
		    output[y*width + x+5] = curOut[5]
		    output[y*width + x+6] = curOut[6]
		    output[y*width + x+7] = curOut[7]

		    output[y*width + x+8] = curOut[8]
		    output[y*width + x+9] = curOut[9]
		    output[y*width + x+10] = curOut[10]
		    output[y*width + x+11] = curOut[11]
		    output[y*width + x+12] = curOut[12]
		    output[y*width + x+13] = curOut[13]
		    output[y*width + x+14] = curOut[14]
		    output[y*width + x+15] = curOut[15]
--		    output[x + width*y] = curOut/curNorm
		end
	end

	var t2 = C.clock()
	C.printf("\n\nVectorized 16wide lin combo %dx%d blur Terra:\n", boundingBox, boundingBox)
	C.printf("output[boundingBox*width + boundingBox] = %f, computation took: %f ms\n", output[boundingBox*width + boundingBox],  (t2-t1)/1000.0f)
	C.printf("C.CLOCKS_PER_SEC = %d\n", C.CLOCKS_PER_SEC)

	for i=boundingBox,boundingBox+10 do
		for j=boundingBox,boundingBox+10 do
			C.printf("%f\t", output[i*width + j])
		end
		C.printf("\n")
	end

	C.printf(":)\n")
    C.free(input)
	C.free(output)
end


--deal with image, mask, and variance planes
local terra blurMaskedImageLinCombo()
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

--	var xSplit : int = 32
--	for xOuter = boundingBox, (width-boundingBox-xSplit), xSplit do
		for y = boundingBox, height-boundingBox do
			for x = boundingBox, (width-boundingBox) do
--			for x = xOuter, xOuter+xSplit, 4 do
				var curKernelVal : float

				var curImOut : float = 0.0f
				var curVarOut : float = 0.0f
				var curMaskOut : uint16 = 0

				var curNorm : float = 0.0f

				escape 
				    for j = -boundingBox, boundingBox do
				        for i = -boundingBox, boundingBox do
				        	emit quote
				        		curKernelVal = polynomial1(x, y)*kernel1(i, j) +
					                polynomial2(x, y)*kernel2(i, j) + polynomial3(x, y)*kernel3(i, j) + 
					                polynomial4(x, y)*kernel4(i, j) + polynomial5(x, y)*kernel5(i, j);

					            var curImInVec = inputIm[(y+j)*width + x+i]
					            var curVarInVec = inputVar[(y+j)*width + x+i]
					            var curMaskInVec = inputMask[(y+j)*width + x+i]


					            curImOut = curImOut + curImInVec*curKernelVal; 
					            curVarOut = curVarOut + curVarInVec*curKernelVal*curKernelVal;

					            if curKernelVal ~= 0.0f then
					            	curMaskOut = curMaskOut or curMaskInVec 
					            end

					            curNorm = curNorm + curKernelVal;
					        end
				        end
				    end
				end
			    curImOut = curImOut/curNorm
			    outputIm[y*width + x] = curImOut

			    curVarOut = curVarOut/(curNorm*curNorm)
			    outputVar[y*width + x] = curVarOut

			    outputMask[y*width + x] = curMaskOut

			end
		end
	--end

	var t2 = C.clock()
	C.printf("\n\nVectorized masked image lin combo %dx%d blur Terra:\n", boundingBox, boundingBox)
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



terra vectormask(a : vector(float,4), b : vector(float,4))
    return terralib.select(a ~= b, [vector(uint16,4)](0xFFFFULL),[vector(uint16,4)](0) ) 
end
--deal with image, mask, and variance planes
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

	for i = 0, width*height do
		inputIm[i] = [float](i + 1)
		inputVar[i] = [float](i + 2)
		inputMask[i] = [uint16](i%width + i/width)


		outputIm[i] = 0.0f
		outputVar[i] = 0.0f
		outputMask[i] = 0
	end

	var t1 = C.clock()

	var xSplit : int = 32
	for xOuter = boundingBox, (width-boundingBox-xSplit), xSplit do
		for y = boundingBox, height-boundingBox do
--			for x = boundingBox, (width-boundingBox), 4 do
			for x = xOuter, xOuter+xSplit, 4 do
				--var curOut : float = 0
				--var curNorm : float = 0
				var curKernelVal : vector(float,4)
				var curKernelValTemp : float[4]
				var curImOut : vector(float,4) = vector(0.f,0.f,0.f,0.f)
				var curVarOut : vector(float,4) = vector(0.f,0.f,0.f,0.f)

	--			var curMaskOut : vector(uint16,4) = vector(0,0,0,0)

				var curMaskOut : vector(uint16, 4) = vector(0, 0, 0, 0)
				var zeroVec : vector(float, 4) = vector(0.f, 0.f, 0.f, 0.f)

				var curNorm : vector(float,4) = vector(0.f,0.f,0.f,0.f)
	--			var curImOutVec = @[&vector(float,4)](&outputIm[y*width + x])


				escape 
				    for j = -boundingBox, boundingBox do
				        for i = -boundingBox, boundingBox do
				        	emit quote
				        		curKernelValTemp[0] = polynomial1(x, y)*kernel1(i, j) +
					                polynomial2(x, y)*kernel2(i, j) + polynomial3(x, y)*kernel3(i, j) + 
					                polynomial4(x, y)*kernel4(i, j) + polynomial5(x, y)*kernel5(i, j);
				        		curKernelValTemp[1] = polynomial1(x+1, y)*kernel1(i, j) +
					                polynomial2(x+1, y)*kernel2(i, j) + polynomial3(x+1, y)*kernel3(i, j) + 
					                polynomial4(x+1, y)*kernel4(i, j) + polynomial5(x+1, y)*kernel5(i, j);
				        		curKernelValTemp[2] = polynomial1(x+2, y)*kernel1(i, j) +
					                polynomial2(x+2, y)*kernel2(i, j) + polynomial3(x+2, y)*kernel3(i, j) + 
					                polynomial4(x+2, y)*kernel4(i, j) + polynomial5(x+2, y)*kernel5(i, j);
				        		curKernelValTemp[3] = polynomial1(x+3, y)*kernel1(i, j) +
					                polynomial2(x+3, y)*kernel2(i, j) + polynomial3(x+3, y)*kernel3(i, j) + 
					                polynomial4(x+3, y)*kernel4(i, j) + polynomial5(x+3, y)*kernel5(i, j);
					            
					            curKernelVal = @[&vector(float,4)](&curKernelValTemp[0])

					            var curImInVec = @[&vector(float,4)](&inputIm[(y+j)*width + x+i])
					            var curVarInVec = @[&vector(float,4)](&inputVar[(y+j)*width + x+i])
					            var curMaskInVec = @[&vector(uint16,4)](&inputMask[(y+j)*width + x+i])


					            curImOut = curImOut + curImInVec*curKernelVal; 
					            curVarOut = curVarOut + curVarInVec*curKernelVal*curKernelVal;

					            var bitMask : vector(uint16, 4) = vectormask(curKernelVal, zeroVec)
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
	--		    outputIm[x + width*y] = curImOut/curNorm

			    curVarOut = curVarOut/(curNorm*curNorm)
			    outputVar[y*width + x] = curVarOut[0]
			    outputVar[y*width + x+1] = curVarOut[1]
			    outputVar[y*width + x+2] = curVarOut[2]
			    outputVar[y*width + x+3] = curVarOut[3]

			    outputMask[y*width + x] = curMaskOut[0]
			    outputMask[y*width + x+1] = curMaskOut[1]
			    outputMask[y*width + x+2] = curMaskOut[2]
			    outputMask[y*width + x+3] = curMaskOut[3]

			end
		end
	end

	var t2 = C.clock()
	C.printf("\n\nVectorized masked image lin combo %dx%d blur Terra:\n", boundingBox, boundingBox)
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


terra vectormask8(a : vector(float,8), b : vector(float,8))
    return terralib.select(a ~= b, [vector(uint16,8)](0xFFFFULL),[vector(uint16,8)](0) ) 
end
local terra blurMaskedImageLinComboVectorize8()
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


			escape 
			    for j = -boundingBox, boundingBox do
			        for i = -boundingBox, boundingBox do
			        	emit quote
			        		curKernelValTemp[0] = polynomial1(x, y)*kernel1(i, j) +
				                polynomial2(x, y)*kernel2(i, j) + polynomial3(x, y)*kernel3(i, j) + 
				                polynomial4(x, y)*kernel4(i, j) + polynomial5(x, y)*kernel5(i, j);
			        		curKernelValTemp[1] = polynomial1(x+1, y)*kernel1(i, j) +
				                polynomial2(x+1, y)*kernel2(i, j) + polynomial3(x+1, y)*kernel3(i, j) + 
				                polynomial4(x+1, y)*kernel4(i, j) + polynomial5(x+1, y)*kernel5(i, j);
			        		curKernelValTemp[2] = polynomial1(x+2, y)*kernel1(i, j) +
				                polynomial2(x+2, y)*kernel2(i, j) + polynomial3(x+2, y)*kernel3(i, j) + 
				                polynomial4(x+2, y)*kernel4(i, j) + polynomial5(x+2, y)*kernel5(i, j);
			        		curKernelValTemp[3] = polynomial1(x+3, y)*kernel1(i, j) +
				                polynomial2(x+3, y)*kernel2(i, j) + polynomial3(x+3, y)*kernel3(i, j) + 
				                polynomial4(x+3, y)*kernel4(i, j) + polynomial5(x+3, y)*kernel5(i, j);
			        		curKernelValTemp[4] = polynomial1(x+4, y)*kernel1(i, j) +
				                polynomial2(x+4, y)*kernel2(i, j) + polynomial3(x+4, y)*kernel3(i, j) + 
				                polynomial4(x+4, y)*kernel4(i, j) + polynomial5(x+4, y)*kernel5(i, j);		
			        		curKernelValTemp[5] = polynomial1(x+5, y)*kernel1(i, j) +
				                polynomial2(x+5, y)*kernel2(i, j) + polynomial3(x+5, y)*kernel3(i, j) + 
				                polynomial4(x+5, y)*kernel4(i, j) + polynomial5(x+5, y)*kernel5(i, j);
			        		curKernelValTemp[6] = polynomial1(x+6, y)*kernel1(i, j) +
				                polynomial2(x+6, y)*kernel2(i, j) + polynomial3(x+6, y)*kernel3(i, j) + 
				                polynomial4(x+6, y)*kernel4(i, j) + polynomial5(x+6, y)*kernel5(i, j);
			        		curKernelValTemp[7] = polynomial1(x+7, y)*kernel1(i, j) +
				                polynomial2(x+7, y)*kernel2(i, j) + polynomial3(x+7, y)*kernel3(i, j) + 
				                polynomial4(x+7, y)*kernel4(i, j) + polynomial5(x+7, y)*kernel5(i, j);
				            
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


--deal with image, mask, and variance planes
local terra blurMaskedImageOrAllMaskLinComboVectorize()
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
		for x = boundingBox, (width-boundingBox), 4 do
			--var curOut : float = 0
			--var curNorm : float = 0
			var curKernelVal : vector(float,4)
			var curKernelValTemp : float[4]
			var curImOut : vector(float,4) = vector(0.f,0.f,0.f,0.f)
			var curVarOut : vector(float,4) = vector(0.f,0.f,0.f,0.f)

			var curMaskOut : vector(uint16,4) = vector(0,0,0,0)

			var curNorm : vector(float,4) = vector(0.f,0.f,0.f,0.f)
			var zeroVec1 : vector(float,4) = vector(0.f,0.f,0.f,0.f)
--			var curImOutVec = @[&vector(float,4)](&outputIm[y*width + x])
			escape 
			    for j = -boundingBox, boundingBox do
			        for i = -boundingBox, boundingBox do
			        	emit quote
			        		curKernelValTemp[0] = polynomial1(x, y)*kernel1(i, j) +
				                polynomial2(x, y)*kernel2(i, j) + polynomial3(x, y)*kernel3(i, j) + 
				                polynomial4(x, y)*kernel4(i, j) + polynomial5(x, y)*kernel5(i, j);
			        		curKernelValTemp[1] = polynomial1(x+1, y)*kernel1(i, j) +
				                polynomial2(x+1, y)*kernel2(i, j) + polynomial3(x+1, y)*kernel3(i, j) + 
				                polynomial4(x+1, y)*kernel4(i, j) + polynomial5(x+1, y)*kernel5(i, j);
			        		curKernelValTemp[2] = polynomial1(x+2, y)*kernel1(i, j) +
				                polynomial2(x+2, y)*kernel2(i, j) + polynomial3(x+2, y)*kernel3(i, j) + 
				                polynomial4(x+2, y)*kernel4(i, j) + polynomial5(x+2, y)*kernel5(i, j);
			        		curKernelValTemp[3] = polynomial1(x+3, y)*kernel1(i, j) +
				                polynomial2(x+3, y)*kernel2(i, j) + polynomial3(x+3, y)*kernel3(i, j) + 
				                polynomial4(x+3, y)*kernel4(i, j) + polynomial5(x+3, y)*kernel5(i, j);
				            
				            curKernelVal = @[&vector(float,4)](&curKernelValTemp[0])

				            var curImInVec = @[&vector(float,4)](&inputIm[(y+j)*width + x+i])
				            var curVarInVec = @[&vector(float,4)](&inputVar[(y+j)*width + x+i])
				            var curMaskInVec = @[&vector(uint16,4)](&inputMask[(y+j)*width + x+i])



				            curImOut = curImOut + curImInVec*curKernelVal; 
				            curVarOut = curVarOut + curVarInVec*curKernelVal*curKernelVal;
				        	curMaskOut = curMaskOut or curMaskInVec

--				            curMaskOut = curMaskOut or curMaskInVec 
--				            curMaskOut = curMaskOut or (curMaskInVec and [vector(uint16,4)]([vector(bool,4)](curKernelVal))) 
--				            curMaskOut = curMaskOut or (curMaskInVec * (curKernelVal)) 


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
--		    outputIm[x + width*y] = curImOut/curNorm

		    curVarOut = curVarOut/(curNorm*curNorm)
		    outputVar[y*width + x] = curVarOut[0]
		    outputVar[y*width + x+1] = curVarOut[1]
		    outputVar[y*width + x+2] = curVarOut[2]
		    outputVar[y*width + x+3] = curVarOut[3]

		    outputMask[y*width + x] = curMaskOut[0]
		    outputMask[y*width + x+1] = curMaskOut[1]
		    outputMask[y*width + x+2] = curMaskOut[2]
		    outputMask[y*width + x+3] = curMaskOut[3]

		end
	end

	var t2 = C.clock()
	C.printf("\n\nVectorized masked image lin combo, ORing all mask pixels %dx%d blur Terra:\n", boundingBox, boundingBox)
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


--the first term is the name in the .o and .h files, the second name is the name in this file
--terralib.saveobj("blurExampleTerra.o",{ funcNameInC = blurImageLinComboVectorize4 })

--blurArrayTerra()
--blurArrayTerraVectorizedAcrossX()

--blurImageLinCombo(true)
--blurImageLinCombo(false)
blurImageLinCombo("testingRefactor")
blurImageLinCombo("normal")
--blurImageLinComboVectorize2()
--blurImageLinCombo(false)
--blurImageLinComboVectorize2()
--blurImageLinComboVectorize4()
--blurImageLinComboVectorize8()
--blurImageLinComboVectorize8Wide1()
--blurImageLinComboVectorize16Wide()

--blurMaskedImageLinComboVectorize()
--blurMaskedImageLinComboVectorize8()
----blurMaskedImageOrAllMaskLinComboVectorize()