local C = terralib.includecstring [[
    #include<stdio.h>
    #include<stdlib.h>
    #include<time.h>
    #include<string.h>
    #include<math.h>
    #include <unistd.h>
]]

local boundingBox = 5
local interpDist = 32
local minKernelVal = 0.0000000000000000000000000001
local kernelScaleFactor = 1
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
 	return 0.001f + 0.002f*x + 0.003f*y + 0.004f*x*x + 0.005f*x*y
                     + 0.006f*y*y +  0.0007f*x*x*x + 0.0008f*x*x*y + 0.0009f*x*y*y
                     + 0.01f*y*y*y;
end

local terra polynomial2(x:float, y:float)
	return 1.001f + 1.002f*x + 1.003f*y + 1.004f*x*x + 1.005f*x*y
                     + 1.006f*y*y +  1.0007f*x*x*x + 1.0008f*x*x*y + 1.0009f*x*y*y
                     + 1.01f*y*y*y;
end

local terra polynomial3(x:float, y:float)
	return 2.001f + 2.002f*x + 2.003f*y + 2.004f*x*x + 2.005f*x*y
                     + 2.006f*y*y +  2.0007f*x*x*x + 2.0008f*x*x*y + 2.0009f*x*y*y
                     + 2.01f*y*y*y;
end

local terra polynomial4(x:float, y:float)
	return 3.001f + 3.002f*x + 3.003f*y + 3.004f*x*x + 3.005f*x*y
                     + 3.006f*y*y +  3.0007f*x*x*x + 3.0008f*x*x*y + 3.0009f*x*y*y
                     + 3.01f*y*y*y;
end

local terra polynomial5(x:float, y:float)
	return 4.001f + 4.002f*x + 4.003f*y + 4.004f*x*x + 4.005f*x*y
                     + 4.006f*y*y +  4.0007f*x*x*x + 4.0008f*x*x*y + 4.0009f*x*y*y
                     + 4.01f*y*y*y;
end


--double precision
local terra polynomial1Double(x:double, y:double) : double
 	return 0.001 + 0.002*x + 0.003*y + 0.004*x*x + 0.005*x*y
                     + 0.006*y*y +  0.0007*x*x*x + 0.0008*x*x*y + 0.0009*x*y*y
                     + 0.01*y*y*y;
end

local terra polynomial2Double(x:double, y:double) : double
	return 1.001 + 1.002*x + 1.003*y + 1.004*x*x + 1.005*x*y
                     + 1.006*y*y +  1.0007*x*x*x + 1.0008*x*x*y + 1.0009*x*y*y
                     + 1.01*y*y*y;
end

local terra polynomial3Double(x:double, y:double) : double
	return 2.001 + 2.002*x + 2.003*y + 2.004*x*x + 2.005*x*y
                     + 2.006*y*y +  2.0007*x*x*x + 2.0008*x*x*y + 2.0009*x*y*y
                     + 2.01*y*y*y;
end

local terra polynomial4Double(x:double, y:double) : double
	return 3.001 + 3.002*x + 3.003*y + 3.004*x*x + 3.005*x*y
                     + 3.006*y*y +  3.0007*x*x*x + 3.0008*x*x*y + 3.0009*x*y*y
                     + 3.01*y*y*y;
end

local terra polynomial5Double(x:double, y:double) : double
	return 4.001 + 4.002*x + 4.003*y + 4.004*x*x + 4.005*x*y
                     + 4.006*y*y +  4.0007*x*x*x + 4.0008*x*x*y + 4.0009*x*y*y
                     + 4.01*y*y*y;
end

--local terra polynomial1(x:float, y:float)
-- 	return C.cos(x+y);
--end
--
--local terra polynomial2(x:float, y:float)
-- 	return C.cos(3*x+y);
--end
--
--local terra polynomial3(x:float, y:float)
-- 	return C.cos(x*y);
--end
--
--local terra polynomial4(x:float, y:float)
-- 	return C.cos(x+5*y);
--end
--
--local terra polynomial5(x:float, y:float)
-- 	return C.cos(2*x+3*y);
--end

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


--Multiply kernel values by large number (10000)
--local terra kernel1(x:float, y:float)
--    var sigmaX1 : float = 2.0f;
--    var sigmaY1 : float = 2.0f;
--    var theta1 : float = 0.0f; --//rotation of sigmaX axis
--    return 10000.0 * (C.exp(-((x*C.cos(theta1) +y*C.sin(theta1))*(x*C.cos(theta1) +y*C.sin(theta1)))
--                    /(2*sigmaX1*sigmaX1)) / (C.sqrtf(2*C.M_PI)*sigmaX1))
--                    *(C.exp(-((y*C.cos(theta1) - x*C.sin(theta1))*(y*C.cos(theta1) - x*C.sin(theta1)))
--                    /(2*sigmaY1*sigmaY1)) / (C.sqrtf(2*C.M_PI)*sigmaY1));
--end
--
--
--
--local terra kernel2(x:float, y:float)
--    var sigmaX2 : float = 0.5f;
--    var sigmaY2 : float = 4.0f;
--    var theta2 : float = 0.0f; --//rotation of sigmaX axis
--    return 10000.0 * (C.exp(-((x*C.cos(theta2) +y*C.sin(theta2))*(x*C.cos(theta2) +y*C.sin(theta2)))
--                    /(2*sigmaX2*sigmaX2)) / (C.sqrtf(2*C.M_PI)*sigmaX2))
--                    *(C.exp(-((y*C.cos(theta2) - x*C.sin(theta2))*(y*C.cos(theta2) - x*C.sin(theta2)))
--                    /(2*sigmaY2*sigmaY2)) / (C.sqrtf(2*C.M_PI)*sigmaY2));
--end
--
--local terra kernel3(x:float, y:float)
--    var sigmaX3 : float = 0.5f;
--    var sigmaY3 : float = 4.0f;
--    var theta3 : float = 3.14159f/4; --//rotation of sigmaX axis
--    return 10000.0 * (C.exp(-((x*C.cos(theta3) +y*C.sin(theta3))*(x*C.cos(theta3) +y*C.sin(theta3)))
--                    /(2*sigmaX3*sigmaX3)) / (C.sqrtf(2*C.M_PI)*sigmaX3))
--                    *(C.exp(-((y*C.cos(theta3) - x*C.sin(theta3))*(y*C.cos(theta3) - x*C.sin(theta3)))
--                    /(2*sigmaY3*sigmaY3)) / (C.sqrtf(2*C.M_PI)*sigmaY3));
--end
--
--local terra kernel4(x:float, y:float)
--    var sigmaX4 : float = 0.5f;
--    var sigmaY4 : float = 4.0f;
--    var theta4 : float = 3.14159f/2; --//rotation of sigmaX axis
--    return 10000.0 * (C.exp(-((x*C.cos(theta4) +y*C.sin(theta4))*(x*C.cos(theta4) +y*C.sin(theta4)))
--                    /(2*sigmaX4*sigmaX4)) / (C.sqrtf(2*C.M_PI)*sigmaX4))
--                    *(C.exp(-((y*C.cos(theta4) - x*C.sin(theta4))*(y*C.cos(theta4) - x*C.sin(theta4)))
--                    /(2*sigmaY4*sigmaY4)) / (C.sqrtf(2*C.M_PI)*sigmaY4));
--end
--
--
--local terra kernel5(x:float, y:float)
--    var sigmaX5 : float = 4.0f;
--    var sigmaY5 : float = 4.0f;
--    var theta5 : float = 0.0; --//rotation of sigmaX axis
--    return 10000.0 * (C.exp(-((x*C.cos(theta5) +y*C.sin(theta5))*(x*C.cos(theta5) +y*C.sin(theta5)))
--                    /(2*sigmaX5*sigmaX5)) / (C.sqrtf(2*C.M_PI)*sigmaX5))
--                    *(C.exp(-((y*C.cos(theta5) - x*C.sin(theta5))*(y*C.cos(theta5) - x*C.sin(theta5)))
--                    /(2*sigmaY5*sigmaY5)) / (C.sqrtf(2*C.M_PI)*sigmaY5));
--end
--


--Original
--local terra kernel1(x:float, y:float)
--    var sigmaX1 : float = 2.0f;
--    var sigmaY1 : float = 2.0f;
--    var theta1 : float = 0.0f; --//rotation of sigmaX axis
--    return (C.exp(-((x*C.cos(theta1) +y*C.sin(theta1))*(x*C.cos(theta1) +y*C.sin(theta1)))
--                    /(2*sigmaX1*sigmaX1)) / (C.sqrtf(2*C.M_PI)*sigmaX1))
--                    *(C.exp(-((y*C.cos(theta1) - x*C.sin(theta1))*(y*C.cos(theta1) - x*C.sin(theta1)))
--                    /(2*sigmaY1*sigmaY1)) / (C.sqrtf(2*C.M_PI)*sigmaY1));
--end
--
--
--
--local terra kernel2(x:float, y:float)
--    var sigmaX2 : float = 0.5f;
--    var sigmaY2 : float = 4.0f;
--    var theta2 : float = 0.0f; --//rotation of sigmaX axis
--    return (C.exp(-((x*C.cos(theta2) +y*C.sin(theta2))*(x*C.cos(theta2) +y*C.sin(theta2)))
--                    /(2*sigmaX2*sigmaX2)) / (C.sqrtf(2*C.M_PI)*sigmaX2))
--                    *(C.exp(-((y*C.cos(theta2) - x*C.sin(theta2))*(y*C.cos(theta2) - x*C.sin(theta2)))
--                    /(2*sigmaY2*sigmaY2)) / (C.sqrtf(2*C.M_PI)*sigmaY2));
--end
--
--local terra kernel3(x:float, y:float)
--    var sigmaX3 : float = 0.5f;
--    var sigmaY3 : float = 4.0f;
--    var theta3 : float = 3.14159f/4; --//rotation of sigmaX axis
--    return (C.exp(-((x*C.cos(theta3) +y*C.sin(theta3))*(x*C.cos(theta3) +y*C.sin(theta3)))
--                    /(2*sigmaX3*sigmaX3)) / (C.sqrtf(2*C.M_PI)*sigmaX3))
--                    *(C.exp(-((y*C.cos(theta3) - x*C.sin(theta3))*(y*C.cos(theta3) - x*C.sin(theta3)))
--                    /(2*sigmaY3*sigmaY3)) / (C.sqrtf(2*C.M_PI)*sigmaY3));
--end
--
--local terra kernel4(x:float, y:float)
--    var sigmaX4 : float = 0.5f;
--    var sigmaY4 : float = 4.0f;
--    var theta4 : float = 3.14159f/2; --//rotation of sigmaX axis
--    return (C.exp(-((x*C.cos(theta4) +y*C.sin(theta4))*(x*C.cos(theta4) +y*C.sin(theta4)))
--                    /(2*sigmaX4*sigmaX4)) / (C.sqrtf(2*C.M_PI)*sigmaX4))
--                    *(C.exp(-((y*C.cos(theta4) - x*C.sin(theta4))*(y*C.cos(theta4) - x*C.sin(theta4)))
--                    /(2*sigmaY4*sigmaY4)) / (C.sqrtf(2*C.M_PI)*sigmaY4));
--end


--local terra kernel5(x:float, y:float)
--    var sigmaX5 : float = 4.0f;
--    var sigmaY5 : float = 4.0f;
--    var theta5 : float = 0.0; --//rotation of sigmaX axis
--    return (C.exp(-((x*C.cos(theta5) +y*C.sin(theta5))*(x*C.cos(theta5) +y*C.sin(theta5)))
--                    /(2*sigmaX5*sigmaX5)) / (C.sqrtf(2*C.M_PI)*sigmaX5))
--                    *(C.exp(-((y*C.cos(theta5) - x*C.sin(theta5))*(y*C.cos(theta5) - x*C.sin(theta5)))
--                    /(2*sigmaY5*sigmaY5)) / (C.sqrtf(2*C.M_PI)*sigmaY5));
--end


--double precision kernels
local terra kernel1Double(x:float, y:float)
    var sigmaX1 : float = 2.0f;
    var sigmaY1 : float = 2.0f;
    var theta1 : float = 0.0f; --//rotation of sigmaX axis
    return (C.exp(-((x*C.cos(theta1) +y*C.sin(theta1))*(x*C.cos(theta1) +y*C.sin(theta1)))
                    /(2*sigmaX1*sigmaX1)) / (C.sqrt(2*C.M_PI)*sigmaX1))
                    *(C.exp(-((y*C.cos(theta1) - x*C.sin(theta1))*(y*C.cos(theta1) - x*C.sin(theta1)))
                    /(2*sigmaY1*sigmaY1)) / (C.sqrt(2*C.M_PI)*sigmaY1));
end

local terra kernel2Double(x:float, y:float)
    var sigmaX2 : float = 0.5f;
    var sigmaY2 : float = 4.0f;
    var theta2 : float = 0.0f; --//rotation of sigmaX axis
    return (C.exp(-((x*C.cos(theta2) +y*C.sin(theta2))*(x*C.cos(theta2) +y*C.sin(theta2)))
                    /(2*sigmaX2*sigmaX2)) / (C.sqrt(2*C.M_PI)*sigmaX2))
                    *(C.exp(-((y*C.cos(theta2) - x*C.sin(theta2))*(y*C.cos(theta2) - x*C.sin(theta2)))
                    /(2*sigmaY2*sigmaY2)) / (C.sqrt(2*C.M_PI)*sigmaY2));
end

local terra kernel3Double(x:float, y:float)
    var sigmaX3 : float = 0.5f;
    var sigmaY3 : float = 4.0f;
    var theta3 : float = 3.14159f/4; --//rotation of sigmaX axis
    return (C.exp(-((x*C.cos(theta3) +y*C.sin(theta3))*(x*C.cos(theta3) +y*C.sin(theta3)))
                    /(2*sigmaX3*sigmaX3)) / (C.sqrt(2*C.M_PI)*sigmaX3))
                    *(C.exp(-((y*C.cos(theta3) - x*C.sin(theta3))*(y*C.cos(theta3) - x*C.sin(theta3)))
                    /(2*sigmaY3*sigmaY3)) / (C.sqrt(2*C.M_PI)*sigmaY3));
end

local terra kernel4Double(x:float, y:float)
    var sigmaX4 : float = 0.5f;
    var sigmaY4 : float = 4.0f;
    var theta4 : float = 3.14159f/2; --//rotation of sigmaX axis
    return (C.exp(-((x*C.cos(theta4) +y*C.sin(theta4))*(x*C.cos(theta4) +y*C.sin(theta4)))
                    /(2*sigmaX4*sigmaX4)) / (C.sqrt(2*C.M_PI)*sigmaX4))
                    *(C.exp(-((y*C.cos(theta4) - x*C.sin(theta4))*(y*C.cos(theta4) - x*C.sin(theta4)))
                    /(2*sigmaY4*sigmaY4)) / (C.sqrt(2*C.M_PI)*sigmaY4));
end

local terra kernel5Double(x:float, y:float)
    var sigmaX5 : float = 4.0f;
    var sigmaY5 : float = 4.0f;
    var theta5 : float = 0.0; --//rotation of sigmaX axis
    return (C.exp(-((x*C.cos(theta5) +y*C.sin(theta5))*(x*C.cos(theta5) +y*C.sin(theta5)))
                    /(2*sigmaX5*sigmaX5)) / (C.sqrt(2*C.M_PI)*sigmaX5))
                    *(C.exp(-((y*C.cos(theta5) - x*C.sin(theta5))*(y*C.cos(theta5) - x*C.sin(theta5)))
                    /(2*sigmaY5*sigmaY5)) / (C.sqrt(2*C.M_PI)*sigmaY5));
end






--round to 0 if kernel value is less than minKernelVal
local terra kernel1(x:float, y:float)
    var sigmaX1 : float = 2.0f;
    var sigmaY1 : float = 2.0f;
    var theta1 : float = 0.0f; --//rotation of sigmaX axis
    var returnVal = kernelScaleFactor * (C.exp(-((x*C.cos(theta1) +y*C.sin(theta1))*(x*C.cos(theta1) +y*C.sin(theta1)))
                    /(2*sigmaX1*sigmaX1)) / (C.sqrtf(2*C.M_PI)*sigmaX1))
                    *(C.exp(-((y*C.cos(theta1) - x*C.sin(theta1))*(y*C.cos(theta1) - x*C.sin(theta1)))
                    /(2*sigmaY1*sigmaY1)) / (C.sqrtf(2*C.M_PI)*sigmaY1));
    if(returnVal > minKernelVal) then
    	return returnVal
   	else
   		return 0.0
   	end
end


local terra kernel2(x:float, y:float)
    var sigmaX2 : float = 0.5f;
    var sigmaY2 : float = 4.0f;
    var theta2 : float = 0.0f; --//rotation of sigmaX axis
    var returnVal = kernelScaleFactor * (C.exp(-((x*C.cos(theta2) +y*C.sin(theta2))*(x*C.cos(theta2) +y*C.sin(theta2)))
                    /(2*sigmaX2*sigmaX2)) / (C.sqrtf(2*C.M_PI)*sigmaX2))
                    *(C.exp(-((y*C.cos(theta2) - x*C.sin(theta2))*(y*C.cos(theta2) - x*C.sin(theta2)))
                    /(2*sigmaY2*sigmaY2)) / (C.sqrtf(2*C.M_PI)*sigmaY2));
    if(returnVal > minKernelVal) then
    	return returnVal
   	else
   		return 0.0
   	end
end


local terra kernel3(x:float, y:float)
    var sigmaX3 : float = 0.5f;
    var sigmaY3 : float = 4.0f;
    var theta3 : float = 3.14159f/4; --//rotation of sigmaX axis
    var returnVal = kernelScaleFactor * (C.exp(-((x*C.cos(theta3) +y*C.sin(theta3))*(x*C.cos(theta3) +y*C.sin(theta3)))
                    /(2*sigmaX3*sigmaX3)) / (C.sqrtf(2*C.M_PI)*sigmaX3))
                    *(C.exp(-((y*C.cos(theta3) - x*C.sin(theta3))*(y*C.cos(theta3) - x*C.sin(theta3)))
                    /(2*sigmaY3*sigmaY3)) / (C.sqrtf(2*C.M_PI)*sigmaY3));
    if(returnVal > minKernelVal) then
    	return returnVal
   	else
   		return 0.0
   	end
end


local terra kernel4(x:float, y:float)
    var sigmaX4 : float = 0.5f;
    var sigmaY4 : float = 4.0f;
    var theta4 : float = 3.14159f/2; --//rotation of sigmaX axis
    var returnVal = kernelScaleFactor * (C.exp(-((x*C.cos(theta4) +y*C.sin(theta4))*(x*C.cos(theta4) +y*C.sin(theta4)))
                    /(2*sigmaX4*sigmaX4)) / (C.sqrtf(2*C.M_PI)*sigmaX4))
                    *(C.exp(-((y*C.cos(theta4) - x*C.sin(theta4))*(y*C.cos(theta4) - x*C.sin(theta4)))
                    /(2*sigmaY4*sigmaY4)) / (C.sqrtf(2*C.M_PI)*sigmaY4));
    if(returnVal > minKernelVal) then
    	return returnVal
   	else
   		return 0.0
   	end
end


local terra kernel5(x:float, y:float)
    var sigmaX5 : float = 4.0f;
    var sigmaY5 : float = 4.0f;
    var theta5 : float = 0.0; --//rotation of sigmaX axis
    var returnVal = kernelScaleFactor * (C.exp(-((x*C.cos(theta5) +y*C.sin(theta5))*(x*C.cos(theta5) +y*C.sin(theta5)))
                    /(2*sigmaX5*sigmaX5)) / (C.sqrtf(2*C.M_PI)*sigmaX5))
                    *(C.exp(-((y*C.cos(theta5) - x*C.sin(theta5))*(y*C.cos(theta5) - x*C.sin(theta5)))
                    /(2*sigmaY5*sigmaY5)) / (C.sqrtf(2*C.M_PI)*sigmaY5));
    if(returnVal > minKernelVal) then
    	return returnVal
   	else
   		return 0.0
   	end
end


--test with kernels that are more orthogonal
--local terra kernel1(x:float, y:float)
--	return C.fabs(C.cos(x+y)/10.0f)
--end
--
--local terra kernel2(x:float, y:float)
--	return C.fabs(C.sin(x+y)/50.0f)
--end
--
--local terra kernel3(x:float, y:float)
--	return C.fabs(C.sin(7*x*y)/30.0f)
--end
--
--local terra kernel4(x:float, y:float)
--	return C.fabs(C.cos(x*y)/20.0f)
--end
--
--
--local terra kernel5(x:float, y:float)
--	return C.fabs(C.cos(5*x+11*y)/80.0f)
--end


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
		C.printf("\n\nImage plane ONLY lin combo %dx%d blur Terra, method = %s:\n", 2*boundingBox+1, 2*boundingBox+1, method)
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



local function blurImageLinComboNoUnroll()
	local terra terraBlur()
		var width : int = 2048
		var height : int = 1489

		var input : &float
		input = [&float](C.malloc(sizeof(float)*width*height))

		var output : &float
		output = [&float](C.malloc(sizeof(float)*width*height))

		var kernelArea = (boundingBox*2+1)*(boundingBox*2+1)
		var kernelWidth = boundingBox*2+1
		var kernelVals = [&float](C.malloc(kernelArea*5*sizeof(float)))


	for x = 0, width do
		for y = 0, height do
			input[y*width + x] = [float](y*x*C.cos(x/(y+1)))

			output[y*width + x] = 0.0f
		end
	end

		var t1 = C.clock()

		for y = boundingBox, height-boundingBox do
			for x = boundingBox, width-boundingBox do
				var curOut : float = 0
				var curNorm : float = 0
				var curKernelVal : float


				escape
					for j = -boundingBox, boundingBox do
						for i = -boundingBox, boundingBox do
							emit quote
								kernelVals[(kernelWidth*(j+boundingBox) + (i+boundingBox))*5] = [luaKernel1(i, j)]
								kernelVals[(kernelWidth*(j+boundingBox) + (i+boundingBox))*5+1] = [luaKernel2(i, j)]
								kernelVals[(kernelWidth*(j+boundingBox) + (i+boundingBox))*5+2] = [luaKernel3(i, j)]
								kernelVals[(kernelWidth*(j+boundingBox) + (i+boundingBox))*5+3] = [luaKernel4(i, j)]
								kernelVals[(kernelWidth*(j+boundingBox) + (i+boundingBox))*5+4] = [luaKernel5(i, j)]
							end
						end
					end
				end

				--escape 
				
				    for j = -boundingBox, boundingBox + 1 do
				        for i = -boundingBox, boundingBox + 1 do
				        	--emit quote
				        		var curKernelLocation = (kernelWidth*(j+boundingBox) + (i+boundingBox))*5
				        		curKernelVal = polynomial1(x, y)*kernelVals[curKernelLocation] +
					                polynomial2(x, y)*kernelVals[curKernelLocation+1] + polynomial3(x, y)*kernelVals[curKernelLocation+2] + 
					                polynomial4(x, y)*kernelVals[curKernelLocation+3] + polynomial5(x, y)*kernelVals[curKernelLocation+4];
					            curOut = curOut + input[(x + i) + width*(y + j)]*curKernelVal ; 
					            curNorm = curNorm + curKernelVal;
					        --end
				        end
				    end
					
				--end
			    output[x + width*y] = curOut/curNorm
			end
		end


		var t2 = C.clock()


		for i=boundingBox,boundingBox+10 do
			for j=boundingBox,boundingBox+10 do
				C.printf("%f\t", output[i*width + j])
			end
			C.printf("\n")
		end

		C.printf("\n\nImage only, no unroll, lin combo %dx%d blur Terra:\n", 2*boundingBox+1, 2*boundingBox+1)
		C.printf("output[boundingBox*width + boundingBox] = %f, computation took: %f ms\n", output[boundingBox*width + boundingBox],  (t2-t1)/1000.0f)
		C.printf("C.CLOCKS_PER_SEC = %d\n", C.CLOCKS_PER_SEC)

	    C.free(input)
		C.free(output)
		C.free(kernelVals)
	end
--	terraBlur:disas()
	terraBlur()
end

local function blurImageLinComboSplitX()
	local terra terraBlur()
		var width : int = 2048
		var height : int = 1489

		var input : &float
		input = [&float](C.malloc(sizeof(float)*width*height))

		var output : &float
		output = [&float](C.malloc(sizeof(float)*width*height))

	for x = 0, width do
		for y = 0, height do
			input[y*width + x] = [float](y*x*C.cos(x/(y+1)))

			output[y*width + x] = 0.0f
		end
	end

		var t1 = C.clock()

		var xSplit : int = 50
		for xOuter = boundingBox, width - boundingBox - xSplit, xSplit do
			for y = boundingBox, height-boundingBox do
				for x = xOuter, xOuter+xSplit do
					var curOut : float = 0
					var curNorm : float = 0
					var curKernelVal : float

					escape 
					
					    for j = -boundingBox, boundingBox do
					        for i = -boundingBox, boundingBox do
					        	emit quote
					        		curKernelVal = polynomial1(x, y)*kernel1(i, j) +
						                polynomial2(x, y)*kernel2(i, j) + polynomial3(x, y)*kernel3(i, j) + 
						                polynomial4(x, y)*kernel4(i, j) + polynomial5(x, y)*kernel5(i, j);
						            curOut = curOut + input[(x + i) + width*(y + j)]*curKernelVal; 

						            curNorm = curNorm + curKernelVal;
						        end
					        end
					    end
						
					end
				    output[x + width*y] = curOut/(curNorm)
				end
			end
		end

		var t2 = C.clock()



		C.printf("\n\nImage only, x dimension split by %d, lin combo %dx%d blur Terra:\n", xSplit, 2*boundingBox+1, 2*boundingBox+1)
		C.printf("output[boundingBox*width + boundingBox] = %f, computation took: %f ms\n", output[boundingBox*width + boundingBox],  (t2-t1)/1000.0f)
		C.printf("C.CLOCKS_PER_SEC = %d\n", C.CLOCKS_PER_SEC)

		for i=boundingBox,boundingBox+10 do
			for j=boundingBox,boundingBox+10 do
				C.printf("%f\t", output[i*width + j])
			end
			C.printf("\n")
		end

	    C.free(input)
		C.free(output)
	end
--	terraBlur:disas()
	terraBlur()
end

local function blurVarianceLinCombo()
	local terra terraBlur()
		var width : int = 2048
		var height : int = 1489

		var input : &float
		input = [&float](C.malloc(sizeof(float)*width*height))

		var output : &float
		output = [&float](C.malloc(sizeof(float)*width*height))

	for x = 0, width do
		for y = 0, height do
			input[y*width + x] = [float](y*x*C.cos(x/(y+1)))

			output[y*width + x] = 0.0f
		end
	end

		var t1 = C.clock()

		for y = boundingBox, height-boundingBox do
			for x = boundingBox, width-boundingBox do
				var curOut : float = 0
				var curNorm : float = 0
				var curKernelVal : float

				escape 
				
				    for j = -boundingBox, boundingBox do
				        for i = -boundingBox, boundingBox do
				        	emit quote
				        		curKernelVal = polynomial1(x, y)*kernel1(i, j) +
					                polynomial2(x, y)*kernel2(i, j) + polynomial3(x, y)*kernel3(i, j) + 
					                polynomial4(x, y)*kernel4(i, j) + polynomial5(x, y)*kernel5(i, j);
					            curOut = curOut + input[(x + i) + width*(y + j)]*curKernelVal*curKernelVal ; 

					            curNorm = curNorm + curKernelVal;
					        end
				        end
				    end
					
				end
			    output[x + width*y] = curOut/(curNorm*curNorm)
			end
		end


		var t2 = C.clock()



		C.printf("\n\nVariance only, no unroll, lin combo %dx%d blur Terra:\n", 2*boundingBox+1, 2*boundingBox+1)
		C.printf("output[boundingBox*width + boundingBox] = %f, computation took: %f ms\n", output[boundingBox*width + boundingBox],  (t2-t1)/1000.0f)
		C.printf("C.CLOCKS_PER_SEC = %d\n", C.CLOCKS_PER_SEC)

		for i=boundingBox,boundingBox+10 do
			for j=boundingBox,boundingBox+10 do
				C.printf("%f\t", output[i*width + j])
			end
			C.printf("\n")
		end

	    C.free(input)
		C.free(output)
	end
--	terraBlur:disas()
	terraBlur()
end

local function blurMaskPlaneLinCombo()
	local terra terraBlur()
		var width : int = 2048
		var height : int = 1489

		var input : &uint16
		input = [&uint16](C.malloc(sizeof(uint16)*width*height))

		var output : &uint16
		output = [&uint16](C.malloc(sizeof(uint16)*width*height))

	for x = 0, width do
		for y = 0, height do
			input[y*width + x] = [uint16](x*C.cos(x*y) + y)

			output[y*width + x] = 0
		end
	end

		var t1 = C.clock()

		for y = boundingBox, height-boundingBox do
			for x = boundingBox, width-boundingBox do
				var curOut : uint16 = 0
				var curKernelVal : float

				escape 
				
				    for j = -boundingBox, boundingBox do
				        for i = -boundingBox, boundingBox do
				        	emit quote
				        		curKernelVal = polynomial1(x, y)*kernel1(i, j) +
					                polynomial2(x, y)*kernel2(i, j) + polynomial3(x, y)*kernel3(i, j) + 
					                polynomial4(x, y)*kernel4(i, j) + polynomial5(x, y)*kernel5(i, j);

					            if curKernelVal ~= 0.0f then
					            	curOut = curOut or input[(x + i) + width*(y + j)] 
					            end

					        end
				        end
				    end
					
				end
			    output[x + width*y] = curOut
			end
		end


		var t2 = C.clock()


		for i=boundingBox,boundingBox+10 do
			for j=boundingBox,boundingBox+10 do
				C.printf("%d\t", output[i*width + j])
			end
			C.printf("\n")
		end

		C.printf("\n\nMask only, lin combo %dx%d blur Terra:\n", 2*boundingBox+1, 2*boundingBox+1)
		C.printf("output[boundingBox*width + boundingBox] = %d, computation took: %f ms\n", output[boundingBox*width + boundingBox],  (t2-t1)/1000.0f)
		C.printf("C.CLOCKS_PER_SEC = %d\n", C.CLOCKS_PER_SEC)

	    C.free(input)
		C.free(output)
	end
--	terraBlur:disas()
	terraBlur()
end



--try only changing the mask plane correctly, but adding one to every
--other plane
local symbolTable = {}
local function blurMaskPlaneLinComboTest()
	local terra terraBlur()
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
			inputVar[y*width + x] = [float](y*x*C.cos(x/(y+1)))
			inputMask[y*width + x] = [uint16](x*C.cos(x*y) + y)


			outputIm[y*width + x] = 0.0f
			outputVar[y*width + x] = 0.0f
			outputMask[y*width + x] = 0
		end
	end

		var t1 = C.clock()


		escape

			local luaKernelArea = (boundingBox*2+1)*(boundingBox*2+1)
			local luaKernelWidth = boundingBox*2+1

	        for i = 1, luaKernelArea*2 do
	            local cur_kernelValue_symbol = symbol(float, "kVal"..i)
	            table.insert(symbolTable, cur_kernelValue_symbol)
	        end

			for j = -boundingBox, boundingBox do
				for i = -boundingBox, boundingBox do
					emit quote
						var [symbolTable[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*2+1]] = inputIm[( i) + width*( j)]
						var [symbolTable[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*2+2]] = inputVar[( i) + width*( j)]
					end
				end
			end
		end



		for y = boundingBox, height-boundingBox do
			for x = boundingBox, width-boundingBox do
				var curOut : uint16 = 0
				var curKernelVal : float

				var imageOut = inputIm[x + width*y]
				var varOut = inputVar[x + width*y]
				escape 
					local luaKernelWidth = boundingBox*2+1

				    for j = -boundingBox, boundingBox do
				        for i = -boundingBox, boundingBox do
				        	emit quote
				        		curKernelVal = polynomial1(x, y)*kernel1(i, j) +
					                polynomial2(x, y)*kernel2(i, j) + polynomial3(x, y)*kernel3(i, j) + 
					                polynomial4(x, y)*kernel4(i, j) + polynomial5(x, y)*kernel5(i, j);

					            if curKernelVal ~= 0.0f then
					            	curOut = curOut or inputMask[(x + i) + width*(y + j)] 
					            end

					            --imageOut = imageOut + inputIm[( i) + width*( j)]
					            --varOut = varOut + inputVar[( i) + width*( j)]
					            --imageOut = imageOut + [symbolTable[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*2+1]]
					            --varOut = varOut + [symbolTable[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*2+2]]

					        end
				        end
				    end
					
				end
			    outputMask[x + width*y] = curOut
			    outputIm[x + width*y] = imageOut
			    outputVar[x + width*y] = varOut
			end
		end


		var t2 = C.clock()

		C.printf("mask:\n")
		for i=boundingBox,boundingBox+10 do
			for j=boundingBox,boundingBox+10 do
				C.printf("%d\t", outputMask[i*width + j])
			end
			C.printf("\n")
		end
		C.printf("image:\n")
		for i=boundingBox,boundingBox+10 do
			for j=boundingBox,boundingBox+10 do
				C.printf("%d\t", outputIm[i*width + j])
			end
			C.printf("\n")
		end
		C.printf("Variance:\n")
		for i=boundingBox,boundingBox+10 do
			for j=boundingBox,boundingBox+10 do
				C.printf("%d\t", outputVar[i*width + j])
			end
			C.printf("\n")
		end

		C.printf("\n\nMask only, ADD TO OthER PLANES, lin combo %dx%d blur Terra:\n", 2*boundingBox+1, 2*boundingBox+1)
		C.printf("outputMask[boundingBox*width + boundingBox] = %d, computation took: %f ms\n", outputMask[boundingBox*width + boundingBox],  (t2-t1)/1000.0f)
		C.printf("C.CLOCKS_PER_SEC = %d\n", C.CLOCKS_PER_SEC)

		C.free(inputIm)
		C.free(outputIm)
		C.free(inputVar)
		C.free(outputVar)
	    C.free(inputMask)
		C.free(outputMask)
	end
--	terraBlur:disas()
	terraBlur()
end



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
	C.printf("\n\nVectorized2 lin combo %dx%d blur Terra:\n", 2*boundingBox+1, 2*boundingBox+1)
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
	C.printf("\n\nVectorized84 lin combo %dx%d blur Terra:\n", 2*boundingBox+1, 2*boundingBox+1)
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
	C.printf("\n\nVectorized8 lin combo %dx%d blur Terra:\n", 2*boundingBox+1, 2*boundingBox+1)
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
	C.printf("\n\nVectorized 8wide lin combo %dx%d blur Terra:\n", 2*boundingBox+1, 2*boundingBox+1)
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
	C.printf("\n\nVectorized 16wide lin combo %dx%d blur Terra:\n", 2*boundingBox+1, 2*boundingBox+1)
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

	for x = 0, width do
		for y = 0, height do
			inputIm[y*width + x] = [float](y*x*C.cos(x/(y+1)))
			inputVar[y*width + x] = [float](y*x*C.cos(x/(y+1)))
			inputMask[y*width + x] = [uint16](x*C.cos(x*y) + y)


			outputIm[y*width + x] = 0.0f
			outputVar[y*width + x] = 0.0f
			outputMask[y*width + x] = 0
		end
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
	C.printf("\n\nMasked image lin combo %dx%d blur Terra:\n", 2*boundingBox+1, 2*boundingBox+1)
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


--deal with image, mask, and variance planes
local terra blurMaskedImageLinComboSplitX()
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
			inputVar[y*width + x] = [float](y*x*C.cos(x/(y+1)))
			inputMask[y*width + x] = [uint16](x*C.cos(x*y) + y)


			outputIm[y*width + x] = 0.0f
			outputVar[y*width + x] = 0.0f
			outputMask[y*width + x] = 0
		end
	end

	var t1 = C.clock()

	var xSplit : int = 50
	for xOuter = boundingBox, (width-boundingBox-xSplit), xSplit do
		for y = boundingBox, height-boundingBox do
			for x = xOuter, xOuter+xSplit do
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
	end

	var t2 = C.clock()
	C.printf("\n\nMasked image lin combo, x split by %d, %dx%d blur Terra:\n", xSplit, 2*boundingBox+1, 2*boundingBox+1)
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


--playing around with aligning arrays to cache lines
--deal with image, mask, and variance planes
local terra blurMaskedImageLinComboNoUnroll()
--	var boundingBox : int = 2
	var width : int = 2048
	var height : int = 1489


--	var inputIm : &float
--	inputIm = [&float](C.malloc(sizeof(float)*width*height))
--	var inputVar : &float
--	inputVar = [&float](C.malloc(sizeof(float)*width*height))
--	var inputMask : &uint16
--	inputMask = [&uint16](C.malloc(sizeof(uint16)*width*height))
--
--	var outputIm : &float
--	outputIm = [&float](C.malloc(sizeof(float)*width*height))
--	var outputVar : &float
--	outputVar = [&float](C.malloc(sizeof(float)*width*height))
--	var outputMask : &uint16
--	outputMask = [&uint16](C.malloc(sizeof(uint16)*width*height))

	var tempPointer : &opaque

	var inputIm : &float
	C.posix_memalign(&tempPointer, 64, sizeof(float)*width*height)
	inputIm = [&float](tempPointer)

	var inputVar : &float
	C.posix_memalign(&tempPointer, 64, sizeof(float)*width*height)
	inputVar = [&float](tempPointer)

	var inputMask : &uint16
	C.posix_memalign(&tempPointer, 64, sizeof(uint16)*width*height)
	inputMask = [&uint16](tempPointer)



	var outputIm : &float
	C.posix_memalign(&tempPointer, 64, sizeof(float)*width*height)
	outputIm = [&float](tempPointer)

	var outputVar : &float
	C.posix_memalign(&tempPointer, 64, sizeof(float)*width*height)
	outputVar = [&float](tempPointer)

	var outputMask : &uint16
	C.posix_memalign(&tempPointer, 64, sizeof(uint16)*width*height)
	outputMask = [&uint16](tempPointer)

	var kernelVals : &float

	for x = 0, width do
		for y = 0, height do
			inputIm[y*width + x] = [float](y*x*C.cos(x/(y+1)))
			inputVar[y*width + x] = [float](y*x*C.cos(x/(y+1)))
			inputMask[y*width + x] = [uint16](x*C.cos(x*y) + y)


			outputIm[y*width + x] = 0.0f
			outputVar[y*width + x] = 0.0f
			outputMask[y*width + x] = 0
		end
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

				var kernelArea = (boundingBox*2+1)*(boundingBox*2+1)
				var kernelWidth = boundingBox*2+1
				--kernelVals = [&float](C.malloc(kernelArea*5*sizeof(float)))



				--C.posix_memalign(&tempPointer, C.sysconf(C._SC_PAGESIZE), kernelArea*5*sizeof(float))
				C.posix_memalign(&tempPointer, 64, kernelArea*5*sizeof(float))
				kernelVals = [&float](tempPointer)

				escape
					for j = -boundingBox, boundingBox do
						for i = -boundingBox, boundingBox do
							emit quote
								kernelVals[(kernelWidth*(j+boundingBox) + (i+boundingBox))*5] = [kernel1(i, j)]
								kernelVals[(kernelWidth*(j+boundingBox) + (i+boundingBox))*5+1] = [kernel2(i, j)]
								kernelVals[(kernelWidth*(j+boundingBox) + (i+boundingBox))*5+2] = [kernel3(i, j)]
								kernelVals[(kernelWidth*(j+boundingBox) + (i+boundingBox))*5+3] = [kernel4(i, j)]
								kernelVals[(kernelWidth*(j+boundingBox) + (i+boundingBox))*5+4] = [kernel5(i, j)]
							end
						end
					end
				end

				escape 
				    for j = -boundingBox, boundingBox do
				        for i = -boundingBox, boundingBox do
				        	emit quote
				--    for j = -boundingBox, boundingBox+1 do
				--        for i = -boundingBox, boundingBox+1 do
				        		var curKernelLocation = (kernelWidth*(j+boundingBox) + (i+boundingBox))*5
				        		curKernelVal = polynomial1(x, y)*kernelVals[curKernelLocation] +
					                polynomial2(x, y)*kernelVals[curKernelLocation+1] + polynomial3(x, y)*kernelVals[curKernelLocation+2] + 
					                polynomial4(x, y)*kernelVals[curKernelLocation+3] + polynomial5(x, y)*kernelVals[curKernelLocation+4];

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
	C.printf("\n\nMasked image lin combo, NO UNROLLING %dx%d blur Terra:\n", 2*boundingBox+1, 2*boundingBox+1)
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
	C.free(kernelVals)
end



--deal with image, mask, and variance planes
local kernelValSymbols = {}

local function luaBlurMaskedImageLinComboUnrollSymbols(_boundingBox)
	local terra blurMaskedImageLinComboUnrollSymbols()
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

		var kernelVals : &float

		for x = 0, width do
			for y = 0, height do
				inputIm[y*width + x] = [float](y*x*C.cos(x/(y+1)))
				inputVar[y*width + x] = [float](y*x*C.cos(x/(y+1)))
				inputMask[y*width + x] = [uint16](x*C.cos(x*y) + y)


				outputIm[y*width + x] = 0.0f
				outputVar[y*width + x] = 0.0f
				outputMask[y*width + x] = 0
			end
		end

		var t1 = C.clock()

		escape

			local luaKernelArea = (_boundingBox*2+1)*(_boundingBox*2+1)
			local luaKernelWidth = _boundingBox*2+1

	        for i = 1, luaKernelArea*5 do
	            local cur_kernelValue_symbol = symbol(float, "kVal"..i)
	            table.insert(kernelValSymbols, cur_kernelValue_symbol)
	        end


			for j = -_boundingBox, _boundingBox do
				for i = -_boundingBox, _boundingBox do
					emit quote
						var [kernelValSymbols[(luaKernelWidth*(j+_boundingBox) + (i+_boundingBox))*5+1]] = [kernel1(i, j)]
						var [kernelValSymbols[(luaKernelWidth*(j+_boundingBox) + (i+_boundingBox))*5+2]] = [kernel2(i, j)]
						var [kernelValSymbols[(luaKernelWidth*(j+_boundingBox) + (i+_boundingBox))*5+3]] = [kernel3(i, j)]
						var [kernelValSymbols[(luaKernelWidth*(j+_boundingBox) + (i+_boundingBox))*5+4]] = [kernel4(i, j)]
						var [kernelValSymbols[(luaKernelWidth*(j+_boundingBox) + (i+_boundingBox))*5+5]] = [kernel5(i, j)]
					end
				end
			end
		end

		var tB = C.clock()
		C.printf("time spent calculating kernel = %f ms!!!!!!!!!!!!!!@@@@@@@@@@@@", (float)(tB-t1)/C.CLOCKS_PER_SEC*1000.0f)



	--	var xSplit : int = 32
	--	for xOuter = _boundingBox, (width-_boundingBox-xSplit), xSplit do
			for y = _boundingBox, height-_boundingBox do
				for x = _boundingBox, (width-_boundingBox) do
	--			for x = xOuter, xOuter+xSplit, 4 do
					var curKernelVal : float

					var curImOut : float = 0.0f
					var curVarOut : float = 0.0f
					var curMaskOut : uint16 = 0

					var curNorm : float = 0.0f

					var kernelArea = (_boundingBox*2+1)*(_boundingBox*2+1)
					var kernelWidth = _boundingBox*2+1
	--				kernelVals = [&float](C.malloc(kernelArea*5*sizeof(float)))


					escape 
						local luaKernelWidth = _boundingBox*2+1

					    for j = -_boundingBox, _boundingBox do
					        for i = -_boundingBox, _boundingBox do
					        	emit quote
					--    for j = -_boundingBox, _boundingBox+1 do
					--        for i = -_boundingBox, _boundingBox+1 do
					--        		var curKernelLocation = (kernelWidth*(j+_boundingBox) + (i+_boundingBox))*5
					        		curKernelVal = polynomial1(x, y)*[kernelValSymbols[(luaKernelWidth*(j+_boundingBox) + (i+_boundingBox))*5+1]] +
						                polynomial2(x, y)*[kernelValSymbols[(luaKernelWidth*(j+_boundingBox) + (i+_boundingBox))*5+2]] + polynomial3(x, y)*[kernelValSymbols[(luaKernelWidth*(j+_boundingBox) + (i+_boundingBox))*5+3]] + 
						                polynomial4(x, y)*[kernelValSymbols[(luaKernelWidth*(j+_boundingBox) + (i+_boundingBox))*5+4]] + polynomial5(x, y)*[kernelValSymbols[(luaKernelWidth*(j+_boundingBox) + (i+_boundingBox))*5+5]];

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
		C.printf("\n\nMasked image lin combo, UNROLLING w/ SYMBOLS %dx%d blur Terra:\n", 2*_boundingBox+1, 2*_boundingBox+1)
		C.printf("outputIm[_boundingBox*width + _boundingBox] = %f, computation took: %f ms\n", outputIm[_boundingBox*width + _boundingBox],  (t2-t1)/1000.0)
		C.printf("C.CLOCKS_PER_SEC = %d\n", C.CLOCKS_PER_SEC)

		C.printf("Image plane, 10x10 box begining at (_boundingBox,_boundingBox)\n")
		for i=_boundingBox,_boundingBox+10 do
			for j=_boundingBox,_boundingBox+10 do
				C.printf("%f\t", outputIm[i*width + j])
			end
			C.printf("\n")
		end

		C.printf("Variance plane, 10x10 box begining at (_boundingBox,_boundingBox)\n")
		for i=_boundingBox,_boundingBox+10 do
			for j=_boundingBox,_boundingBox+10 do
				C.printf("%f\t", outputVar[i*width + j])
			end
			C.printf("\n")
		end

		C.printf("Mask plane, 10x10 box begining at (_boundingBox,_boundingBox)\n")
		for i=_boundingBox,_boundingBox+10 do
			for j=_boundingBox,_boundingBox+10 do
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
	--	C.free(kernelVals)
	end
	blurMaskedImageLinComboUnrollSymbols()
end


--deal with image, mask, and variance planes
local kernelValSymbols6 = {}

local terra blurMaskedImageLinComboUnrollSymbolsDoublePrecisionKernel()
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

	var kernelVals : &float

	for x = 0, width do
		for y = 0, height do
			inputIm[y*width + x] = [float](y*x*C.cos(x/(y+1)))
			inputVar[y*width + x] = [float](y*x*C.cos(x/(y+1)))
			inputMask[y*width + x] = [uint16](x*C.cos(x*y) + y)


			outputIm[y*width + x] = 0.0f
			outputVar[y*width + x] = 0.0f
			outputMask[y*width + x] = 0
		end
	end

	var t1 = C.clock()

	escape

		local luaKernelArea = (boundingBox*2+1)*(boundingBox*2+1)
		local luaKernelWidth = boundingBox*2+1

        for i = 1, luaKernelArea*5 do
            local cur_kernelValue_symbol = symbol(double, "kVal"..i)
            table.insert(kernelValSymbols6, cur_kernelValue_symbol)
        end


		for j = -boundingBox, boundingBox do
			for i = -boundingBox, boundingBox do
				emit quote
					var [kernelValSymbols6[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+1]] = [kernel1Double(i, j)]
					var [kernelValSymbols6[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+2]] = [kernel2Double(i, j)]
					var [kernelValSymbols6[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+3]] = [kernel3Double(i, j)]
					var [kernelValSymbols6[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+4]] = [kernel4Double(i, j)]
					var [kernelValSymbols6[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+5]] = [kernel5Double(i, j)]
				end
			end
		end
	end

	var tB = C.clock()
	C.printf("time spent calculating kernel = %f ms!!!!!!!!!!!!!!@@@@@@@@@@@@", (float)(tB-t1)/C.CLOCKS_PER_SEC*1000.0f)



--	var xSplit : int = 32
--	for xOuter = boundingBox, (width-boundingBox-xSplit), xSplit do
		for y = boundingBox, height-boundingBox do
			for x = boundingBox, (width-boundingBox) do
--			for x = xOuter, xOuter+xSplit, 4 do
				var curKernelVal : double

				var curImOut : float = 0.0f
				var curVarOut : float = 0.0f
				var curMaskOut : uint16 = 0

				var curNorm : float = 0.0f

				var kernelArea = (boundingBox*2+1)*(boundingBox*2+1)
				var kernelWidth = boundingBox*2+1
--				kernelVals = [&float](C.malloc(kernelArea*5*sizeof(float)))


				escape 
					local luaKernelWidth = boundingBox*2+1

				    for j = -boundingBox, boundingBox do
				        for i = -boundingBox, boundingBox do
				        	emit quote
				--    for j = -boundingBox, boundingBox+1 do
				--        for i = -boundingBox, boundingBox+1 do
				--        		var curKernelLocation = (kernelWidth*(j+boundingBox) + (i+boundingBox))*5
				        		curKernelVal = polynomial1Double(x, y)*[kernelValSymbols6[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+1]] +
					                polynomial2Double(x, y)*[kernelValSymbols6[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+2]] + polynomial3Double(x, y)*[kernelValSymbols6[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+3]] + 
					                polynomial4Double(x, y)*[kernelValSymbols6[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+4]] + polynomial5Double(x, y)*[kernelValSymbols6[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+5]];

					            var curImInVec : double = [double](inputIm[(y+j)*width + x+i])
					            var curVarInVec : double = [double](inputVar[(y+j)*width + x+i])
					            var curMaskInVec = inputMask[(y+j)*width + x+i]


					            curImOut = curImOut + [float](curImInVec*curKernelVal); 
					            curVarOut = curVarOut + [float](curVarInVec*curKernelVal*curKernelVal);

					            if curKernelVal ~= 0.0f then
					            	curMaskOut = curMaskOut or curMaskInVec 
					            end

					            curNorm = curNorm + [float](curKernelVal);
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
	C.printf("\n\nMasked image lin combo, UNROLLING w/ SYMBOLS, Double precision kernel %dx%d blur Terra:\n", 2*boundingBox+1, 2*boundingBox+1)
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
--	C.free(kernelVals)
end



terra vectormask8(a : vector(float,8), b : vector(float,8))
    return terralib.select(a ~= b, [vector(uint16,8)](0xFFFFULL),[vector(uint16,8)](0) ) 
end

local kernelValSymbols1 = {}
local terra blurMaskedImageLinComboVectorize8UnrollSymbols()
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
			inputVar[y*width + x] = [float](y*x*C.cos(x/(y+1)))
			inputMask[y*width + x] = [uint16](x*C.cos(x*y) + y)


			outputIm[y*width + x] = 0.0f
			outputVar[y*width + x] = 0.0f
			outputMask[y*width + x] = 0
		end
	end

	var t1 = C.clock()

	escape

		local luaKernelArea = (boundingBox*2+1)*(boundingBox*2+1)
		local luaKernelWidth = boundingBox*2+1

        for i = 1, luaKernelArea*5 do
            local cur_kernelValue_symbol = symbol(float, "kVal"..i)
            table.insert(kernelValSymbols1, cur_kernelValue_symbol)
        end

		for j = -boundingBox, boundingBox do
			for i = -boundingBox, boundingBox do
				emit quote
					var [kernelValSymbols1[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+1]] = [kernel1(i, j)]
					var [kernelValSymbols1[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+2]] = [kernel2(i, j)]
					var [kernelValSymbols1[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+3]] = [kernel3(i, j)]
					var [kernelValSymbols1[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+4]] = [kernel4(i, j)]
					var [kernelValSymbols1[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+5]] = [kernel5(i, j)]
				end
			end
		end
	end


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
				local luaKernelWidth = boundingBox*2+1

			    for j = -boundingBox, boundingBox do
			        for i = -boundingBox, boundingBox do
			        	emit quote
			        		curKernelValTemp[0] = polynomial1(x, y)*[kernelValSymbols1[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+1]] +
				                polynomial2(x, y)*[kernelValSymbols1[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+2]]
				                 + polynomial3(x, y)*[kernelValSymbols1[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+3]]
				                 + polynomial4(x, y)*[kernelValSymbols1[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+4]]
				                 + polynomial5(x, y)*[kernelValSymbols1[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+5]];
			        		curKernelValTemp[1] = polynomial1(x+1, y)*[kernelValSymbols1[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+1]] +
				                polynomial2(x+1, y)*[kernelValSymbols1[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+2]]
				                 + polynomial3(x+1, y)*[kernelValSymbols1[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+3]]
				                 + polynomial4(x+1, y)*[kernelValSymbols1[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+4]]
				                 + polynomial5(x+1, y)*[kernelValSymbols1[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+5]];
			        		curKernelValTemp[2] = polynomial1(x+2, y)*[kernelValSymbols1[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+1]] +
				                polynomial2(x+2, y)*[kernelValSymbols1[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+2]]
				                 + polynomial3(x+2, y)*[kernelValSymbols1[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+3]]
				                 + polynomial4(x+2, y)*[kernelValSymbols1[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+4]]
				                 + polynomial5(x+2, y)*[kernelValSymbols1[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+5]];
			        		curKernelValTemp[3] = polynomial1(x+3, y)*[kernelValSymbols1[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+1]] +
				                polynomial2(x+3, y)*[kernelValSymbols1[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+2]]
				                 + polynomial3(x+3, y)*[kernelValSymbols1[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+3]]
				                 + polynomial4(x+3, y)*[kernelValSymbols1[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+4]]
				                 + polynomial5(x+3, y)*[kernelValSymbols1[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+5]];
			        		curKernelValTemp[4] = polynomial1(x+4, y)*[kernelValSymbols1[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+1]] +
				                polynomial2(x+4, y)*[kernelValSymbols1[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+2]]
				                 + polynomial3(x+4, y)*[kernelValSymbols1[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+3]]
				                 + polynomial4(x+4, y)*[kernelValSymbols1[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+4]]
				                 + polynomial5(x+4, y)*[kernelValSymbols1[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+5]];		
			        		curKernelValTemp[5] = polynomial1(x+5, y)*[kernelValSymbols1[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+1]] +
				                polynomial2(x+5, y)*[kernelValSymbols1[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+2]]
				                 + polynomial3(x+5, y)*[kernelValSymbols1[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+3]]
				                 + polynomial4(x+5, y)*[kernelValSymbols1[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+4]]
				                 + polynomial5(x+5, y)*[kernelValSymbols1[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+5]];
			        		curKernelValTemp[6] = polynomial1(x+6, y)*[kernelValSymbols1[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+1]] +
				                polynomial2(x+6, y)*[kernelValSymbols1[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+2]]
				                 + polynomial3(x+6, y)*[kernelValSymbols1[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+3]]
				                 + polynomial4(x+6, y)*[kernelValSymbols1[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+4]]
				                 + polynomial5(x+6, y)*[kernelValSymbols1[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+5]];
			        		curKernelValTemp[7] = polynomial1(x+7, y)*[kernelValSymbols1[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+1]] +
				                polynomial2(x+7, y)*[kernelValSymbols1[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+2]]
				                 + polynomial3(x+7, y)*[kernelValSymbols1[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+3]]
				                 + polynomial4(x+7, y)*[kernelValSymbols1[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+4]]
				                 + polynomial5(x+7, y)*[kernelValSymbols1[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+5]];
				            
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
	C.printf("\n\nVectorized 8wide masked image lin combo, unrolled w/Symbols %dx%d blur Terra:\n", 2*boundingBox+1, 2*boundingBox+1)
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




--local kernelValSymbols2 = {}
--local terra blurMaskedImageLinComboVectorizeAcrossKernelUnrollSymbols()
----	var boundingBox : int = 2
--	var width : int = 2048
--	var height : int = 1489
--
--	var inputIm : &float
--	inputIm = [&float](C.malloc(sizeof(float)*width*height))
--	var inputVar : &float
--	inputVar = [&float](C.malloc(sizeof(float)*width*height))
--	var inputMask : &uint16
--	inputMask = [&uint16](C.malloc(sizeof(uint16)*width*height))
--
--	var outputIm : &float
--	outputIm = [&float](C.malloc(sizeof(float)*width*height))
--	var outputVar : &float
--	outputVar = [&float](C.malloc(sizeof(float)*width*height))
--	var outputMask : &uint16
--	outputMask = [&uint16](C.malloc(sizeof(uint16)*width*height))
--
--	for x = 0, width do
--		for y = 0, height do
--			inputIm[y*width + x] = [float](y*x*C.cos(x/(y+1)))
--			inputVar[y*width + x] = [float](y*x*C.cos(x/(y+1)))
--			inputMask[y*width + x] = [uint16](x*C.cos(x*y) + y)
--
--
--			outputIm[y*width + x] = 0.0f
--			outputVar[y*width + x] = 0.0f
--			outputMask[y*width + x] = 0
--		end
--	end
--
--	var t1 = C.clock()
--
--	escape
--
--		--local luaKernelArea = 49
--		local luaKernelWidth = 7
--
--        for i = 1, (luaKernelWidth*5) do
--            local cur_kernelValue_symbol = symbol(vector(float,8), "kVal"..i)
--            table.insert(kernelValSymbols2, cur_kernelValue_symbol)
--        end
--
--		for j = -boundingBox, boundingBox do
--			emit quote
--				var [kernelValSymbols2[5*(j+boundingBox) + 1]] = vector([kernel1(-3, j)], [kernel1(-2, j)], [kernel1(-1, j)], [kernel1(0, j)], [kernel1(1, j)], [kernel1(2, j)], [kernel1(3, j)], 0.0f)
--				var [kernelValSymbols2[5*(j+boundingBox) + 2]] = vector([kernel2(-3, j)], [kernel2(-2, j)], [kernel2(-1, j)], [kernel2(0, j)], [kernel2(1, j)], [kernel2(2, j)], [kernel2(3, j)], 0.0f)
--				var [kernelValSymbols2[5*(j+boundingBox) + 3]] = vector([kernel3(-3, j)], [kernel3(-2, j)], [kernel3(-1, j)], [kernel3(0, j)], [kernel3(1, j)], [kernel3(2, j)], [kernel3(3, j)], 0.0f)
--				var [kernelValSymbols2[5*(j+boundingBox) + 4]] = vector([kernel4(-3, j)], [kernel4(-2, j)], [kernel4(-1, j)], [kernel4(0, j)], [kernel4(1, j)], [kernel4(2, j)], [kernel4(3, j)], 0.0f)
--				var [kernelValSymbols2[5*(j+boundingBox) + 5]] = vector([kernel5(-3, j)], [kernel5(-2, j)], [kernel5(-1, j)], [kernel5(0, j)], [kernel5(1, j)], [kernel5(2, j)], [kernel5(3, j)], 0.0f)
--
--			end
--		end
--	end
--
--
--	for y = boundingBox, height-boundingBox do
--		for x = boundingBox, (width-boundingBox) do
--			--var curOut : float = 0
--			--var curNorm : float = 0
--			var curKernelVal : vector(float,8)
--			var curKernelValTemp : float[8]
--			var curImOut : vector(float,8) = vector(0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f)
--			var curVarOut : vector(float,8) = vector(0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f)
--
----			var curMaskOut : vector(uint16,4) = vector(0,0,0,0)
--
--			var curMaskOut : vector(uint16, 8) = vector(0, 0, 0, 0,0, 0, 0, 0)
--			var zeroVec : vector(float, 8) = vector(0.f, 0.f, 0.f, 0.f,0.f,0.f,0.f,0.f)
--
--			var curNorm : vector(float,8) = vector(0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f)
----			var curImOutVec = @[&vector(float,4)](&outputIm[y*width + x])
--		
--			var p1 = polynomial1(x, y)
--			var p2 = polynomial2(x, y)
--			var p3 = polynomial3(x, y)
--			var p4 = polynomial4(x, y)
--			var p5 = polynomial5(x, y)
--
--
----        	curKernelVal = p1*[kernelValSymbols2[5*(0+boundingBox) + 1]] +
----	            p2*[kernelValSymbols2[5*(0+boundingBox) + 2]]
----	             + p3*[kernelValSymbols2[5*(0+boundingBox) + 3]]
----	             + p4*[kernelValSymbols2[5*(0+boundingBox) + 4]]
----	             + p5*[kernelValSymbols2[5*(0+boundingBox) + 5]];
----			curKernelVal = vector(1.3f, 4.5f, 3.5f, 42.23f,1.3f, 4.5f, 3.5f, 0.0f)
--
--			escape 
--				local luaKernelWidth = boundingBox*2+1
--			    for j = -boundingBox, boundingBox do
--		        	emit quote
--		        		curKernelVal = p1*[kernelValSymbols2[5*(j+boundingBox) + 1]] +
--			                p2*[kernelValSymbols2[5*(j+boundingBox) + 2]]
--			                 + p3*[kernelValSymbols2[5*(j+boundingBox) + 3]]
--			                 + p4*[kernelValSymbols2[5*(j+boundingBox) + 4]]
--			                 + p5*[kernelValSymbols2[5*(j+boundingBox) + 5]];
--
--			            var curImInVec = @[&vector(float,8)](&inputIm[(y+j)*width + x-3])
--			            var curVarInVec = @[&vector(float,8)](&inputVar[(y+j)*width + x-3])
--			            var curMaskInVec = @[&vector(uint16,8)](&inputMask[(y+j)*width + x-3])
--
--
--			            curImOut = curImOut + curImInVec*curKernelVal; 
--			            curVarOut = curVarOut + curVarInVec*curKernelVal*curKernelVal;
--
--			            var bitMask : vector(uint16, 8) = vectormask8(curKernelVal, zeroVec)
--			            curMaskOut = curMaskOut or (curMaskInVec and bitMask)
--
--			            curNorm = curNorm + curKernelVal;
--			        end
--			    end
--			end
--			var totalNorm = curNorm[0] + curNorm[1] + curNorm[2] + curNorm[3] + curNorm[4] + curNorm[5] + curNorm[6]
--		    var totalCurImOut = curImOut[0] + curImOut[1] + curImOut[2] + curImOut[3] + curImOut[4] + curImOut[5] + curImOut[6]
--		    var totalCurVarOut = curVarOut[0] + curVarOut[1] + curVarOut[2] + curVarOut[3] + curVarOut[4] + curVarOut[5] + curVarOut[6]
--		    var totalCurMaskOut = curMaskOut[0] or curMaskOut[1] or curMaskOut[2] or curMaskOut[3] or curMaskOut[4] or curMaskOut[5] or curMaskOut[6]
--
--		    totalCurImOut = totalCurImOut/totalNorm
--		    totalCurVarOut = totalCurVarOut/(totalNorm*totalNorm)
--
--		    outputIm[y*width + x] = totalCurImOut
--		    outputVar[y*width + x] = totalCurVarOut
--		    outputMask[y*width + x] = totalCurMaskOut
--		end
--	end
--
--	var t2 = C.clock()
--	C.printf("\n\nVectorized Across Kernel (7wide HARD) wide masked image lin combo, unrolled w/Symbols %dx%d blur Terra:\n", 2*boundingBox+1, 2*boundingBox+1)
--	C.printf("outputIm[boundingBox*width + boundingBox] = %f, computation took: %f ms\n", outputIm[boundingBox*width + boundingBox],  (t2-t1)/1000.0)
--	C.printf("C.CLOCKS_PER_SEC = %d\n", C.CLOCKS_PER_SEC)
--
--	C.printf("Image plane, 10x10 box begining at (boundingBox,boundingBox)\n")
--	for i=boundingBox,boundingBox+10 do
--		for j=boundingBox,boundingBox+10 do
--			C.printf("%f\t", outputIm[i*width + j])
--		end
--		C.printf("\n")
--	end
--
--	C.printf("Variance plane, 10x10 box begining at (boundingBox,boundingBox)\n")
--	for i=boundingBox,boundingBox+10 do
--		for j=boundingBox,boundingBox+10 do
--			C.printf("%f\t", outputVar[i*width + j])
--		end
--		C.printf("\n")
--	end
--
--	C.printf("Mask plane, 10x10 box begining at (boundingBox,boundingBox)\n")
--	for i=boundingBox,boundingBox+10 do
--		for j=boundingBox,boundingBox+10 do
--			C.printf("%d\t", outputMask[i*width + j])
--		end
--		C.printf("\n")
--	end
--
--    C.free(inputIm)
--	C.free(outputIm)
--	C.free(inputVar)
--	C.free(outputVar)
--	C.free(inputMask)
--	C.free(outputMask)
--end





--deal with image, mask, and variance planes
local kernelValSymbols0 = {}
local terra blurMaskedImageLinComboUnrollSymbolsSplitX()
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

	var kernelVals : &float

	for x = 0, width do
		for y = 0, height do
			inputIm[y*width + x] = [float](y*x*C.cos(x/(y+1)))
			inputVar[y*width + x] = [float](y*x*C.cos(x/(y+1)))
			inputMask[y*width + x] = [uint16](x*C.cos(x*y) + y)


			outputIm[y*width + x] = 0.0f
			outputVar[y*width + x] = 0.0f
			outputMask[y*width + x] = 0
		end
	end

	var t1 = C.clock()

	var xSplit : int = 50
	for xOuter = boundingBox, (width-boundingBox-xSplit), xSplit do
		for y = boundingBox, height-boundingBox do
			for x = xOuter, xOuter+xSplit do
				var curKernelVal : float

				var curImOut : float = 0.0f
				var curVarOut : float = 0.0f
				var curMaskOut : uint16 = 0

				var curNorm : float = 0.0f

				var kernelArea = (boundingBox*2+1)*(boundingBox*2+1)
				var kernelWidth = boundingBox*2+1
				kernelVals = [&float](C.malloc(kernelArea*5*sizeof(float)))



				escape

					local luaKernelArea = (boundingBox*2+1)*(boundingBox*2+1)
					local luaKernelWidth = boundingBox*2+1

			        for i = 1, luaKernelArea*5 do
			            local cur_kernelValue_symbol = symbol(float, "kVal"..i)
			            table.insert(kernelValSymbols0, cur_kernelValue_symbol)
			        end

					for j = -boundingBox, boundingBox do
						for i = -boundingBox, boundingBox do
							emit quote
								var [kernelValSymbols0[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+1]] = [kernel1(i, j)]
								var [kernelValSymbols0[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+2]] = [kernel2(i, j)]
								var [kernelValSymbols0[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+3]] = [kernel3(i, j)]
								var [kernelValSymbols0[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+4]] = [kernel4(i, j)]
								var [kernelValSymbols0[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+5]] = [kernel5(i, j)]
							end
						end
					end
				end

				escape 
					local luaKernelWidth = boundingBox*2+1

				    for j = -boundingBox, boundingBox do
				        for i = -boundingBox, boundingBox do
				        	emit quote
				--    for j = -boundingBox, boundingBox+1 do
				--        for i = -boundingBox, boundingBox+1 do
				--        		var curKernelLocation = (kernelWidth*(j+boundingBox) + (i+boundingBox))*5
				        		curKernelVal = polynomial1(x, y)*[kernelValSymbols0[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+1]] +
					                polynomial2(x, y)*[kernelValSymbols0[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+2]] + polynomial3(x, y)*[kernelValSymbols0[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+3]] + 
					                polynomial4(x, y)*[kernelValSymbols0[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+4]] + polynomial5(x, y)*[kernelValSymbols0[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+5]];

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
	end

	var t2 = C.clock()
	C.printf("\n\nMasked image lin combo, UNROLLING w/ SYMBOLS, x split by %d, %dx%d blur Terra:\n", xSplit, 2*boundingBox+1, 2*boundingBox+1)
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
	C.free(kernelVals)
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

	for x = 0, width do
		for y = 0, height do
			inputIm[y*width + x] = [float](y*x*C.cos(x/(y+1)))
			inputVar[y*width + x] = [float](y*x*C.cos(x/(y+1)))
			inputMask[y*width + x] = [uint16](x*C.cos(x*y) + y)


			outputIm[y*width + x] = 0.0f
			outputVar[y*width + x] = 0.0f
			outputMask[y*width + x] = 0
		end
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
	C.printf("\n\nVectorized masked image lin combo %dx%d blur Terra:\n", 2*boundingBox+1, 2*boundingBox+1)
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

	for x = 0, width do
		for y = 0, height do
			inputIm[y*width + x] = [float](y*x*C.cos(x/(y+1)))
			inputVar[y*width + x] = [float](y*x*C.cos(x/(y+1)))
			inputMask[y*width + x] = [uint16](x*C.cos(x*y) + y)


			outputIm[y*width + x] = 0.0f
			outputVar[y*width + x] = 0.0f
			outputMask[y*width + x] = 0
		end
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
	C.printf("\n\nVectorized 8wide masked image lin combo %dx%d blur Terra:\n", 2*boundingBox+1, 2*boundingBox+1)
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
	C.printf("\n\nVectorized masked image lin combo, ORing all mask pixels %dx%d blur Terra:\n", 2*boundingBox+1, 2*boundingBox+1)
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


local kernelValSymbols3 = {}
local terra blurMaskedImageLinComboVectorize8UnrollSymbolsInterp()
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
			inputVar[y*width + x] = [float](y*x*C.cos(x/(y+1)))
			inputMask[y*width + x] = [uint16](x*C.cos(x*y) + y)


			outputIm[y*width + x] = 0.0f
			outputVar[y*width + x] = 0.0f
			outputMask[y*width + x] = 0
		end
	end

	var t1 = C.clock()

	escape

		local luaKernelArea = (boundingBox*2+1)*(boundingBox*2+1)
		local luaKernelWidth = boundingBox*2+1

        for i = 1, luaKernelArea*5 do
            local cur_kernelValue_symbol = symbol(float, "kVal"..i)
            table.insert(kernelValSymbols3, cur_kernelValue_symbol)
        end

		for j = -boundingBox, boundingBox do
			for i = -boundingBox, boundingBox do
				emit quote
					var [kernelValSymbols3[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+1]] = [kernel1(i, j)]
					var [kernelValSymbols3[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+2]] = [kernel2(i, j)]
					var [kernelValSymbols3[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+3]] = [kernel3(i, j)]
					var [kernelValSymbols3[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+4]] = [kernel4(i, j)]
					var [kernelValSymbols3[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+5]] = [kernel5(i, j)]
				end
			end
		end
	end

	for yOuter = boundingBox, (height-boundingBox)-interpDist, interpDist do
		for x = boundingBox, (width-boundingBox), 8 do
			--var curOut : float = 0
			--var curNorm : float = 0
			var curKernelVal : vector(float,8)
			var curKernelVal_yOuter : vector(float,8)
			var curKernelVal_plusInterp : vector(float,8)
			var curKernelVal_dy : vector(float,8) --change in the y direction

			var curKernelValTemp : float[8]

			var curImOut : vector(float,8)[interpDist]
			var curVarOut : vector(float,8)[interpDist]
			var curMaskOut : vector(uint16, 8)[interpDist]
			var curNorm : vector(float,8)[interpDist]


			for yInner = 0, interpDist do
				curImOut[yInner] = vector(0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f)
				curVarOut[yInner] = vector(0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f)
				curMaskOut[yInner] = vector(0, 0, 0, 0,0, 0, 0, 0)
				curNorm[yInner] = vector(0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f)
			end


			var zeroVec : vector(float, 8) = vector(0.f, 0.f, 0.f, 0.f,0.f,0.f,0.f,0.f)

			escape 
				local luaKernelWidth = boundingBox*2+1

			    for j = -boundingBox, boundingBox do
			        for i = -boundingBox, boundingBox do
			        	for k = 0, 7 do
			        		emit quote
				        		curKernelValTemp[k] = polynomial1(x+k, yOuter)*[kernelValSymbols3[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+1]] +
					                polynomial2(x+k, yOuter)*[kernelValSymbols3[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+2]]
					                 + polynomial3(x+k, yOuter)*[kernelValSymbols3[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+3]]
					                 + polynomial4(x+k, yOuter)*[kernelValSymbols3[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+4]]
					                 + polynomial5(x+k, yOuter)*[kernelValSymbols3[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+5]];
			        		end
			        	end
			        	emit quote
			        		curKernelVal_yOuter = @[&vector(float,8)](&curKernelValTemp[0])
			        	end
			        	for k = 0, 7 do
			        		emit quote
				        		curKernelValTemp[k] = polynomial1(x+k, yOuter+(interpDist - 1))*[kernelValSymbols3[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+1]] +
					                polynomial2(x+k, yOuter+(interpDist - 1))*[kernelValSymbols3[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+2]]
					                 + polynomial3(x+k, yOuter+(interpDist - 1))*[kernelValSymbols3[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+3]]
					                 + polynomial4(x+k, yOuter+(interpDist - 1))*[kernelValSymbols3[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+4]]
					                 + polynomial5(x+k, yOuter+(interpDist - 1))*[kernelValSymbols3[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+5]];
			        		end
			        	end
			        	emit quote
			        		curKernelVal_plusInterp = @[&vector(float,8)](&curKernelValTemp[0])
			        		curKernelVal_dy = (curKernelVal_plusInterp - curKernelVal_yOuter)/(interpDist - 1)
			        	end
			        		
			        	for yInner = 0, interpDist - 1 do
				        	emit quote							            
					            curKernelVal = curKernelVal_yOuter + curKernelVal_dy*yInner

					            var curImInVec = @[&vector(float,8)](&inputIm[(yOuter+yInner+j)*width + x+i])
					            var curVarInVec = @[&vector(float,8)](&inputVar[(yOuter+yInner+j)*width + x+i])
					            var curMaskInVec = @[&vector(uint16,8)](&inputMask[(yOuter+yInner+j)*width + x+i])


					            curImOut[yInner] = curImOut[yInner] + curImInVec*curKernelVal; 
					            curVarOut[yInner] = curVarOut[yInner] + curVarInVec*curKernelVal*curKernelVal;

					            var bitMask : vector(uint16, 8) = vectormask8(curKernelVal, zeroVec)
					            curMaskOut[yInner] = curMaskOut[yInner] or (curMaskInVec and bitMask)

					            curNorm[yInner] = curNorm[yInner] + curKernelVal;
					        end
					    end
			        end
			    end
			end

			for yInner = 0, interpDist do
				curImOut[yInner] = curImOut[yInner]/curNorm[yInner]
		    	curVarOut[yInner] = curVarOut[yInner]/(curNorm[yInner]*curNorm[yInner])

				for vecOffset = 0, 8 do
		    		outputIm[(yOuter+yInner)*width + x+vecOffset] = curImOut[yInner][vecOffset]
		    		outputVar[(yOuter+yInner)*width + x+vecOffset] = curVarOut[yInner][vecOffset]
		    		outputMask[(yOuter+yInner)*width + x+vecOffset] = curMaskOut[yInner][vecOffset]

				end
			end
		end
	end

	var t2 = C.clock()
	C.printf("\n\nVectorized 8wide masked image lin combo, unrolled w/Symbols, Interpolation Experiment%dx%d blur Terra:\n", 2*boundingBox+1, 2*boundingBox+1)
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



--deal with image, mask, and variance planes
local kernelValSymbols4 = {}
local curKernelValSymbols = {}
local terra blurMaskedImageLinComboUnrollSymbolsInterpolate()
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

	var kernelVals : &float

	for x = 0, width do
		for y = 0, height do
			inputIm[y*width + x] = [float](y*x*C.cos(x/(y+1)))
			inputVar[y*width + x] = [float](y*x*C.cos(x/(y+1)))
			inputMask[y*width + x] = [uint16](x*C.cos(x*y) + y)


			outputIm[y*width + x] = 0.0f
			outputVar[y*width + x] = 0.0f
			outputMask[y*width + x] = 0
		end
	end

	var t1 = C.clock()

	escape

		local luaKernelArea = (boundingBox*2+1)*(boundingBox*2+1)
		local luaKernelWidth = boundingBox*2+1

        for i = 1, luaKernelArea*5 do
            local cur_kernelValue_symbol = symbol(float, "kVal"..i)
            table.insert(kernelValSymbols4, cur_kernelValue_symbol)
        end


		for j = -boundingBox, boundingBox do
			for i = -boundingBox, boundingBox do
				emit quote
					var [kernelValSymbols4[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+1]] = [kernel1(i, j)]
					var [kernelValSymbols4[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+2]] = [kernel2(i, j)]
					var [kernelValSymbols4[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+3]] = [kernel3(i, j)]
					var [kernelValSymbols4[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+4]] = [kernel4(i, j)]
					var [kernelValSymbols4[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+5]] = [kernel5(i, j)]
				end
			end
		end
	end

	var tB = C.clock()
	C.printf("time spent calculating kernel = %f ms!!!!!!!!!!!!!!@@@@@@@@@@@@", (float)(tB-t1)/C.CLOCKS_PER_SEC*1000.0f)



--	var xSplit : int = 32
--	for xOuter = boundingBox, (width-boundingBox-xSplit), xSplit do
		for yOuter = boundingBox, height-boundingBox - interpDist, interpDist do
			for xOuter = boundingBox, (width-boundingBox) - interpDist, interpDist do

				escape		
					local luaKernelArea = (boundingBox*2+1)*(boundingBox*2+1)
					local luaKernelWidth = boundingBox*2+1

			        for i = 1, luaKernelArea do
			            local cur_kernelValue_symbol = symbol(float, "kVal"..i)
			            table.insert(curKernelValSymbols, cur_kernelValue_symbol)
			        end

					for j = -boundingBox, boundingBox do
						for i = -boundingBox, boundingBox do
							emit quote
								var [curKernelValSymbols[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))+1]] = 
										polynomial1(xOuter, yOuter)*[kernelValSymbols4[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+1]] +
						                polynomial2(xOuter, yOuter)*[kernelValSymbols4[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+2]] +
						                polynomial3(xOuter, yOuter)*[kernelValSymbols4[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+3]] + 
						                polynomial4(xOuter, yOuter)*[kernelValSymbols4[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+4]] +
						                polynomial5(xOuter, yOuter)*[kernelValSymbols4[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+5]]

							end
						end
					end
				end

				for yInner = 0, interpDist do
					for xInner = 0, interpDist do
						var curKernelVal : float
						var curImOut : float = 0.0f
						var curVarOut : float = 0.0f
						var curMaskOut : uint16 = 0

						var curNorm : float = 0.0f

						var kernelArea = (boundingBox*2+1)*(boundingBox*2+1)
						var kernelWidth = boundingBox*2+1
						kernelVals = [&float](C.malloc(kernelArea*5*sizeof(float)))




						escape 
							local luaKernelWidth = boundingBox*2+1

						    for j = -boundingBox, boundingBox do
						        for i = -boundingBox, boundingBox do
						        	emit quote
						        		curKernelVal = [curKernelValSymbols[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))+1]]

							            var curImInVec = inputIm[((yOuter+yInner)+j)*width + (xOuter+xInner)+i]
							            var curVarInVec = inputVar[((yOuter+yInner)+j)*width + (xOuter+xInner)+i]
							            var curMaskInVec = inputMask[((yOuter+yInner)+j)*width + (xOuter+xInner)+i]


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
					    outputIm[(yOuter+yInner)*width + (xOuter+xInner)] = curImOut

					    curVarOut = curVarOut/(curNorm*curNorm)
					    outputVar[(yOuter+yInner)*width + (xOuter+xInner)] = curVarOut

					    outputMask[(yOuter+yInner)*width + (xOuter+xInner)] = curMaskOut
					end
				end
			end
		end
	--end

	var t2 = C.clock()
	C.printf("\n\nMasked image lin combo, UNROLLING w/ SYMBOLS, InterpolationExperiment no explicit Vectorization %dx%d blur Terra:\n", 2*boundingBox+1, 2*boundingBox+1)
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
	C.free(kernelVals)
end



--deal with image, mask, and variance planes
local kernelValSymbols5 = {}

local terra blurMaskedImageLinComboUnrollSymbolsTile()
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

	var kernelVals : &float

	for x = 0, width do
		for y = 0, height do
			inputIm[y*width + x] = [float](y*x*C.cos(x/(y+1)))
			inputVar[y*width + x] = [float](y*x*C.cos(x/(y+1)))
			inputMask[y*width + x] = [uint16](x*C.cos(x*y) + y)


			outputIm[y*width + x] = 0.0f
			outputVar[y*width + x] = 0.0f
			outputMask[y*width + x] = 0
		end
	end

	var t1 = C.clock()

	escape

		local luaKernelArea = (boundingBox*2+1)*(boundingBox*2+1)
		local luaKernelWidth = boundingBox*2+1

        for i = 1, luaKernelArea*5 do
            local cur_kernelValue_symbol = symbol(float, "kVal"..i)
            table.insert(kernelValSymbols5, cur_kernelValue_symbol)
        end


		for j = -boundingBox, boundingBox do
			for i = -boundingBox, boundingBox do
				emit quote
					var [kernelValSymbols5[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+1]] = [kernel1(i, j)]
					var [kernelValSymbols5[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+2]] = [kernel2(i, j)]
					var [kernelValSymbols5[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+3]] = [kernel3(i, j)]
					var [kernelValSymbols5[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+4]] = [kernel4(i, j)]
					var [kernelValSymbols5[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+5]] = [kernel5(i, j)]
				end
			end
		end
	end

	var tB = C.clock()
	C.printf("time spent calculating kernel = %f ms!!!!!!!!!!!!!!@@@@@@@@@@@@", (float)(tB-t1)/C.CLOCKS_PER_SEC*1000.0f)



--	var xSplit : int = 32
--	for xOuter = boundingBox, (width-boundingBox-xSplit), xSplit do
		for yOuter = boundingBox, height-boundingBox - interpDist, interpDist do
			for xOuter = boundingBox, (width-boundingBox) - interpDist, interpDist do
				for yInner = 0, interpDist do
					for xInner = 0, interpDist do
						var y = yOuter + yInner
						var x = xOuter + xInner
		--			for x = xOuter, xOuter+xSplit, 4 do
						var curKernelVal : float

						var curImOut : float = 0.0f
						var curVarOut : float = 0.0f
						var curMaskOut : uint16 = 0

						var curNorm : float = 0.0f

						var kernelArea = (boundingBox*2+1)*(boundingBox*2+1)
						var kernelWidth = boundingBox*2+1
		--				kernelVals = [&float](C.malloc(kernelArea*5*sizeof(float)))


						escape 
							local luaKernelWidth = boundingBox*2+1

						    for j = -boundingBox, boundingBox do
						        for i = -boundingBox, boundingBox do
						        	emit quote
						--    for j = -boundingBox, boundingBox+1 do
						--        for i = -boundingBox, boundingBox+1 do
						--        		var curKernelLocation = (kernelWidth*(j+boundingBox) + (i+boundingBox))*5
						        		curKernelVal = polynomial1(x, y)*[kernelValSymbols5[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+1]] +
							                polynomial2(x, y)*[kernelValSymbols5[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+2]] + polynomial3(x, y)*[kernelValSymbols5[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+3]] + 
							                polynomial4(x, y)*[kernelValSymbols5[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+4]] + polynomial5(x, y)*[kernelValSymbols5[(luaKernelWidth*(j+boundingBox) + (i+boundingBox))*5+5]];

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
			end
		end
	--end

	var t2 = C.clock()
	C.printf("\n\nMasked image lin combo, UNROLLING w/ SYMBOLS, TILED %dx%d blur Terra:\n", 2*boundingBox+1, 2*boundingBox+1)
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
--	C.free(kernelVals)
end


--the first term is the name in the .o and .h files, the second name is the name in this file
--terralib.saveobj("blurExampleTerra.o",{ funcNameInC = blurImageLinComboVectorize4 })

--blurArrayTerra()
--blurArrayTerraVectorizedAcrossX()

--blurImageLinCombo(true)
--blurImageLinCombo(false)
--blurImageLinCombo("testingRefactor")
--blurImageLinCombo("normal")
--blurImageLinComboVectorize2()
--blurImageLinCombo(false)
--blurImageLinComboVectorize2()
--blurImageLinComboVectorize4()
--blurImageLinComboVectorize8()
--blurImageLinComboVectorize8Wide1()
--blurImageLinComboVectorize16Wide()


--@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

print("luaKernel1:")
for j = -boundingBox, boundingBox do
	for i = -boundingBox, boundingBox do
		io.write(luaKernel1(i, j), "\t")
	end
	io.write("\n")
end

local terra printKernel1()
	C.printf("terraKernel1:\n")
	for j = -boundingBox, boundingBox + 1 do
		for i = -boundingBox, boundingBox + 1 do
			C.printf("%E\t", kernel1(i, j))
		end
		C.printf("\n")
	end 
end

printKernel1()

--blurImageLinComboSplitX()

--blurImageLinComboNoUnroll()

--blurMaskedImageLinComboNoUnroll()
--
--!!!!blurMaskPlaneLinComboTest()
--!!!!blurMaskPlaneLinCombo()
--!!!!blurVarianceLinCombo()
--
--blurMaskedImageLinComboUnrollSymbolsSplitX()

--blurMaskedImageLinComboNoUnroll()

--blurMaskedImageLinComboUnrollSymbolsTile()

--blurMaskedImageLinComboVectorize8UnrollSymbolsInterp()
--blurMaskedImageLinComboUnrollSymbolsInterpolate()
--blurMaskedImageLinComboVectorizeAcrossKernelUnrollSymbols()
luaBlurMaskedImageLinComboUnrollSymbols(5)


blurMaskedImageLinComboUnrollSymbolsDoublePrecisionKernel()

blurMaskedImageLinComboVectorize8UnrollSymbols()

--
blurMaskedImageLinCombo()
--blurMaskedImageLinComboSplitX()

--blurMaskedImageLinComboVectorize()
blurMaskedImageLinComboVectorize8()



--@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
----blurMaskedImageOrAllMaskLinComboVectorize()