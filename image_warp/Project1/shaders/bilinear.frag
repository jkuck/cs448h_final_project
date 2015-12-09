# version 330
#define PI 3.1415926535897932384626433832795
in vec2 texCoord;

layout(location=0) out float outImg;
layout(location=1) out float outVar;
layout(location=2) out uvec4 outMask;

uniform float W;
uniform float H;
uniform int kernel_half;

uniform sampler2D imgTex;
uniform sampler2D varTex;
uniform usampler2D maskTex;

float lanczos(float x) {
    if (abs(x) < 0.001) {
        return 1;
    } else if (abs(x) > float(kernel_half)) {
        return 0;
    }
    return 1.0*kernel_half * sin(PI*x) * sin(PI*x/kernel_half) / (PI*PI*x*x);
}

void computeKernel(float pos, float stepSize,inout float[64] kernel) {
    for (int i = 1-kernel_half; i <= kernel_half; ++i) {
        kernel[31+i] = lanczos(pos/stepSize - (floor(pos/stepSize)+i+0.5));
    }
}

float computeAreaNorm(vec2 texcoord) {
    vec2 xDeriv = dFdx(texcoord);
    vec2 yDeriv = dFdy(texcoord);
    float area = xDeriv.x * yDeriv.y - xDeriv.y * yDeriv.x;
    if (area < 0) area = - area;
    return W * H * area;
}

void main() {
    vec2 texcoord = texCoord;
    texcoord[1] = 1.0 - texCoord[1];
	if (texcoord[0] < 0.0 || texcoord[1] < 0.0 || texcoord[0] > 1.0 || texcoord[1] < 0.0) {
        outImg = 0.0;
        outVar = 0.0;
        outMask = uvec4(0xff, 0u, 0u, 0u);
	} else {
        float area = computeAreaNorm(texcoord);
        float[64] kernelX;
        float[64] kernelY;
        computeKernel(texcoord[0], 1/W, kernelX);
        computeKernel(texcoord[1], 1/H, kernelY);
        float intensity = 0.0;
        float var = 0.0;
        uint mask = uint(0u);
        for (int i = 1-kernel_half; i < kernel_half; i += 2) {
            for (int j = 1-kernel_half; j < kernel_half; j += 2) {
			    float coef = 0.0;
				vec2 coord = vec2(0.0, 0.0);
				vec2 coord_00 =  vec2(floor(texcoord[0]*W+i)/W, floor(texcoord[1]*H+j)/H);
				for (int x = 0; x < 2; ++x) {
				    for (int y = 0; y < 2; ++y) {
					   coef += kernelX[31+i+x] * kernelY[31+j+y];
					   coord += kernelX[31+i+x] * kernelY[31+j+y] * (coord_00 + vec2(1.0*x/W, 1.0*y/H));
					}
				}
				coord /= coef;
                intensity += coef * texture(imgTex, coord).r;
                var += coef * coef * texture(varTex, coord).r;
    			if (coef > 0.0) {
                    mask |= texture(maskTex, coord).r;
                }
            }
        }
        outImg = intensity * area;
        outVar = var * area;
        outMask = uvec4(mask, 0u, 0u, 0u);
	}
}