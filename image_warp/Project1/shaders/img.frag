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
        for (int i = 1-kernel_half; i <= kernel_half; ++i) {
            for (int j = 1-kernel_half; j <=kernel_half; ++j) {
                vec2 coord = texcoord + vec2(1.0*i/W, 1.0*j/H);
                intensity += kernelX[31+i] * kernelY[31+j] * texture(imgTex, coord).r;
                var += kernelX[31+i]*kernelX[31+i]*kernelY[31+j]*kernelY[31+j] * texture(varTex, coord).r;
    			if (kernelX[31+i]*kernelY[31+j] > 0.0) {
                    mask |= texture(maskTex, coord).r;
                }
            }
        }
        outImg = intensity * area;
        outVar = var * area;
        outMask = uvec4(mask, 0u, 0u, 0u);
	}
}