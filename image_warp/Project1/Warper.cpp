#include <assert.h>
#include <float.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <sstream>
#include <stdio.h>
#include <time.h>

#include <GL\glew.h>
#include <GL\GL.h>
#include <GL\GLU.h>
#include "Warper.h"

Warper::Warper(WarpData input) : input(input), warped(false) {
    if (input.fetchMethod == DIRECT) {
        vShaderFile = "shaders/img.vert";
        fShaderFile = "shaders/img.frag";
    } else if (input.fetchMethod == BILINEAR) {
        vShaderFile = "shaders/bilinear.vert";
        fShaderFile = "shaders/bilinear.frag";
    }
    loadInputs();
    constructArrays();
}

Warper::~Warper() {
    for (int i = 0; i < 2; ++i) {
        free(textures[i]);
        free(targets[i]);
    }
    free(maskTexture);
    free(maskTarget);
    free(templateCoords);
    free(scienceCoords);
}

float Warper::calc_PSNR(string refOutputFile) {
    readInput((void**)&targetReference, refOutputFile);
    float val = psnr(targetReference, targets[0]);
    free(targetReference);
    return val;
}

float Warper::psnr(float* ref, float* result) {
	float resultMax = FLT_MIN;
	float refMax = FLT_MIN;
	float resultMin = FLT_MAX;
	float refMin = FLT_MAX;
	for (int i = 0; i < input.w*input.h; ++i) {
		if (ref[i] > refMax) refMax = ref[i];
		if (ref[i] < refMin) refMin = ref[i];
		if (result[i] > resultMax) resultMax = result[i];
		if (result[i] < resultMin) resultMin = result[i];
	}
	float mse = 0.0;
	int nanCount = 0;
	for (int i = 0; i < input.w*input.h; ++i) {
		float sol = (ref[i] - refMin) / (refMax - refMin);
		float noise = (result[i] - resultMin) / (resultMax - resultMin);
		float diff = noise - sol;
		float msePlus = diff * diff;
		if (!(msePlus == msePlus)) {
			nanCount++;
		} else {
			mse += diff * diff;
		}
	}
	mse /= (input.w * input.h - nanCount);
	return -10 * log(mse);
}

void Warper::warp(ofstream& fs) {
    // glew and context setup
    GLFWwindow* window = initGLFW();
    glfwMakeContextCurrent(window);
    glewExperimental = GL_TRUE;
    glewInit();
    // bind buffer objects
    bindBuffers();
    // load shader
    createShaderProgram();
    // va array
    setupVertexAttribArray();
    // textures
    bindTextures();
    setUniforms();
    setupFramebuffer();
	clock_t t = clock();

    glBindFramebuffer(GL_FRAMEBUFFER, fbo);
    GLenum buf[3] = { GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1, GL_COLOR_ATTACHMENT2 };
    glDrawBuffers(3, buf);
    glViewport(0, 0, input.w, input.h);
    glClearColor(0.0f, 0.0f, 0.0f, 1.f);
    glClear(GL_COLOR_BUFFER_BIT);
    glDrawElements(GL_TRIANGLES, 6 * input.grid_w*input.grid_h, GL_UNSIGNED_INT, 0);
	glFinish();
	t = clock() - t;
	fs << (double)t / CLOCKS_PER_SEC << ",";
	storeWarpedResults();
	cleanup();

    warped = true;
}

GLFWwindow* Warper::initGLFW() {
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_COMPAT_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwWindowHint(GLFW_RESIZABLE, GL_TRUE);

    return glfwCreateWindow(input.w, input.h, "Warper", nullptr, nullptr);
}

// modified from cplusplus.com
void Warper::readInput(void **buffer, const string& filename) {
    FILE *pFile;
    fopen_s(&pFile, filename.c_str(), "rb");
    if (pFile==NULL) {fputs ("File error",stderr); exit (1);}

    // obtain file size:
    fseek (pFile , 0 , SEEK_END);
    long lSize = ftell (pFile);
    rewind (pFile);

    // allocate memory to contain the whole file:
    *buffer = (char*) malloc (sizeof(char)*lSize);
    if (*buffer == NULL) {fputs ("Memory error",stderr); exit (2);}

    // copy the file into the buffer:
    size_t result = fread (*buffer,1,lSize,pFile);
    if (result != lSize) {fputs ("Reading error",stderr); exit (3);}
    fclose (pFile);
}

void Warper::loadInputs() {
    readInput((void**)&textures[0], input.inputFITS[0]);
    readInput((void**)&textures[1], input.inputFITS[1]);
    readInput((void**)&maskTexture, input.inputFITS[2]);
    readInput((void**)&templateCoords, input.inputCoords[0]);
    readInput((void**)&scienceCoords, input.inputCoords[1]);

	coords = (double*) malloc(4 * sizeof(double) * input.grid_w * input.grid_h);
	gridElements = (GLuint*) malloc(6 * sizeof(GLuint) * input.grid_w * input.grid_h);

    targets[0] = (float*)malloc(sizeof(float)*input.w*input.h);
    targets[1] = (float*)malloc(sizeof(float)*input.w*input.h);
    maskTarget = (unsigned short*)malloc(sizeof(unsigned short)*input.w*input.h);
}

void Warper::constructArrays() {
    for (int i = 0; i < input.grid_w; ++i) {
        for (int j = 0; j < input.grid_h; ++j) {
            double sourceX = templateCoords[2 * (i + j*input.grid_w)] / input.w;
            double sourceY = 1 - templateCoords[2 * (i + j*input.grid_w) + 1] / input.h;
            double destX = scienceCoords[2 * (i + j*input.grid_w)] / (input.w/2) - 1.0;
            double destY = scienceCoords[2 * (i + j*input.grid_w) + 1] / (input.h/2) - 1.0;
            coords[4 * (i + j*input.grid_w)] = destX;
            coords[4 * (i + j*input.grid_w)+1] = destY;
            coords[4 * (i + j*input.grid_w)+2] = sourceX;
            coords[4 * (i + j*input.grid_w)+3] = sourceY;
            if (i < input.grid_w - 1 && j < input.grid_h - 1) {
                gridElements[6 * (i + j*input.grid_w)] = (i + j*input.grid_w);
                gridElements[6 * (i + j*input.grid_w) + 1] = (i + 1 + j*input.grid_w);
                gridElements[6 * (i + j*input.grid_w) + 2] = (i + 1 + (j + 1)*input.grid_w);
                gridElements[6 * (i + j*input.grid_w) + 3] = (i + j*input.grid_w);
                gridElements[6 * (i + j*input.grid_w) + 4] = (i + 1 + (j + 1)*input.grid_w);
                gridElements[6 * (i + j*input.grid_w) + 5] = (i + (j + 1)*input.grid_w);
            }
        }
    }
}

void Warper::bindBuffers() {
    // vao
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);
    // vbo
    glGenBuffers(1, &vbo);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, 4*sizeof(double)*input.grid_w*input.grid_h, coords, GL_STATIC_DRAW);
    //EBO
    glGenBuffers(1, &ebo);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, 6*sizeof(GLuint)*input.grid_w*input.grid_h, gridElements, GL_STATIC_DRAW);
}

void Warper::loadShader(const string& file, GLenum vOrF) {
    // load shader file
    ifstream in(file.c_str());
    if (!in) {
        fprintf(stderr, "Failed to open shader: %s\n", file);
        assert(false);
    }
    ostringstream oss;
    oss << in.rdbuf();
    string fileStr = oss.str();
    // load shader to GPU
    const char *source = fileStr.c_str();
    GLint status = GL_FALSE;
    char buffer[512];
    if (vOrF == GL_VERTEX_SHADER) {
        vertShader = glCreateShader(GL_VERTEX_SHADER);
        glShaderSource(vertShader, 1, &source, NULL);
        glCompileShader(vertShader);
        glGetShaderiv(vertShader, GL_COMPILE_STATUS, &status);
        glGetShaderInfoLog(vertShader, 512, NULL, buffer);
    } else {
        fragShader = glCreateShader(GL_FRAGMENT_SHADER);
        glShaderSource(fragShader, 1, &source, NULL);
        glCompileShader(fragShader);
        glGetShaderiv(fragShader, GL_COMPILE_STATUS, &status);
        glGetShaderInfoLog(fragShader, 512, NULL, buffer);
    }
    if (status != GL_TRUE) {
        fprintf(stderr, "Compiling shader failed: %s\n", file.c_str());
        fprintf(stderr, "Compiling shader failed: %s\n", buffer);
        assert(false);
    }
}

void Warper::createShaderProgram() {
    shaderProgram = glCreateProgram();
    loadShader(vShaderFile, GL_VERTEX_SHADER);
    loadShader(fShaderFile, GL_FRAGMENT_SHADER);
    glAttachShader(shaderProgram, vertShader);
    glAttachShader(shaderProgram, fragShader);
    glLinkProgram(shaderProgram);
    glUseProgram(shaderProgram);
}

void Warper::setUniforms() {
    GLuint texWUni = glGetUniformLocation(shaderProgram, "W");
    glUniform1f(texWUni, input.w);
    GLuint texHUni = glGetUniformLocation(shaderProgram, "H");
    glUniform1f(texHUni, input.h);
    GLuint kernelSizeUni = glGetUniformLocation(shaderProgram, "kernel_half");
    glUniform1i(kernelSizeUni, input.kernelSize/2);
    GLuint texLoc = glGetUniformLocation(shaderProgram, "imgTex");
    glUniform1i(texLoc, 0);
    texLoc = glGetUniformLocation(shaderProgram, "varTex");
    glUniform1i(texLoc, 1);
    texLoc = glGetUniformLocation(shaderProgram, "maskTex");
    glUniform1i(texLoc, 2);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, tex[0]);
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, tex[1]);
    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_2D, tex[2]);
}

void Warper::setupVertexAttribArray() {
    GLint posAttrib = glGetAttribLocation(shaderProgram, "position");
    glEnableVertexAttribArray(posAttrib);
    glVertexAttribPointer(posAttrib, 2, GL_DOUBLE, GL_FALSE, 4 * sizeof(double), 0);
    GLint texAttrib = glGetAttribLocation(shaderProgram, "texcoord");
    glEnableVertexAttribArray(texAttrib);
    glVertexAttribPointer(texAttrib, 2, GL_DOUBLE, GL_FALSE, 4 * sizeof(double), (void*)(2 * sizeof(double)));
}

void Warper::bindTextures() {
    for (int i = 0; i < 3; ++i) {
        glGenTextures(1, &tex[i]);
        glBindTexture(GL_TEXTURE_2D, tex[i]);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        if (i < 2) {
            glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, input.w, input.h, 0, GL_RED, GL_FLOAT, textures[i]);
        } else {
            glTexImage2D(GL_TEXTURE_2D, 0, GL_R16UI, input.w, input.h, 0, GL_RED_INTEGER, GL_UNSIGNED_SHORT, maskTexture);
        }
    }
}

void Warper::setupFramebuffer() {
    // FBO
    glGenFramebuffers(1, &fbo);
    glGenRenderbuffers(1, &rbo[0]);
    glBindRenderbuffer(GL_RENDERBUFFER, rbo[0]);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_R32F, input.w, input.h);
    glBindFramebuffer(GL_FRAMEBUFFER, fbo);
    glFramebufferRenderbuffer(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_RENDERBUFFER, rbo[0]);

    glGenRenderbuffers(1, &rbo[1]);
    glBindRenderbuffer(GL_RENDERBUFFER, rbo[1]);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_R32F, input.w, input.h);
    glBindFramebuffer(GL_FRAMEBUFFER, fbo);
    glFramebufferRenderbuffer(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, GL_RENDERBUFFER, rbo[1]);

    glGenRenderbuffers(1, &rbo[2]);
    glBindRenderbuffer(GL_RENDERBUFFER, rbo[2]);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_R32UI, input.w, input.h);
    glBindFramebuffer(GL_FRAMEBUFFER, fbo);
    glFramebufferRenderbuffer(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT2, GL_RENDERBUFFER, rbo[2]);
}

void Warper::cleanup() {
    glDeleteProgram(shaderProgram);
	
    glDeleteBuffers(1, &ebo);
    glDeleteBuffers(1, &vbo);
    glDeleteVertexArrays(1, &vao);
    glDeleteFramebuffers(1, &fbo);
    glDeleteRenderbuffers(3, rbo);
    glfwTerminate();
}

// modified from cplusplus.com
void Warper::writeToFile(void** data, const string& filename, size_t elemSize, int imgSize) {
    FILE *pFile;
    fopen_s(&pFile, filename.c_str(), "wb");
    fwrite(*data, elemSize, imgSize, pFile);
    fclose(pFile);
}

void Warper::storeWarpedResults() {
    glReadBuffer(GL_COLOR_ATTACHMENT0);
    glReadPixels(0, 0, input.w, input.h, GL_RED, GL_FLOAT, targets[0]);
	writeToFile((void**)&targets[0], input.outputFITS[0], sizeof(float), input.w*input.h);
	glReadBuffer(GL_COLOR_ATTACHMENT1);
    glReadPixels(0, 0, input.w, input.h, GL_RED, GL_FLOAT, targets[1]);
	writeToFile((void**)&targets[1], input.outputFITS[1], sizeof(float), input.w*input.h);
	glReadBuffer(GL_COLOR_ATTACHMENT2);
    glReadPixels(0, 0, input.w, input.h, GL_RED_INTEGER, GL_UNSIGNED_SHORT, maskTarget);
	writeToFile((void**)&maskTarget, input.outputFITS[2], sizeof(unsigned short), input.w*input.h);
}