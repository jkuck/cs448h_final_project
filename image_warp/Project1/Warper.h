#include "GLFW\glfw3.h"

using namespace std;

enum VERTEX_FETCH {
    DIRECT, BILINEAR
};

typedef struct WarpData {
    int w, h; // width and height of image
    int grid_w, grid_h; // coordinate inputs are grid_w x grid_h grid points
    int kernelSize; // kernel size
    VERTEX_FETCH fetchMethod; // method of fetching, can be direct or using bilinear
    string inputFITS[3]; // input filenames for the three planes
    string outputFITS[3]; // output filenames for the three planes
    string inputCoords[2]; // grid point (i, j) would be stored in the i+j*grid_w entry
} WarpData;

class Warper {
public:
    Warper(WarpData input);
    ~Warper();
    void warp(ofstream& fs);
    float calc_PSNR(string refOutputFile);
    void storeWarpedResults();

    WarpData get_input() { return input; }
    bool hasWarped() { return warped; }
    int get_kernel_size() { return input.kernelSize; };
	int set_kernel_size(int size) {
        input.kernelSize = size;
        warped = false; 
    }
    int get_image_width() { return input.w; };
    int get_image_height() { return input.h; };

private:
    bool warped;
    WarpData input;
	// functions and variables used for calculating PSNR
    float psnr(float* ref, float* result);
    float *targetReference;

	// functions and variables used for reading data / warping
    GLFWwindow* initGLFW();
    void readInput(void **buffer, const string& filename);
    void loadInputs();
    void constructArrays();
    void bindBuffers();
    void loadShader(const string& file, GLenum vOrF);
    void createShaderProgram();
    void setUniforms();
    void setupVertexAttribArray();
    void bindTextures();
    void setupFramebuffer();
    void writeToFile(void** data, const string& filename, size_t elemSize, int imgSize);
    void cleanup();


    float* textures[2];
    unsigned short* maskTexture;
    float* targets[2];
    unsigned short *maskTarget;

    double *templateCoords;
    double *scienceCoords;
    double* coords;
    GLuint* gridElements;

    GLuint vao, vbo, ebo;
    GLuint fbo;
    GLuint rbo[3];

    char *vShaderFile, *fShaderFile; // shader filenames, dependent on the type of vertex fetching
    GLuint shaderProgram;
    GLuint vertShader, fragShader;
    GLuint tex[3];
};