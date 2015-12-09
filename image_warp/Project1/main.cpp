#include <fstream>
#include <iostream>
#include "Warper.h"

using namespace std;

extern "C" {
    _declspec(dllexport) int NvOptimusEnablement = 0x00000001;
}

int main(int argc, char **argv) {
	FILE *stats;
	fopen_s(&stats, "testdata/stats", "a");
	ofstream fs(stats);
    WarpData in;
    in.w = 2046;
    in.h = 4094;
    in.grid_w = 206;
    in.grid_h = 411;
    in.kernelSize = 64;
    in.fetchMethod = BILINEAR;
    in.inputCoords[0] = "testdata/sourceCoords.raw";
    in.inputCoords[1] = "testdata/destinationCoords.raw";
    in.inputFITS[0] = "testdata/input.img";
    in.inputFITS[1] = "testdata/input.var";
    in.inputFITS[2] = "testdata/input.mask";
    in.outputFITS[0] = "testdata/outImg";
    in.outputFITS[1] = "testdata/outVar";
    in.outputFITS[2] = "testdata/outMask";

    Warper warper(in);
	for (int i = 0; i < 30; ++i) {
	
        warper.warp(fs);

	}
    float psnr = warper.calc_PSNR("testdata/output.img");
    fs << psnr << endl;
	fclose(stats);
	return 0;
}