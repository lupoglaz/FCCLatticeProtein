#include "GL/GlutFramework.h"
#include "cLatticeModel.h"
#include "cMonteCarloMover.h"
using namespace glutFramework;

void runGLUT(int argc, char** argv){
    GlutFramework framework;

    std::vector<int> protein{1,2,1,1,       //HPH2
                            2,1,1,1,1,      //PH4
                            2,1,1,1,2,2,     //PH3P2
                            1,2,1,2,1,2,     //HPHPHP
                            1,2,1,2,1,2,     //HPHPHP
                            1,2,2,2,2,2,2,2,2,//HP8
                            1};             //H
    std::vector<int> proteinRR{1, 14, 15, 11, 16, 19, 9,15,1,1,16,9,1};
    cMonteCarloMover MMC(proteinRR, 10.0);
    //cMonteCarloMover MMC("1ATZ_distMat.dat", 10.0);
    //cLatticeModel lMdl(protein);
    framework.addObject(&MMC);
    Vector<double> lookAtPos = MMC.getProteinCenterPos();
    framework.setLookAt(0.0, 2.0, 10.0,lookAtPos.x, lookAtPos.y, lookAtPos.z, 0.0, 1.0, 0.0);
    framework.startFramework(argc, argv);
}

int main(int argc, char** argv){
    runGLUT(argc,argv);
    return 0;
}
