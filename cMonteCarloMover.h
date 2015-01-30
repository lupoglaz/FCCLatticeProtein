#pragma once

#include "cLatticeModel.h"
#include <random>

class cMonteCarloMover: public Object{
public:
    cMonteCarloMover(std::vector<int> protein, double T);
    cMonteCarloMover(std::string filename, double T);
    ~cMonteCarloMover();

    void display();
    Vector<double> getProteinCenterPos(){return minModel->getProteinCenterPos();};

private:
    cLatticeModel *curModel;
    cLatticeModel *trialModel;
    cLatticeModel *minModel;

    std::uniform_real_distribution<double> *distr;
    std::mt19937 *mt;

    double E, Emin;
    double T;

    void move();
};