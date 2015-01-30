#pragma once
#include "GlutFramework.h"

using namespace glutFramework;

class cAminoAcid{
public:
    cAminoAcid(const Vector<int> &pos, int ind, int type){
      this->pos=pos;this->ind = ind;this->type = type;
    };
    ~cAminoAcid(){;};

    Vector<int> pos;
    int ind;
    int type;

    cAminoAcid &operator=( const cAminoAcid &other ){
        pos=other.pos;
        ind=other.ind;
        type=other.type;
    }
};

class cLatticeModel: public Object{

public:
    cLatticeModel(std::vector<int> protein);
    cLatticeModel(std::string &filename);
    ~cLatticeModel();
    void display();

    int endMove();
    int crankShaftMove();
    int snakeMove();
    

    cLatticeModel &operator=( const cLatticeModel &other );

    double computeEnergy();
    
    Vector<double> getProteinCenterPos();

private:
    int N;

    std::vector<Vector<double>> baseVectors;
    std::vector<Vector<int>> baseDirections;

    std::vector<cAminoAcid> protein;
    int ***occupationTable;
    
    int **contactMap;
    
    void initLattice();

    bool isInTable(int i, int j, int k);
    bool isInTable(const Vector<int> &pos);
    int getResidueFromIdx(int i, int j, int k);
    int getResidueFromIdx(const Vector<int> &pos);
    void setResidueFromIdx(int i, int j, int k, int idx);
    void setResidueFromIdx(const Vector<int> &pos, int idx);
    bool areAdjacent(int i1, int j1, int k1, int i2, int j2, int k2);
    bool areAdjacent(const Vector<int> &pos1, const Vector<int> &pos2);
    void moveResidue(const int idx, const Vector<int> &pos);
    //corrects position of protein so that it fits occupation table and leaves at least one cell from the border
    Vector<int> getProteinCM();
    void moveProtein(Vector<int> &dPos);
    void checkProteinPosition();
    
    double computeEnergyHP();
    double computeEnergyContacts();
};