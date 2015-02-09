#include "cLatticeModel.h"
#include <random>
#include <iostream>
#include <fstream>
std::default_random_engine generator;

cLatticeModel::cLatticeModel(std::vector<int> in_protein) {
    this->N = in_protein.size()+6;
    

    initLattice();
    
    for(int i=0;i<in_protein.size();i++){
        cAminoAcid aa(Vector<int>(i,0,0),i,in_protein[i]);
        protein.push_back(aa);
        setResidueFromIdx(i,0,0,i);
    }
    
    checkProteinPosition();
    RREnergy = NULL;
    loadRREnergy("IRRMatrix.dat");

}

cLatticeModel::cLatticeModel(std::string& filename){
    std::ifstream myfile (filename.c_str());
    int protSize, protSizeY;
    myfile>>protSize;
    myfile>>protSizeY;
    if(protSize!=protSizeY){
        std::cout<<"Error!\n"<<std::endl;
        exit(1);
    }
    
    this->N = protSize+6;
    initLattice();
    
    contactMap = new int *[N];
    for( int i=0;i<protSize;i++){
        contactMap[i] = new int [N];
        for( int j=0;j<protSize;j++){
            contactMap[i][j] = 0;
            myfile>>contactMap[i][j];
        }
    }
    myfile.close();
    
    for(int i=0;i<protSize;i++){
        cAminoAcid aa(Vector<int>(i,0,0),i,i);
        protein.push_back(aa);
        setResidueFromIdx(i,0,0,i);
    }
    
    checkProteinPosition();
    
}

void cLatticeModel::initLattice(){
    //Basis vectors for FCC
    baseVectors.push_back(Vector<double>(1,0,0));
    baseVectors.push_back(Vector<double>(0.5,sqrt(3)/2.0,0));
    baseVectors.push_back(Vector<double>(0.5,sqrt(3)/6.0,sqrt(6)/3.0));

    //Neighbour directions for FCC
    baseDirections.push_back(Vector<int>(1,0,0));
    baseDirections.push_back(Vector<int>(0,1,0));
    baseDirections.push_back(Vector<int>(0,0,1));
    baseDirections.push_back(Vector<int>(-1,0,0));
    baseDirections.push_back(Vector<int>(0,-1,0));
    baseDirections.push_back(Vector<int>(0,0,-1));

    baseDirections.push_back(Vector<int>(1,-1,0));
    baseDirections.push_back(Vector<int>(1,0,-1));
    baseDirections.push_back(Vector<int>(0,1,-1));
    baseDirections.push_back(Vector<int>(-1,1,0));
    baseDirections.push_back(Vector<int>(-1,0,1));
    baseDirections.push_back(Vector<int>(0,-1,1));

    occupationTable = new int **[N];
    for( int i=0;i<N;i++){
        occupationTable[i] = new int *[N];
        for( int j=0;j<N;j++){
            occupationTable[i][j] = new int [N];
            for(int k=0;k<N;k++){
                occupationTable[i][j][k] = -1;
            }
        }
    }
    
}
cLatticeModel::~cLatticeModel() {
    for( int i=0;i<N;i++){
        for( int j=0;j<N;j++)
            delete [] occupationTable[i][j];
        delete [] occupationTable[i];
    }
    delete [] occupationTable;
    if(RREnergy!=NULL){
        for( int i=0;i<21;i++){
            delete [] RREnergy[i];
        }
        delete [] RREnergy;
    }
}

void cLatticeModel::loadRREnergy(const std::string filename){
    
    RREnergy = new double *[21];
    for( int i=0;i<21;i++){
        RREnergy[i] = new double [21];
        for( int j=0;j<21;j++){
            RREnergy[i][j] = 0;
        }
    }
    std::ifstream myfile (filename.c_str());
    double en;
    int a,b;
    for( int i=0;i<210;i++){
        myfile>>a;
        myfile>>b;
        myfile>>en;
        RREnergy[a][b]=en;
        RREnergy[b][a]=en;
    }
    myfile.close();
}

void cLatticeModel::display(){
    static int time;
    glPointSize(10);
    glLineWidth(4);
    /*glBegin(GL_POINTS);
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            for(int k=0;k<N;k++){
                if(getResidueFromIdx(i, j, k)==-1){
                    Vector<double> p = baseVectors[0]*i + baseVectors[1]*j+ baseVectors[2]*k;
                    glVertex3f(p.x,p.y,p.z);
                }
            }
        }
    }
    glEnd();*/
    
    glBegin(GL_LINES);
 
    if(time%20==0){
        //endMove();
        //crankShaftMove();
        //snakeMove();
    }
    glColor4f(1.f, 1.f, 1.f, 1.0f);
    for(int i=0;i<protein.size()-1;i++){
        Vector<double> p1 = baseVectors[0]*(protein[i].pos.x)+baseVectors[1]*(protein[i].pos.y)+baseVectors[2]*(protein[i].pos.z);
        Vector<double> p2 = baseVectors[0]*(protein[i+1].pos.x)+baseVectors[1]*(protein[i+1].pos.y)+baseVectors[2]*(protein[i+1].pos.z);
        glVertex3f(p1.x,p1.y,p1.z);
        glVertex3f(p2.x,p2.y,p2.z);
    }

    glEnd();
    
    glBegin(GL_POINTS);
    for(int i=0;i<protein.size();i++){
        /*if(protein[i].type==1)
            glColor4f(0.0f, 0.0f, 0.0f, 1.0f);
        else
            glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
        */
        protein[i].setColor();
        Vector<double> p1 = baseVectors[0]*(protein[i].pos.x)+baseVectors[1]*(protein[i].pos.y)+baseVectors[2]*(protein[i].pos.z);
        glVertex3f(p1.x,p1.y,p1.z);
    }
    glEnd();
    
    time+=1;
    //std::cout<<time<<std::endl;
}
bool cLatticeModel::isInTable(int i, int j, int k) {
    return ( i>=0 && i<N ) && ( j>=0 && j<N ) && ( k>=0 && k<N );
}
bool cLatticeModel::isInTable(const Vector<int> &pos) {
    return ( pos.x>=0 && pos.x<N ) && ( pos.y>=0 && pos.y<N ) && ( pos.z>=0 && pos.z<N );
}
int cLatticeModel::getResidueFromIdx(int i, int j, int k) {
    return occupationTable[i][j][k];
}
int cLatticeModel::getResidueFromIdx(const Vector<int>& pos){
    return occupationTable[pos.x][pos.y][pos.z];
}
void cLatticeModel::setResidueFromIdx(int i, int j, int k, int idx) {
    occupationTable[i][j][k] = idx;
}
void cLatticeModel::setResidueFromIdx(const Vector<int>& pos, int idx){
    occupationTable[pos.x][pos.y][pos.z] = idx;
}
bool cLatticeModel::areAdjacent(int i1, int j1, int k1, int i2, int j2, int k2) {
    for(int i=0;i<baseDirections.size();i++){
        if (Vector<int>(i1-i2,j1-j2,k1-k2) == baseDirections[i]){
            return true;
        }
    }
    return false;
}
bool cLatticeModel::areAdjacent(const Vector<int>& pos1, const Vector<int>& pos2){
    Vector<int> d = pos2-pos1;
    for(int i=0;i<baseDirections.size();i++){
        if (d==baseDirections[i])
            return true;
    }
    return false;
}

void cLatticeModel::moveResidue(const int idx, const Vector<int>& pos){
    setResidueFromIdx(protein[idx].pos.x, protein[idx].pos.y,protein[idx].pos.z, -1);
    setResidueFromIdx(pos.x, pos.y, pos.z, idx);
    protein[idx].pos=pos;
}

int cLatticeModel::endMove() {
    std::uniform_int_distribution<int> distribution(1,2);
    int end_choice = distribution(generator);
    int idx=0, idx_next=1;
    if(end_choice==1){ //moving the beginning residue
        idx=0;
        idx_next=1;
    }else{ //moving ending residue
        idx=protein.size()-1;
        idx_next=protein.size()-2;
    }
    std::vector<Vector<int>> emptyPositions;
    for(int i=0;i<baseDirections.size();i++){
        Vector<int> n_pos = protein[idx].pos+baseDirections[i];
        if(getResidueFromIdx(n_pos)==-1){
            if(areAdjacent(n_pos, protein[idx_next].pos))
                emptyPositions.push_back(n_pos);
        }
    }
    if(emptyPositions.size()==0){
        return -1;
    }
    std::uniform_int_distribution<int> neighbour_distribution(0,emptyPositions.size()-1);
    int neighbour_choice = neighbour_distribution(generator);
    moveResidue(idx, emptyPositions[neighbour_choice]);    
    checkProteinPosition();
}

int cLatticeModel::crankShaftMove() {
    std::uniform_int_distribution<int> distribution(1,protein.size()-1);
    int idx1 = distribution(generator);
    int idx0 = idx1-1;
    int idx2 = idx1+1;

    std::vector<Vector<int>> emptyPositions;
    for(int i=0;i<baseDirections.size();i++){
        Vector<int> n_pos = protein[idx0].pos+baseDirections[i];
        if(getResidueFromIdx(n_pos)==-1){
            if(areAdjacent(n_pos, protein[idx2].pos))
                emptyPositions.push_back(n_pos);
        }
    }
    if(emptyPositions.size()==0){
        return -1;
    }
    std::uniform_int_distribution<int> neighbour_distribution(0,emptyPositions.size()-1);
    int neighbour_choice = neighbour_distribution(generator);
    moveResidue(idx1, emptyPositions[neighbour_choice]);
    checkProteinPosition();
}

int cLatticeModel::snakeMove() {
    std::uniform_int_distribution<int> distribution(1,2);
    int direction = distribution(generator);
    int idx;
    if( direction == 1){
        idx=0;
    }else{
        idx=protein.size()-1;
    }
    std::vector<Vector<int>> emptyPositions;
    for(int i=0;i<baseDirections.size();i++){
        Vector<int> n_pos = protein[idx].pos+baseDirections[i];
        if(getResidueFromIdx(n_pos)==-1){
           emptyPositions.push_back(n_pos);
        }
    }
    if(emptyPositions.size()==0){
        return -1;
    }
    std::uniform_int_distribution<int> neighbour_distribution(0,emptyPositions.size()-1);
    int neighbour_choice = neighbour_distribution(generator);
    Vector<int> nextPos = protein[idx].pos;
    moveResidue(idx,emptyPositions[neighbour_choice]);
    
    if(direction==1) {
        for (int i = 1; i < protein.size(); i++) {
            Vector<int> resPosOld = protein[i].pos;
            moveResidue(i,nextPos);
            nextPos = resPosOld;
        }
    }else{
        for (int i = protein.size()-2; i > -1; i--) {
            Vector<int> resPosOld = protein[i].pos;
            moveResidue(i,nextPos);
            nextPos = resPosOld;
        }
    }
    checkProteinPosition();
}

int cLatticeModel::pullMove() {
    std::uniform_int_distribution<int> distributionDirections(1,2);
    std::uniform_int_distribution<int> distributionNodes(1,protein.size()-2);
    int direction = 2;//distributionDirections(generator);
    int idx = distributionNodes(generator);
    int idx1;
    if(direction==1){
        idx1=idx-1;
    }else{
        idx1=idx+1;
    }
    std::vector<Vector<int>> emptyPositions;
    for(int i=0;i<baseDirections.size();i++){
        Vector<int> n_pos = protein[idx].pos+baseDirections[i];
        if(getResidueFromIdx(n_pos)==-1){
            if(areAdjacent(n_pos, protein[idx1].pos)){
                emptyPositions.push_back(n_pos);
            }
        }
    }
    if(emptyPositions.size()==0){
        return -1;
    }
    std::uniform_int_distribution<int> neighbour_distribution(0,emptyPositions.size()-1);
    int neighbour_choice = neighbour_distribution(generator);
    Vector<int> nextPos = protein[idx].pos;
    moveResidue(idx,emptyPositions[neighbour_choice]);
    
    if(direction==1) {
        for (int i = 1; i < protein.size(); i++) {
            Vector<int> resPosOld = protein[i].pos;
            moveResidue(i,nextPos);
            if(i>2){
                if(areAdjacent(nextPos,protein[i+1].pos)){
                    break;
                }
            }
            nextPos = resPosOld;
            
        }
    }else{
        for (int i = protein.size()-2; i > -1; i--) {
            Vector<int> resPosOld = protein[i].pos;
            moveResidue(i,nextPos);
            if(i>0){
                if(areAdjacent(nextPos,protein[i-1].pos)){
                    break;
                }
            }
            nextPos = resPosOld;
        }
    }
    checkProteinPosition();
}

Vector<int> cLatticeModel::getProteinCM(){
    Vector<int> mCenter(0,0,0);
    for(int i=0;i<protein.size();i++){
        mCenter+=protein[i].pos;
    }
    mCenter/=protein.size();
    return mCenter;
}

void cLatticeModel::moveProtein(Vector<int>& dPos){
    for(int i=0;i<protein.size();i++){
        moveResidue(i,protein[i].pos+dPos);
    }
}

void cLatticeModel::checkProteinPosition(){
    Vector<int> cMass = getProteinCM();
    Vector<int> latCenter(N/2,N/2,N/2);
    Vector<int> dR = latCenter-cMass;
    moveProtein(dR);
}

cLatticeModel& cLatticeModel::operator=(const cLatticeModel &other) {
    for( int i=0;i<N;i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                occupationTable[i][j][k]=other.occupationTable[i][j][k];
            }
        }
    }
    protein = other.protein;
}

double cLatticeModel::computeEnergy() {
    
    //return computeEnergyContacts();
    return computeEnergyRR();
}

double cLatticeModel::computeEnergyRR(){
    double energy=0.0;
    for(int idx=0;idx<protein.size();idx++){
        if(protein[idx].type==2) continue;
        for(int i=0;i<baseDirections.size();i++){
            Vector<int> n_pos = protein[idx].pos + baseDirections[i];
            int n_idx = getResidueFromIdx(n_pos);
            if(n_idx==-1 || abs(n_idx-idx)<2) continue;
            energy+=RREnergy[protein[idx].type][protein[n_idx].type];
        }
    }
    return energy;
}

double cLatticeModel::computeEnergyHP(){
    double energy=0.0;
    for(int idx=0;idx<protein.size();idx++){
        if(protein[idx].type==2) continue;
        for(int i=0;i<baseDirections.size();i++){
            Vector<int> n_pos = protein[idx].pos + baseDirections[i];
            int n_idx = getResidueFromIdx(n_pos);
            if(n_idx==-1 || abs(n_idx-idx)<2) continue;
            if( protein[n_idx].type==1 ){
                energy+=-1;
            }
        }
    }
    return energy;
}

double cLatticeModel::computeEnergyContacts(){
    double energy=0.0;
    for(int idx=0;idx<protein.size();idx++){
        for(int i=0;i<baseDirections.size();i++){
            Vector<int> n_pos = protein[idx].pos + baseDirections[i];
            int n_idx = getResidueFromIdx(n_pos);
            if(n_idx==-1 || abs(n_idx-idx)<2) continue;
            if( contactMap[idx][n_idx]==1){
                energy+=-1;
            }
        }
    }
    return energy;
}

Vector<double> cLatticeModel::getProteinCenterPos(){
    Vector<int> cM = getProteinCM();
    return baseVectors[0]*cM.x+baseVectors[1]*cM.y+baseVectors[2]*cM.z;
}

#define AA_ORANGE 0.9f, 0.6f, 0.0f, 1.0f
#define AA_BLUE 0.0f, 0.0f, 0.9f, 1.0f
#define AA_RED 0.9f, 0.0f, 0.0f, 1.0f
#define AA_MAGNETTA 0.9f, 0.0f, 0.9f, 1.0f
#define AA_GREEN 0.0f, 0.9f, 0.0f, 1.0f
#define AA_BLACK 0.0f, 0.0f, 0.0f, 1.0f
void cAminoAcid::setColor() const{
    switch(this->type){
        case 0:
            glColor4f(AA_BLACK);
        case 1: case 6: case 16: case 17:
            glColor4f(AA_ORANGE);
            return;
        case 3: case 4:
            glColor4f(AA_BLUE);
            return;
        case 9: case 15:
            glColor4f(AA_RED);
            return;
        case 7: case 12: case 14:
            glColor4f(AA_MAGNETTA);
            return;
        default:
            glColor4f(AA_GREEN);
            return;
    }
    
}