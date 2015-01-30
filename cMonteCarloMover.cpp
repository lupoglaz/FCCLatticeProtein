#include "cMonteCarloMover.h"
extern std::default_random_engine generator;
cMonteCarloMover::cMonteCarloMover(std::vector<int> protein, double T){
    curModel = new cLatticeModel(protein);
    trialModel = new cLatticeModel(protein);
    minModel = new cLatticeModel(protein);
    E = curModel->computeEnergy();
    Emin=E;
    this->T = T;

    std::random_device rd;
    mt = new std::mt19937(rd());
    distr = new std::uniform_real_distribution<double>(0.0, 1.0);

}

cMonteCarloMover::cMonteCarloMover(std::string filename, double T){
    curModel = new cLatticeModel(filename);
    trialModel = new cLatticeModel(filename);
    minModel = new cLatticeModel(filename);
    E = curModel->computeEnergy();
    Emin=E;
    this->T = T;

    std::random_device rd;
    mt = new std::mt19937(rd());
    distr = new std::uniform_real_distribution<double>(0.0, 1.0);
}

cMonteCarloMover::~cMonteCarloMover() {
    delete curModel;
    delete trialModel;
    delete mt;
    delete distr;
}

void cMonteCarloMover::move() {
    (*trialModel) = (*curModel);

    std::uniform_int_distribution<int> distribution(1,3);
    switch(distribution(generator)){
        case 1:
            trialModel->endMove();
            break;
        case 2:
            trialModel->crankShaftMove();
            break;
        case 3:
            trialModel->snakeMove();
            break;
    }
    
    double E1 = trialModel->computeEnergy();
    if(E1<=E){
    
        (*curModel)=(*trialModel);
        //std::cout<<"E1 = "<<E1<<" E ="<<E<<std::endl;
        E=E1;
        
    
    }else{
        double dE = E1-E;
        double rndNum = (*distr)(*mt);
        double threshold = exp(-dE/T);
        if( rndNum <  threshold){
            (*curModel)=(*trialModel);
            E=E1;
            //std::cout<<E<<" dE = "<<dE<<" rndNum = "<<rndNum<<" Threshold = "<<threshold<<std::endl;

        }
    }
    if( E<Emin){
        (*minModel)=(*curModel);
        Emin=E;
        std::cout<<"Emin = "<<Emin<<std::endl;
    }
    
    //std::cout<<"Emin = "<<Emin<<std::endl;

}

void cMonteCarloMover::display() {
    move();
    minModel->display();
    //curModel->display();

}