#ifndef KM_H
#define KM_H
#include "in.h"
#include"global_variables.h"

class KM
{
public:
    KM();
    int myNum;
    double Teta,Omega,tempTeta;
    std::vector<int> conns;//connections of network
    std::vector<double> conTetas;//connections tetas of network
    //std::vector<double> localKs;//local Ks
    void Initialize(double omega , double teta, std::vector<int> &cv);
    void Run();
    void Update(std::vector<double> &tetS);
private:
    int state;
    double Sum(double myTet);
    double fTeta(double thisTeta);
    double k1v(double thisTeta, double h);
    double k2v(double thisTeta, double h);
    double k3v(double thisTeta, double h);
    double k4v(double thisTeta, double h);
    double RKteta(double thisTeta);
};

#endif // KM_H
