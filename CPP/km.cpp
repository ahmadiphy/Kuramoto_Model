#include "km.h"
using namespace std;
KM::KM()
{
    state=0;
}
//-----------------------------------------------------------------
//----------------   Kuramoto Model Sum Function   ----------------
//-----------------------------------------------------------------
double KM::Sum(double myTet)
{
    double result=0;
    for(int i=0;i<conTetas.size();++i)
    {
        double tetTar=conTetas[i];
        result+=sin(tetTar-myTet);
    }
    return result;
}
//-----------------------------------------------------------------
//---------------   Kuramoto Model Teta' Function   ---------------
//-----------------------------------------------------------------
double KM::fTeta(double thisTeta)
{
    double result=Omega*(K/cc)*Sum(thisTeta);
    return result;
}
//-----------------------------------------------------------------
//-----------------   step 1 in runge kutta   ---------------------
//-----------------------------------------------------------------
double KM::k1v(double thisTeta,double h)
{
    double result=h*fTeta(thisTeta);
    return result;
}
//------------------------------------------------------------------
//-----------------   Step 2 in Runge Kutta   ----------------------
//------------------------------------------------------------------
double KM::k2v(double thisTeta,double h)
{
    double result=h*fTeta(thisTeta+0.5*k1v(thisTeta,h));
    return result;
}
//------------------------------------------------------------------
//-----------------   Step 3 in Runge Kutta   ----------------------
//------------------------------------------------------------------
double KM::k3v(double thisTeta,double h)
{
    double result=h*fTeta(thisTeta+0.5*k2v(thisTeta,h));
    return result;
}
//-------------------------------------------------------------------
//------------------   Step 4 in Runge Kutta   ----------------------
//-------------------------------------------------------------------
double KM::k4v(double thisTeta,double h)
{
    double result=h*fTeta(thisTeta+k3v(thisTeta,h));
    return result;
}
//-------------------------------------------------------------------
//--------------------   Runge Kutta Result   -----------------------
//-------------------------------------------------------------------
double KM::RKteta(double thisTeta)
{
    double result=thisTeta+0.166666667*(k1v(thisTeta,h)+2*k2v(thisTeta,h)+2*k3v(thisTeta,h)+k4v(thisTeta,h));
    return result;
}
//-------------------------------------------------------------------
//---------------   Initialize and Connect Neurons   ----------------
//-------------------------------------------------------------------
void KM::Initialize(double omega , double teta,std::vector<int> &cv)
{
    Teta=teta;
    tempTeta=teta;
    Omega=omega;
    for(int i=0;i<cv.size();++i)
    {
        conns.push_back(cv[i]);
        conTetas.push_back(0);
    }
    conns.shrink_to_fit();
    conTetas.shrink_to_fit();
    //Omega=conns.size();
    state=1;
}
//--------------------------------------------------------------------
//-------------------   Run Neural Dynamics   ------------------------
//--------------------------------------------------------------------
void KM::Run()
{
    if(state==1)
    {
        tempTeta=RKteta(Teta);//solve v equation RK4
    }else
    {
        cout<<"   ERROR! (initial value problem)  "<<endl;
        exit (EXIT_FAILURE);
    }
}
//--------------------------------------------------------------------
//--------------------   Get Random Current   ------------------------
//--------------------------------------------------------------------
void KM::Update(vector<double> &tetS)
{
    Teta=tempTeta;
    tempTeta=0;
    for(int i=0;i<conTetas.size();++i)
    {
        conTetas[i]=tetS[conns[i]];
    }
}
