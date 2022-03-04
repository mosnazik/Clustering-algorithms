#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <vector>
#include "zada21.h"
//211-Ziatdinov-Nail
//File s function
int Point::newpoint(float ox, float oy){ //smeshaet tochki otnositelno centra
    X1=ox;
    X2=oy;
    return 0;
}
Point::Point(){

    X1=0;
    X2=0;
    Y=0;
    Num=0;
    Cl=0;
}
Point::Point(float x1, float x2, float y, int N, int Cl1){
    X1=x1;
    X2=x2;
    Y=y;
    Num=N;
    Cl=Cl1;
}
void Point::setNum(int N){
    Num=N;
}
float Point::getx1(){ //Vozvrashaet coordinatu 1

    return X1;
}
float Point::getx2(){ //Vozvrashaet coordinatu 2

    return X2;
}
float Point::gety(){ //Vozvrashaet znachenie v tochke
    return Y;
}
int Point::getN(){ //Vozvrashaet znachenie v tochke
    return Num;
}
Point & Point::operator=(Point arg){
    Cl=arg.Cl;
    X1=arg.getx1();
    X2=arg.getx2();
    Y=arg.gety();
    Num=arg.getN();
    return *this;
}
Cloud::Cloud(){
    N_points=0;
    O.newpoint(0,0);
}
Cloud::Cloud(float ox, float oy,float xdisper, float ydisper,int from, int num){
    if(LogOn==true)
        {log1<<"Sozdan cloud "<<num<<std::endl;}
    xdisp=xdisper;
    ydisp=ydisper;
    N_points=num;
    O.newpoint(ox,oy);
    //searchekstr();
}
Point Cloud::createpoint(float ox, float oy,float xdisper, float ydisper,int from, int Cl){
    int i,j,r;
    float x1, x2, y;
    for(i=0;i<1;i++){
        x1=0.0;
        x2=0.0;
        y=0.0;
        for (j = 0; j < 1000; j++)
        {
            r = rand()%1999-990;
            x1 = x1 + (float)(r)/998.5;
        }
        x1 = x1*xdisp/100.0;
        for (j = 0; j < 1000; j++)
        {
            r = rand()%1999-990;
            x2 = x2 + (float)(r)/998.5;
        }
        x2 = x2*ydisp/100.0;
        for (j = 0; j < 1000; j++)
        {
            r = rand()%1999-990;
            y = y + (float)(r)/998.5;
        }
        y = y/100.0;
    }

    Point x_Point(x1+ox,x2+oy,y,from, Cl);
    return x_Point;
}
int Cloud::searchekstr(){ //poisk ekstremum
    int i;
    float min1=100.0,min2=100.0,max1=-100.0,max2=-100.0;
    for(i=0;i<N_points;i++){
        if(m_Point[i]->getx1()<min1){min1=m_Point[i]->getx1();}
        if(m_Point[i]->getx1()>max1){max1=m_Point[i]->getx1();}
        if(m_Point[i]->getx2()<min2){min2=m_Point[i]->getx1();}
        if(m_Point[i]->getx2()>max2){max2=m_Point[i]->getx1();}
    }
    minx1=min1;
    minx2=min2;
    maxx1=max1;
    maxx2=max2;
    return 0;
}
int Cloud::Print_Cloud(std::ofstream& outdata,int CL){             //pechataet cloud
    if(LogOn==true)
        {log1<<"Start print cloud"<<std::endl;}
    int i;
    int r = rand() % 256, g = rand() % 256, b = rand() % 256;
	for(i=0;i<int(m_Point.size());i++)
    {
        //std::cout<<m_Point[i]->Cl<<std::endl;
        outdata<<m_Point[i]->getx1()<<" ";
        outdata<<m_Point[i]->getx2()<<" ";
        outdata<<r<<" "<<g<<" "<<b<<std::endl;
    }
    if(LogOn==true)
        {log1<<"Finish print cloud"<<std::endl;}
	return 0;
}
int Cloud::SetN(int N1){ //redactiruet kol-vo tochek v cloud

    N_points=N1;
    return 0;
}
void Cloud::getekstr(float &max1,float &max2,float &min1,float &min2){ //ustanavlivaet ekstremum
    max1=maxx1;
    max2=maxx2;
    min1=minx1;
    min2=minx2;
}
int Cloud::getN(){ //Vozvrashaet kol-vo tochek
        return N_points;
}
Buffer::Buffer(){
    Status_Work=false;
    N=0;
}
int Buffer::getN(){ //vozvrashaet kol-vo tochek
    return N;
}
int Buffer::upload(){  //zagrugaet oblako iz bufera
    m2_Point.clear();
    Status_Work=false;
    if(LogOn==true)
        {log1<<"Finish upload bufer"<<std::endl;}
    return 0;
}
int Buffer::rotateCl(float ang){     //povorachivaet oblako v buffere na ugol ang
    if(LogOn==true)
        {log1<<"Start bufer rotate"<<std::endl;}
    if(Status_Work==false)
    {
        std::cout << "Error" << std::endl;
        return -1;
    }
    int i;
    float s = sin(ang);
    float c = cos(ang);
    for(i=0;i<N;i++)
    {
        m2_Point[i].newpoint((m2_Point[i].getx1()-O.getx1())*c+(m2_Point[i].getx2()-O.getx2())*(s)+O.getx1(),(m2_Point[i].getx1()-O.getx1())*(-s)+(m2_Point[i].getx2()-O.getx2())*c+O.getx2());
    }
    if(LogOn==true)
        {log1<<"Finish bufer rotate"<<std::endl;}
    return 0;
}
int Buffer::shift(float ox, float oy){  //smeshaet oblako na vector(ox,oy)
    if(LogOn==true)
        {log1<<"Start bufer shift"<<std::endl;}
    if(Status_Work==false)
    {
        std::cout << "Error" << std::endl;
        return -1;
    }
    int i;
    for(i=0;i<N;i++)
    {
        m2_Point[i].newpoint(m2_Point[i].getx1()+ox,m2_Point[i].getx2()+oy);
    }
    O.newpoint(O.getx1()+ox,O.getx2()+oy);
    if(LogOn==true)
        {log1<<"Finish bufer shift"<<std::endl;}
    return 0;
}
int Buffer::zoom(float ly){             //sgimaet oblako v ly raz iz buffera
    if(LogOn==true)
        {log1<<"Start bufer zoom"<<std::endl;}
    if(Status_Work==false)
    {
        std::cout << "Error" << std::endl;
        return -1;
    }
    int i;
    for(i=0;i<N;i++)
    {
        m2_Point[i].newpoint((m2_Point[i].getx1()-O.getx1()),(m2_Point[i].getx2()-O.getx2()));
        m2_Point[i].newpoint((m2_Point[i].getx1()*ly),(m2_Point[i].getx2()*ly));
        m2_Point[i].newpoint((m2_Point[i].getx1()+O.getx1()),(m2_Point[i].getx2()+O.getx2()));
    }
    if(LogOn==true)
        {log1<<"Finish bufer zoom"<<std::endl;}
    return 0;
}
int Buffer::copyCl(Cloud* n_Cloud){      //copy oblako n_Cloud
    if(LogOn==true)
        {log1<<"Start bufer copy"<<std::endl;}
    N=n_Cloud->getN();
    O.newpoint(n_Cloud->O.getx1(),n_Cloud->O.getx2());
    Status_Work=true;
    if(LogOn==true)
        {log1<<"Finish bufer copy"<<std::endl;}
    return 0;
}
Field::Field(){
    N_points = 0;
    N_clouds = 0;
    Status_Work = true;
}
float Field::rastpoint(int n1, int n2){ //vichislyaet rastoyanie megdu tochkami n1 i n2
    int nn1, nn2,i;
    for(i=0;i<N_points;i++){
        if(m1_Point[i].getN()==n1){nn1=i;}
        if(m1_Point[i].getN()==n2){nn2=i;}
    }
    return sqrt((m1_Point[nn1].getx1()-m1_Point[nn2].getx1())*(m1_Point[nn1].getx1()-m1_Point[nn2].getx1())+(m1_Point[nn1].getx2()-m1_Point[nn2].getx2())*(m1_Point[nn1].getx2()-m1_Point[nn2].getx2()));
}
int Field::getCl(){     //vozvrashaet kol-vo cloud
    return N_clouds;
}
int Field::SetP(int N1){ //redactiruet kol-vo tochek v field
    N_points=N_points+N1;
    return 0;
}
int Field::work(){
    if(Status_Work==true){
        return 1;
    }
    else{
        return 0;
    }
}
int Field::createcloud(float ox, float oy,float xdisper, float ydisper, int num){
    int i;
    for(i=0;i<num;i++){
        Point y_Point;
        y_Point=m1_Cloud[m1_Cloud.size()-1].createpoint(ox,oy,xdisper,ydisper,N_points, m1_Cloud.size()-1);
        m1_Point.push_back(Point());
        m1_Point[N_points]=y_Point;
        m1_Cloud[m1_Cloud.size()-1].m_Point.push_back(&m1_Point[m1_Point.size()-1]);
        N_points++;
    }
    for(i=0;i<N_clouds;i++){
        for(int j=0;j<m1_Cloud[i].getN();++j){
            for(int k=0;k<int(m1_Point.size());++k){
                if((m1_Point[k].getN()==j)&&(m1_Point[k].Cl==i)){
                    m1_Cloud[i].m_Point[j]=&m1_Point[k];
                }
            }
        }
    }
    return 0;
}
int Field::Print_Field(std::ofstream& outdata){ //pechataet field
    int i;
    for(i=0; i<N_points;i++){
        log1<<m1_Point[i].getx1()<<'\t'<<m1_Point[i].getx2()<<std::endl;
    }
    for(i=0; i<N_clouds;i++){
        m1_Cloud[i].Print_Cloud(outdata,i);
    }
	return 0;
}
int Field::getN(){  //vozvrashaet kol-vo tochek
    return N_points;
}
int Field::SetCloud(int N1){   //dobavlyaet N1 tochek v chetchik v field
    N_clouds=N_clouds+N1;
    return 0;
}
int Field::savematrix_dist(){ //sozdaet matrix rastoyaniy
    if(LogOn==true)
        {log1<<"Start poisk matrix dist"<<std::endl;}
    int i,j;
    if(Status_Work==false){
            if(LogOn==true)
                {log1<<"Stop poisk"<<std::endl;}
            return 0;
    }
    m_Matrix_dist.setN(N_points);
    for(i=0;i<N_points;i++)
    {
        for(j=0;j<N_points;j++)
        {
            m_Matrix_dist.Mat_dist.push_back(0);
        }
    }

    for(i=0;i<N_points;i++)
    {
        for(j=0;j<N_points;j++)
        {
            m_Matrix_dist.Mat_dist.at(i*N_points+j)=rastpoint(i,j);
        }
    }
    if(LogOn==true){
        for(i=0;i<100;i++)
        {
            for(j=0;j<100;j++)
            {
                if(m_Matrix_dist.Mat_dist[i*N_points+j]<0.00001){
                    log1.width(8);
                    log1.precision(4);
                    log1<<m_Matrix_dist.Mat_dist[i*N_points+j]<<" "<<i<<"+"<<j<<" ";
                }
                else{
                    log1<<m_Matrix_dist.Mat_dist[i*N_points+j]<<" "<<i<<"+"<<j<<" ";
                }
            }
            log1<<'\n';
        }
    }
    Status_Work=false;
    if(LogOn==true)
        {log1<<"Finish poisk matrix dist"<<std::endl;}
	return 0;
}
