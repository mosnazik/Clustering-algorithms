#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <vector>
#include "zada21.h"
#define eps 0.000000001
#define PI 3.14159265
//211-Ziatdinov-Nail
//File s function
arr::arr(){
    N=0;
}
int arr::getN(){
    return N;
}
int arr::setN(int N1){
    N=N1;
    return 0;
}
int arr::operator[]( int k){
        if((k>= N )||(k<0)){
            if(LogOn==true)
                {log1<<"Oshibka v massive a:"<<k<<" size:"<<N<<'\n';}
        }
        return mass[k];
}
void arr::operator+=(const int k){
    mass.push_back(k);
    m.push_back(0);
    N++;
}
Launch::Launch(){
    N_points=0;
    N_cluster=0;
}
Launch::Launch(float delt, int k1){
    N_points=0;
    N_cluster=0;
    k=k1;
    d=delt;
}
int Launch::Print_Cluster(std::ofstream& outdata){  //pechataet cluster
    if(LogOn==true)
        {log1<<"Start print cluster"<<std::endl;}
    int j;
    for(j=0;j<N_cluster;j++){
        m1_Cluster[j]->Print_Cluster(outdata);
    }
    if(LogOn==true)
        {log1<<"Finish print cluster"<<std::endl;}
    return 0;
}
int Launch::Print_tree(std::ofstream& outdata){
    if(LogOn==true)
        {log1<<"Start print tree"<<std::endl;}
    m_Min_tree.Print_tree(outdata);
    if(LogOn==true)
        {log1<<"Finish print tree"<<std::endl;}
    return 0;
}
int Launch::Start_alg(int n){  // zapusk volnovogo algoritma
    if(LogOn==true)
        {log1<<"Start algoritm "<<n<<std::endl;}
    int i;
    if(n==1){
        for(i=0;i<m_Wave_alg.m_Field->getN();i++){
            m_Wave_alg.param.push_back(0);
        };
        m_Wave_alg.savematrix_graf(d);
        m_Wave_alg.Start_alg(d);
        N_cluster=m_Wave_alg.getCl();
        N_points=m_Wave_alg.getN();
        for(i=0;i< N_cluster;i++){
            m1_Cluster.push_back(&m_Wave_alg.m_Cluster[i]);
        }
        return 0;
    }
    if(n==2){
        for(i=0;i<m_DBScan.m1_Wave_alg.m_Field->getN();i++){
            m_DBScan.m1_Wave_alg.param.push_back(-1);
        };
        m_DBScan.Start_alg(d,k);
        N_cluster=m_DBScan.getCl();
        N_points=m_DBScan.getN();
        for(i=0;i< N_cluster;i++){
            m1_Cluster.push_back(&m_DBScan.m1_Wave_alg.m_Cluster[i]);
        }
    }
    if(n==3){
        m_Min_tree.Start_alg();
        N_points=m_Min_tree.getN();
        return 0;
    }
    if(n==4){
        m_Min_tree.ClusterS();
        N_cluster=m_Min_tree.getCl();
        for(i=0;i< N_cluster;i++){
            m1_Cluster.push_back(&m_Min_tree.m_Cluster[i]);
        }
        return 0;
    }
    if(n==5){
        m_k_means.creat_matr(k);
        m_k_means.Start_alg(k);
        N_points=m_k_means.getN();
        N_cluster=m_k_means.getCl();
        for(i=0;i< N_cluster;i++){
            m1_Cluster.push_back(&m_k_means.m_Cluster[i]);
        }
        return 0;
    }
    if(n==6){
        m_EM.creat_matr(k);
        m_EM.Start_alg(k);
        N_points=m_EM.getN();
        N_cluster=m_EM.getCl();
        for(i=0;i< N_cluster;i++){
            m1_Cluster.push_back(&m_EM.m_Cluster[i]);
        }
        return 0;
    }
    if(LogOn==true)
        {log1<<"Finish algoritm"<<std::endl;}
    return 0;
}
int Launch::SetP(int N1){  //redactiruet kol-vo tochek v launch
    N_points=N1;
    return 0;
}
int Launch::getCl(){  //Vozvrashaet kol-vo clusterov
    return N_cluster;
}
int Launch::SetCluster(int N1){  //redactiruet kol-vo clusterov v cloud
    N_cluster=N1;
    return 0;
}
int Launch::getN(){    //Vozvrashaet kol-vo tochek
    return N_points;
}
Cluster::Cluster(){
    N_points=0;
}
int Cluster::Print_Cluster(std::ofstream& outdata){ //pechataet cluster
    int i;
    int r = (rand() % 128)+(rand() % 128), g = rand() % 256, b = rand() % 256;
    for(i=0;i<N_points;i++)
    {
        outdata<<m_Point[i]->getx1()<<" ";
        outdata<<m_Point[i]->getx2()<<" ";
        outdata<<r<<" "<<g<<" "<<b<<std::endl;
    }
    if(N_points>10){
        outdata<<O.getx1()<<" ";
        outdata<<O.getx2()<<" ";
        outdata<<0<<" "<<0<<" "<<0<<std::endl;
    }
    return 0;
}
int Cluster::getN(){  //Vozvrashaet kol-vo tochek
    return N_points;
}
void Cluster::getekstr(float &max1,float &max2,float &min1,float &min2){ //ustanavlivaet ekstremum
    max1=maxx1;
    max2=maxx2;
    min1=minx1;
    min2=minx2;
}
int Cluster::searchekstr(){
    int i;
    float min1,min2,max1,max2,ox1=0.0,ox2=0.0;
    min1=min2=1000.0;
    max1=max2=-1000.0;
    for(i=0;i<N_points;i++){
        if(m_Point[i]->getx1()<min1){min1=m_Point[i]->getx1();}
        if(m_Point[i]->getx1()>max1){max1=m_Point[i]->getx1();}
        if(m_Point[i]->getx2()<min2){min2=m_Point[i]->getx1();}
        if(m_Point[i]->getx2()>max2){max2=m_Point[i]->getx1();}
        ox1=m_Point[i]->getx1()+ox1;
        ox2=ox2+m_Point[i]->getx1();
    }
    minx1=min1;
    minx2=min2;
    maxx1=max1;
    maxx2=max2;
    ox1=ox1/float(N_points);
    ox2=ox2/float(N_points);
    if(N_points>10){O.newpoint(ox1,ox2);}
    return 0;
}
int Cluster::SetP(int N1){ //redactiruet kol-vo tochek v cluster
    N_points=N1;
    return 0;
}
Matrix_dist::Matrix_dist(){

}
int Matrix_dist::getN(){   //Vozvrashaet kol-vo tochek
    return N1;
}
void Matrix_dist::setN(int n){ //redactiruet kol-vo tochek v matrix distance
    N1=n;
}
Matrix::Matrix(){

}
int Matrix::getN(){   //Vozvrashaet kol-vo tochek
    return N1;
}
void Matrix::setN(int n){  //redactiruet kol-vo tochek v matrix grafa
    N1=n;
}
Min_tree::Min_tree(){
    N_points=0;
    N_cluster=0;
}
int Min_tree::Print_tree(std::ofstream& outdata){  //pechataet tree
    int i,k,j, mk,mj;
    for(k=0;k<m_Field->getN();k++){
        for(j=0;j<int(m_Matrix_graf.Mat[k].mass.size());j++){
                for(i=0;i<m_Field->getN();i++){
                    if(m_Field->m1_Point[i].getN()==k){mk=i;}
                    if(m_Field->m1_Point[i].getN()==m_Matrix_graf.Mat[k][j]){mj=i;}
                }
                outdata<<m_Field->m1_Point[mk].getx1()<<" "<<m_Field->m1_Point[mk].getx2()<<" "<<0<<std::endl;
                outdata<<m_Field->m1_Point[mj].getx1()<<" "<<m_Field->m1_Point[mj].getx2()<<" "<<0<<std::endl;
                outdata<<std::endl;
        }
    }
    return 0;
}
int Min_tree::Start_alg(){   //poisk tree
    std::vector<int> param;
    int N=0,i,j,k, mini, minj, maxi, maxj, flag=0;
    float minel, maxel;
    N_points=m_Field->getN();
    m_Matrix_graf.setN(m_Field->getN());
    for(i=0;i<m_Field->getN();i++)
    {
        m_Matrix_graf.Mat.push_back(arr());
    }
    for(i=0;i<m_Field->getN();i++)
    {
        param.push_back(0);
    }
    minel=m_Field->m_Matrix_dist.Mat_dist[1];
    maxel=m_Field->m_Matrix_dist.Mat_dist[1];
    mini=0;
    minj=1;
    maxi=0;
    maxj=1;
    for(i=0;i<m_Field->getN();i++)
    {
        for(j=i+1;j<m_Field->getN();j++)
        {
            if(m_Field->m_Matrix_dist.Mat_dist[i*m_Field->m_Matrix_dist.getN()+j]<minel)
            {
                minel=m_Field->m_Matrix_dist.Mat_dist[i*m_Field->m_Matrix_dist.getN()+j];
                mini=i;
                minj=j;
            }
            if(m_Field->m_Matrix_dist.Mat_dist[i*m_Field->m_Matrix_dist.getN()+j]>maxel)
            {
                maxel=m_Field->m_Matrix_dist.Mat_dist[i*m_Field->m_Matrix_dist.getN()+j];
                maxi=i;
                maxj=j;
            }
        }
    }
    param[minj]=1;
    param[mini]=1;
    if(LogOn==true)
        {log1<<N<<" "<<mini<<" "<<minj<<" "<<minel<<std::endl;}
    m_Matrix_graf.Mat[mini]+=minj;
    N++;
    m_Matrix_graf.Mat[mini].m[m_Matrix_graf.Mat[mini].m.size()-1]=minel;
    m_Matrix_graf.Mat[minj]+=mini;
    m_Matrix_graf.Mat[minj].m[m_Matrix_graf.Mat[minj].m.size()-1]=minel;
    while(1){
        minel=maxel;
        mini=maxi;
        minj=maxj;
        flag=0;
        for(k=0;k<m_Field->getN();k++){
            if(param[k]==1){
                for(j=0;j<m_Field->getN();j++){
                    if((m_Field->m_Matrix_dist.Mat_dist[k*m_Field->m_Matrix_dist.getN()+j]<minel)&&(k!=j)&&(param[j]!=1)){
                        minel=m_Field->m_Matrix_dist.Mat_dist[k*m_Field->m_Matrix_dist.getN()+j];
                        mini=k;
                        minj=j;
                        flag=1;
                    }
                }
            }
        }
        if(flag==0){break;}
        param[minj]=1;
        param[mini]=1;
        if(LogOn==true)
            {log1<<N<<" "<<mini<<" "<<minj<<" "<<minel<<std::endl;}
        m_Matrix_graf.Mat[mini]+=minj;
        N++;
        m_Matrix_graf.Mat[mini].m[m_Matrix_graf.Mat[mini].m.size()-1]=minel;
        m_Matrix_graf.Mat[minj]+=mini;
        m_Matrix_graf.Mat[minj].m[m_Matrix_graf.Mat[minj].m.size()-1]=minel;
        if(flag==0){break;}

    }
    return 0;
}
int Min_tree::getN(){   //Vozvrashaet kol-vo tochek
    return N_points;
}
int Min_tree::getCl(){   //Vozvrashaet kol-vo clusterov
    return N_cluster;
}
int Min_tree::SetP(int N1){ //redactiruet kol-vo tochek
    N_points=N1;
    return 0;
}
int Min_tree::ClusterS(){  //poisk clusterov po tree
    int i,j,k,t,flag=0;
    std::vector<int> param;
    for(i=0;i<m_Field->getN();i++){
        param.push_back(0);
    }
    i=0;
    for(k=0;k<m_Field->getN();k++){
        for(j=0;j<int(m_Matrix_graf.Mat[k].m.size());j++){
                if(m_Matrix_graf.Mat[k].m[j]>m_Histogramm.getQ()){
                    i++;
                    if(j!=int(m_Matrix_graf.Mat[k].mass.size()-1)){
                            m_Matrix_graf.Mat[k].mass[j]=m_Matrix_graf.Mat[k].mass[m_Matrix_graf.Mat[k].getN()-1];
                            m_Matrix_graf.Mat[k].setN(m_Matrix_graf.Mat[k].getN()-1);
                    }
                    else{
                        m_Matrix_graf.Mat[k].setN(m_Matrix_graf.Mat[k].getN()-1);
                    }
                }
        }
    }
    t=1;
    for(i=0;i<m_Field->getN();i++){
        if(param[i]==0){
            m_Cluster.push_back(Cluster());
            N_cluster++;
            t=1;
            param[i]=t;
            while(1){
                flag=0;
                for(k=0;k<m_Field->getN();k++){
                    if(param[k]==t){
                        for(j=0;j<int(m_Matrix_graf.Mat[k].getN());j++){

                            if((param[m_Matrix_graf.Mat[k][j]]==0)){
                                param[m_Matrix_graf.Mat[k][j]]=t+1;
                                flag=1;
                            }
                        }
                    }
                }
                t=t+1;
                if(flag==0){break;}
            }
            for(k=0;k<m_Field->getN();k++){
                int m,j;
                if((param[k]!=0)&&(param[k]!=-1)){
                    for(j=0;j<m_Field->getN();j++){
                        if(m_Field->m1_Point[j].getN()==k){m=j;}
                    }
                    m_Cluster[N_cluster-1].m_Point.push_back(&m_Field->m1_Point[m]);
                    m_Cluster[N_cluster-1].SetP(m_Cluster[N_cluster-1].getN()+1);
                    N_points++;
                    param[k]=-1;
                }
            }
            //m_Cluster[N_cluster-1].searchekstr();
        }
    }
    return 0;
}
Histogramm::Histogramm(){
    NQ=0;
    recQ=0;
}
float Histogramm::getQ(){ //Vozvrashaet recomenduemoe q
    return recQ;
}
int Histogramm::Print_hist(std::ofstream& outdata){  //pechataet histogramm
    for(int i=0;i<NQ;i++){
        outdata<<otr[i]<<":"<<otr[i+1]<<" "<<num[i]<<std::endl;
    }
    return 0;
}
int Histogramm::Start(int Nq){   //poisk histogramm
    int i,j;
    float r;
    qmax=q[0];
    qmin=q[0];
    NQ=Nq;
    for(i=0;i<int(q.size());i++){
        if(q[i]<qmin){qmin=q[i];}
        if(q[i]>qmax){qmax=q[i];}
    }
    otr.push_back(0);
    for(i=0;i<Nq;i++){
        num.push_back(0);
        otr.push_back(0);
    }
    r=(qmax-qmin)/Nq;
    otr[0]=qmin;
    for(i=1;i<=Nq;i++){
        otr[i]=qmin+i*r;
    }
    for(i=0;i<int(q.size());i++){
        for(j=0;j<Nq;j++){
            if(q[i]<otr[j+1]){
                    num[j]++;
                    break;
            }
        }
    }
    for(i=0;i<Nq;i++){
        if(num[i]==0){
            std::cout<<"Recomendation d:"<<(otr[i]+otr[i+1])/2<<std::endl;
            recQ=(otr[i]+otr[i+1])/2;
            break;
        }
    }
    return 0;
}
int Wave_alg::SetCluster(int N1){  //redactiruet kol-vo clusterov v cloud
    N_cluster=N1;
    return 0;
}
int Wave_alg::getN(){    //Vozvrashaet kol-vo tochek
    return N_points;
}
int Wave_alg::SetP(int N1){  //redactiruet kol-vo tochek
    N_points=N1;
    return 0;
}
int Wave_alg::getCl(){  //Vozvrashaet kol-vo clusterov
    return N_cluster;
}
int Wave_alg::savematrix_graf(float d){  //sozdaet matrix grafa s parametrom d
    if(LogOn==true)
        {log1<<"Start poisk matrix graf"<<std::endl;}
    int i,j;
    m_Matrix_graf.setN(m_Field->getN());
    m_Matrix_graf.Mat.clear();
    for(i=0;i<m_Field->getN();i++)
    {
        m_Matrix_graf.Mat.push_back(arr());
    }
    for(i=0;i<m_Field->getN();i++)
    {
        for(j=i+1;j<m_Field->getN();j++)
        {
            if((m_Field->m_Matrix_dist.Mat_dist[i*m_Field->m_Matrix_dist.getN()+j]<d)&&(i!=j))
            {
                m_Matrix_graf.Mat[i]+=j;
                m_Matrix_graf.Mat[j]+=i;
            }
        }
    }
    if(LogOn==true){
        for(i=0;i<m_Field->getN();i++)
        {
            for(j=0;j<int(m_Matrix_graf.Mat[i].getN());j++)
            {
                log1<<m_Matrix_graf.Mat[i][j]<<" ";
            }
            log1<<'\n';
        }
    }
    if(LogOn==true)
        {log1<<"Finish poisk matrix graf"<<std::endl;}
	return 0;
}
Wave_alg::Wave_alg(){
    N_cluster=0;
    N_points=0;
}
int Wave_alg::Start_alg(float d){ // zapusk volnovogo algoritma
    if(LogOn==true)
        {log1<<"Start wave"<<std::endl;}
    int i,j,k,t,flag=0;
    t=1;
    for(i=0;i<m_Field->getN();i++){
        if(param[i]==0){
            m_Cluster.push_back(Cluster());
            N_cluster++;
            t=1;
            param[i]=t;
            while(1){
                flag=0;
                for(k=0;k<m_Field->getN();k++){
                    if(param[k]==t){
                        for(j=0;j<m_Matrix_graf.Mat[k].getN();j++){

                            if((param[m_Matrix_graf.Mat[k].mass[j]]==0)){
                                param[m_Matrix_graf.Mat[k].mass[j]]=t+1;
                                flag=1;
                            }
                        }
                    }
                }
                t=t+1;
                if(flag==0){break;}
            }
            for(k=0;k<m_Field->getN();k++){
                int m,j;
                if((param[k]!=0)&&(param[k]!=-1)){
                    for(j=0;j<m_Field->getN();j++){
                        if(m_Field->m1_Point[j].getN()==k){m=j;break;}
                    }
                    m_Cluster[N_cluster-1].m_Point.push_back(&m_Field->m1_Point[m]);
                    m_Cluster[N_cluster-1].SetP(m_Cluster[N_cluster-1].getN()+1);
                    N_points++;
                    param[k]=-1;
                }
            }
            //m_Cluster[N_cluster-1].searchekstr();
        }
    }
    if(LogOn==true)
        {log1<<"Finish wave"<<std::endl;}
	return 0;
}
DBScan::DBScan(){
    N_cluster=0;
    N_points=0;
}
int DBScan::Start_alg(float d, int k){ //start algoritm DBScan
    int i,j;
    if(LogOn==true)
        {log1<<"Start dbscan"<<std::endl;}
    m1_Wave_alg.savematrix_graf(d);
    for(i=0;i<m1_Wave_alg.m_Field->getN();i++){
        if(int(m1_Wave_alg.m_Matrix_graf.Mat[i].mass.size())>=k){
            m1_Wave_alg.param[i]=1;
        }
    }
    for(i=0;i<m1_Wave_alg.m_Field->getN();i++){
        for(j=0;j<int(m1_Wave_alg.m_Matrix_graf.Mat[i].getN());j++){
            if(m1_Wave_alg.param[i]==1){m1_Wave_alg.param[m1_Wave_alg.m_Matrix_graf.Mat[i][j]]=0;}
        }
        if(m1_Wave_alg.param[i]==1){m1_Wave_alg.param[i]=0;}
    }
    m1_Wave_alg.Start_alg(d);
    N_cluster=m1_Wave_alg.getCl();
    N_points=m1_Wave_alg.getN();
    if(LogOn==true)
        {log1<<"Finish dbscan"<<std::endl;}
    return 0;
}
int DBScan::getN(){  //vozvrashaet kol-vo tochek
    return N_points;
}
int DBScan::getCl(){  //vozvrashaet kol-vo clusterov
    return N_cluster;
}
int DBScan::SetCluster(int N1){ //menyaet kol-vo clusterov
    N_cluster=N1;
    return 0;
}
int DBScan::SetP(int N1){  //menyaet kol-vo tochek
    N_points=N1;
    return 0;
}
int k_means::creat_matr(int k){
    int i,j;
    m_Matrix.setN(m_Field->getN());
    m_Matrix.Mat.clear();
    for(i=0;i<k;i++)
    {
        m_Matrix.Mat.push_back(arr());
    }
    for(i=0;i<k;i++)
    {
        for(j=0;j<m_Field->getN();j++)
        {
            if(i==0){
                m_Matrix.Mat[i]+=1;
            }
            else{m_Matrix.Mat[i]+=0;}
        }
    }
    return 0;
}
k_means::k_means(){
    N_points=0;
    N_cluster=0;
}
int k_means::Start_alg(int k){ //zapusk algoritma k srednih
    int i,j=0,m,m1,c,num;
    float rast,minp,c1,c2,raz=0;
    c=m_Field->getN()/k;
    for(i=0;i<k;i++)
    {
            j=(rand() % c)+i*c;
            for(m=0;m<m_Field->getN();m++)
            {
                if(m_Field->m1_Point[m].getN()==j){
                    centrs.push_back(m_Field->m1_Point[m].getx1());
                    centrs.push_back(m_Field->m1_Point[m].getx2());
                }
            }
    }
    while(1){
        raz=0;
        for(m=0;m<m_Field->getN();m++){
            minp=sqrt((m_Field->m1_Point[m].getx1()-centrs[0])*(m_Field->m1_Point[m].getx1()-centrs[0])+(m_Field->m1_Point[m].getx2()-centrs[1])*(m_Field->m1_Point[m].getx2()-centrs[1]));
            num=0;
            for(i=0;i<k;i++)
            {
                rast=sqrt((m_Field->m1_Point[m].getx1()-centrs[2*i])*(m_Field->m1_Point[m].getx1()-centrs[2*i])+(m_Field->m1_Point[m].getx2()-centrs[2*i+1])*(m_Field->m1_Point[m].getx2()-centrs[2*i+1]));
                if(rast<minp){
                    minp=rast;
                    num=i;
                }
            }
            for(i=0;i<k;i++)
            {
                if((m_Matrix.Mat[i].mass[m]==0)&&(num==i)){
                    m_Matrix.Mat[i].mass[m]=1;
                }
                if((m_Matrix.Mat[i].mass[m]==1)&&(num!=i)){
                    m_Matrix.Mat[i].mass[m]=0;
                }
            }
        }
        for(i=0;i<k;i++)
        {
            j=0;
            c1=0;
            c2=0;
            for(m=0;m<m_Field->getN();m++){
                if(m_Matrix.Mat[i].mass[m]==1){
                    j++;
                    for(m1=0;m1<m_Field->getN();m1++)
                    {
                        if(m_Field->m1_Point[m1].getN()==m){
                            c1=c1+m_Field->m1_Point[m1].getx1();
                            c2=c2+m_Field->m1_Point[m1].getx2();
                        }
                    }
                }
            }
            if(j!=0){
                raz=raz+fabs(centrs[2*i]-c1/j)+fabs(centrs[2*i+1]-c2/j);
                centrs[2*i]=c1/j;
                centrs[2*i+1]=c2/j;
            }
             if(LogOn==true)
                {log1<<i<<" "<<centrs[2*i]<<" "<<centrs[2*i+1]<<" "<<raz<<std::endl;}
        }
        if(raz<eps){break;}
    }
    for(i=0;i<k;i++)
    {
        m_Cluster.push_back(Cluster());
        N_cluster++;
        for(m=0;m<m_Field->getN();m++){
                if(m_Matrix.Mat[i].mass[m]==1){
                    for(m1=0;m1<m_Field->getN();m1++)
                    {
                        if(m_Field->m1_Point[m1].getN()==m){
                            m_Cluster[N_cluster-1].m_Point.push_back(&m_Field->m1_Point[m1]);
                            m_Cluster[N_cluster-1].SetP(m_Cluster[N_cluster-1].getN()+1);
                            N_points++;
                        }
                    }
                }
        }
        m_Cluster[N_cluster-1].searchekstr();
        m_Cluster[N_cluster-1].O.newpoint(centrs[2*i],centrs[2*i+1]);
    }
    return 0;
}
int k_means::getN(){return N_points;} //Vozvrashaet kol-vo tochek
int k_means::getCl(){return N_cluster;} //Vozvrashaet kol-vo clusterov
int k_means::SetCluster(int N1){N_cluster=N1;return 0;} //redactiruet kol-vo clusterov
int k_means::SetP(int N1){N_points=N1;return 0;}  //redactiruet kol-vo tochek
int EM::creat_matr(int k){
    int i,j;
    m_Matrix.setN(m_Field->getN());
    m_Matrix.Mat.clear();
    for(i=0;i<k;i++)
    {
        m_Matrix.Mat.push_back(arr());
    }
    for(i=0;i<k;i++)
    {
        for(j=0;j<m_Field->getN();j++)
        {
            m_Matrix.Mat[i].m.push_back(0);
        }
    }
    return 0;
}
EM::EM(){
    N_points=0;
    N_cluster=0;
}
int EM::Start_alg(int k){ //zapusk EM algoritma
    int i,j=0,m,m1,m2,c, num;
    std::vector<int> color;
    float raz=0, sym=0,corx1, corx2, drob;
    std::string txt=".txt";
    std::string multik ="multik";
    std::ofstream mult;
    c=m_Field->getN()/k;
    k1=k;
    for(i=0;i<k;i++)
    {
            j=(rand() % c)+i*c;
            for(m=0;m<m_Field->getN();m++)
            {
                if(m_Field->m1_Point[m].getN()==j){
                    centrs.push_back(m_Field->m1_Point[m].getx1());
                    centrs.push_back(m_Field->m1_Point[m].getx2());
                }
            }
    }
    drob=1.0/k;
    for(i=0;i<k;i++)
    {
        massa.push_back(drob);
    }
    for(i=0;i<k;i++)
    {
        sigma.push_back(1.0);
        sigma.push_back(0.0);
        sigma.push_back(0.0);
        sigma.push_back(1.0);
    }
    raz=100000.0;
    for(m=0;m<k;m++){
            j=(rand() % 128)+(rand() % 128);
            color.push_back(j);
            j=(rand() % 128)+(rand() % 128);
            color.push_back(j);
            j=(rand() % 128)+(rand() % 128);
            color.push_back(j);
    }
    j=0;
    while(1){
        std::string fail;
        for(i=0;i<m_Field->getN();i++)
        {
            sym=0.0;
            for(m=0;m<k;m++){
                sym=sym+massa[m]*fi(i,m);
            }
            for(m=0;m<k;m++){
                m_Matrix.Mat[m].m[i]=(massa[m]*fi(i,m))/sym;
            }
        }
        for(m=0;m<k;m++){
            massa[m]=0.0;
            centrs[2*m]=0.0;
            centrs[2*m+1]=0.0;
            sigma[m*4]=0.0;
            sigma[m*4+1]=0.0;
            sigma[m*4+2]=0.0;
            sigma[m*4+3]=0.0;
            for(i=0;i<m_Field->getN();i++){
                for(m1=0;m1<m_Field->getN();m1++)
                {
                    if(m_Field->m1_Point[m1].getN()==i){
                        break;
                    }
                }
                massa[m]=massa[m]+m_Matrix.Mat[m].m[i];
                centrs[2*m]=centrs[2*m]+m_Matrix.Mat[m].m[i]*m_Field->m1_Point[m1].getx1();
                centrs[2*m+1]=centrs[2*m+1]+m_Matrix.Mat[m].m[i]*m_Field->m1_Point[m1].getx2();
            }
            massa[m]=massa[m]/m_Field->getN();
            centrs[2*m]=centrs[2*m]/(m_Field->getN()*massa[m]);
            centrs[2*m+1]=centrs[2*m+1]/(m_Field->getN()*massa[m]);
        }
        for(m=0;m<k;m++){
            sigma[m*4]=0.0;
            sigma[m*4+1]=0.0;
            sigma[m*4+2]=0.0;
            sigma[m*4+3]=0.0;
            for(i=0;i<m_Field->getN();i++){
                for(m1=0;m1<m_Field->getN();m1++)
                {
                    if(m_Field->m1_Point[m1].getN()==i){
                        break;
                    }
                }
                corx1=m_Field->m1_Point[m1].getx1()-centrs[2*m];
                corx2=m_Field->m1_Point[m1].getx2()-centrs[2*m+1];
                sigma[m*4]=sigma[m*4]+corx1*corx1*m_Matrix.Mat[m].m[i];
                sigma[m*4+1]=sigma[m*4+1]+corx2*corx1*m_Matrix.Mat[m].m[i];
                sigma[m*4+2]=sigma[m*4+2]+corx1*corx2*m_Matrix.Mat[m].m[i];
                sigma[m*4+3]=sigma[m*4+3]+corx2*corx2*m_Matrix.Mat[m].m[i];
            }
            sigma[m*4]=sigma[m*4]/(m_Field->getN()*massa[m]);
            sigma[m*4+1]=sigma[m*4+1]/(m_Field->getN()*massa[m]);
            sigma[m*4+2]=sigma[m*4+2]/(m_Field->getN()*massa[m]);
            sigma[m*4+3]=sigma[m*4+3]/(m_Field->getN()*massa[m]);
        }
        /*fail=multik+std::__cxx11::to_string(j)+txt;
        mult.open(fail, std::ofstream::out | std::ofstream::trunc);
        for(i=0;i<k;i++)
        {
            for(m=0;m<m_Field->getN();m++){
                num=0;
                for(m2=1;m2<k;m2++){
                    if(m_Matrix.Mat[m2].m[m]>m_Matrix.Mat[num].m[m]){num=m2;}
                }
                if(num==i){
                    for(m1=0;m1<m_Field->getN();m1++)
                    {
                        if(m_Field->m1_Point[m1].getN()==m){
                            mult<<m_Field->m1_Point[m1].getx1()<<" "<<m_Field->m1_Point[m1].getx2()<<" "<<color[3*i]<<" "<<color[3*i+1]<<" "<<color[3*i+2]<<std::endl;
                        }
                    }
                }
            }
            mult<<centrs[2*i]<<" "<<centrs[2*i+1]<<" "<<0<<" "<<0<<" "<<0 <<std::endl;
        }
        mult.close();*/
        j=j+1;
        if(fabs(raz-prov())<eps){break;}
        else{raz=prov();}
    }
    if(LogOn==true)
        {log1<<"index "<<j<<std::endl;}
    for(i=0;i<k;i++)
    {
        m_Cluster.push_back(Cluster());
        N_cluster++;
        for(m=0;m<m_Field->getN();m++){
                num=0;
                for(m2=1;m2<k;m2++){
                    if(m_Matrix.Mat[m2].m[m]>m_Matrix.Mat[num].m[m]){num=m2;}
                }
                if(num==i){
                    for(m1=0;m1<m_Field->getN();m1++)
                    {
                        if(m_Field->m1_Point[m1].getN()==m){
                            m_Cluster[N_cluster-1].m_Point.push_back(&m_Field->m1_Point[m1]);
                            m_Cluster[N_cluster-1].SetP(m_Cluster[N_cluster-1].getN()+1);
                            N_points++;
                        }
                    }
                }
        }
        m_Cluster[N_cluster-1].searchekstr();
        m_Cluster[N_cluster-1].O.newpoint(centrs[2*i],centrs[2*i+1]);
    }
    //razl();

    for(m=0;m<k;m++){
        float dist=0.0;
        for(m1=0;m1<m_Field->getN();m1++){
            num=0;
            for(m2=1;m2<k;m2++){
                if(m_Matrix.Mat[m2].m[m1]>m_Matrix.Mat[num].m[m1]){num=m2;}
            }
            if((num==m)&&(delt(m1,m)>dist)){dist=delt(m1,m);}
        }
        if(LogOn==true)
            {log1<<"f"<<m+1<<"(x,y) = "<<sigma[m*4+3]/deter(m)<<"*(x-("<<centrs[2*m]<<"))"<<"*(x-("<<centrs[2*m]<<"))+("<<-2.0*sigma[m*4+1]/deter(m)<<")*(x-("<<centrs[2*m]<<"))"<<"*(y-("<<centrs[2*m+1]<<"))+("<<sigma[m*4]/deter(m)<<")*(y-("<<centrs[2*m+1]<<"))*(y-("<<centrs[2*m+1]<<"))-("<<dist*dist<<")"<<std::endl;}
    }

    return 0;
}
/*
int EM::razl(){
    int i,j,m;
    std::vector<float> d1;
    std::vector<float> d2;
    float **E1 = new float *[m_Field->getN()];

    for (int i = 0; i < m_Field->getN(); i++){
        E1[i] = new float [m_Field->getN()];
        d1.push_back(0);
    }
    float **E2 = new float *[k1];

    for (int i = 0; i < k1; i++){
        E2[i] = new float [k1];
        d2.push_back(0);
    }
    float **v1 = new float *[m_Field->getN()];

    for (int i = 0; i < m_Field->getN(); i++)
        v1[i] = new float [m_Field->getN()];
    float **v2 = new float *[k1];

    for (int i = 0; i < k1; i++)
        v2[i] = new float [k1];
    for(i=0;i<m_Field->getN();i++){
            for(j=0;j<m_Field->getN();j++){
                E1[i][j]=0.0;
                for(m=0;m<k1;m++){
                    E1[i][j]=E1[i][j]+m_Matrix.Mat[m].m[i]*m_Matrix.Mat[m].m[j];
                }
            }
    }
    for(i=0;i<k1;i++){
            for(j=0;j<k1;j++){
                E2[i][j]=0.0;
                for(m=0;m<m_Field->getN();m++){
                    E2[i][j]=E2[i][j]+m_Matrix.Mat[i].m[m]*m_Matrix.Mat[j].m[m];
                }
            }
    }
//    jacobi(m_Field->getN(),E1,d1,v1);
//    inversion(v1,m_Field->getN());
//    jacobi(k1,E2,d2,v2);
//    inversion(v2,k1);
    for(i=0;i<k1;i++){
        for(j=0;j<k1;j++){
            log1<<v2[j][i]<<" ";
        }
        log1<<std::endl;
    }
    for(j=0;j<m_Field->getN();j++){
        log1<<d1[j]<<std::endl;
    }
    for(i=0;i<=k1;i++){
        delete v2[i];
        delete E2[i];
    }
    for(i=0;i<m_Field->getN();i++){
        delete v1[i];
        delete E1[i];
    }
    delete [] E1;
    delete [] E2;
    delete [] v1;
    delete [] v2;
    return 0;
}*/
float EM::prov(){ //poisk proverochnogo virageniya
    int m, i;
    float prover=0.0, sym1;
    for(i=0;i<m_Field->getN();i++){
            sym1=0.0;
            for(m=0;m<k1;m++){
                sym1=sym1+massa[m]*fi(i,m);
            }
            prover=prover+log(sym1);
    }
    return prover;
}
float EM::deter(int j){ //poisk determine sigma
    return (sigma[j*4]*sigma[j*4+3]-sigma[j*4+1]*sigma[j*4+2]);
}
float EM::delt(int i, int j){ //poisk rastoyaniya Mahanalobisa
    float cor1, cor2, det;
    int m1;
    for(m1=0;m1<m_Field->getN();m1++)
    {
        if(m_Field->m1_Point[m1].getN()==i){
            break;
        }
    }
    cor1=m_Field->m1_Point[m1].getx1()-centrs[2*j];
    cor2=m_Field->m1_Point[m1].getx2()-centrs[2*j+1];
    det=deter(j);
    return sqrt(cor1*cor1*sigma[j*4+3]/det-cor1*cor2*(sigma[j*4+1]+sigma[j*4+2])/det+cor2*cor2*sigma[j*4]/det);
}
float EM::fi(int i, int j){ //poisk fi
    float del;
    del=delt(i,j);
    return(exp(-del*del/2)/(2*PI*sqrt(deter(j))));
}
int EM::getN(){return N_points;} //Vozvrashaet kol-vo tochek
int EM::getCl(){return N_cluster;} //Vozvrashaet kol-vo clusterov
int EM::SetCluster(int N1){N_cluster=N1;return 0;} //redactiruet kol-vo clusterov
int EM::SetP(int N1){N_points=N1;return 0;}  //redactiruet kol-vo tochek
/*void jacobi ( const unsigned int n,float **a, std::vector<float> d,float **v )
{
    if ( n == 0 ) return;
    float * b = new float[n+n];
    float * z = b + n;
    unsigned int i, j;
    for ( i = 0; i < n; ++i )
    {
        z[i] = 0.;
        b[i] = d[i] = a[i][i];
        for ( j = 0; j < n; ++j )
        {
            v[i][j] = i == j ? 1. : 0.;
        }
    }
    for ( i = 0; i < 200; ++i )
    {
        float sm = 0.;
        unsigned int p, q;
        for ( p = 0; p < n - 1; ++p )
        {
            for ( q = p + 1; q < n; ++q )
            {
                sm += fabs ( a[p][q] );
            }
        }
        if ( sm == 0 ) break;
        const float tresh = i < 3 ? 0.2 * sm / ( n*n ) : 0.;
        for ( p = 0; p < n - 1; ++p )
        {
            for ( q = p + 1; q < n; ++q )
            {
                const float g = 1e12 * fabs ( a[p][q] );
                if ( i >= 3 && fabs ( d[p] ) > g && fabs ( d[q] ) > g ) a[p][q] = 0.;
                else
                if ( fabs ( a[p][q] ) > tresh )
                {
                    const float theta = 0.5 * ( d[q] - d[p] ) / a[p][q];
                    float t = 1. / ( fabs(theta) + sqrt(1.+theta*theta) );
                    if ( theta < 0 ) t = - t;
                    const float c = 1. / sqrt ( 1. + t*t );
                    const float s = t * c;
                    const float tau = s / ( 1. + c );
                    const float h = t * a[p][q];
                    z[p] -= h;
                    z[q] += h;
                    d[p] -= h;
                    d[q] += h;
                    a[p][q] = 0.;
                    for ( j = 0; j < p; ++j )
                    {
                        const float g = a[j][p];
                        const float h = a[j][q];
                        a[j][p] = g - s * ( h + g * tau );
                        a[j][q] = h + s * ( g - h * tau );
                    }
                    for ( j = p+1; j < q; ++j )
                    {
                        const float g = a[p][j];
                        const float h = a[j][q];
                        a[p][j] = g - s * ( h + g * tau );
                        a[j][q] = h + s * ( g - h * tau );
                    }
                    for ( j = q+1; j < n; ++j )
                    {
                        const float g = a[p][j];
                        const float h = a[q][j];
                        a[p][j] = g - s * ( h + g * tau );
                        a[q][j] = h + s * ( g - h * tau );
                    }
                    for ( j = 0; j < n; ++j )
                    {
                        const float g = v[j][p];
                        const float h = v[j][q];
                        v[j][p] = g - s * ( h + g * tau );
                        v[j][q] = h + s * ( g - h * tau );
                    }
                }
            }
        }
        for ( p = 0; p < n; ++p )
        {
            d[p] = ( b[p] += z[p] );
            z[p] = 0.;
        }
    }
    delete[] b;
}
void inversion(float **A, int N)
{
    float temp;

    float **E = new float *[N];

    for (int i = 0; i < N; i++)
        E[i] = new float [N];

    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
        {
            E[i][j] = 0.0;

            if (i == j)
                E[i][j] = 1.0;
        }

    for (int k = 0; k < N; k++)
    {
        temp = A[k][k];

        for (int j = 0; j < N; j++)
        {
            A[k][j] /= temp;
            E[k][j] /= temp;
        }

        for (int i = k + 1; i < N; i++)
        {
            temp = A[i][k];

            for (int j = 0; j < N; j++)
            {
                A[i][j] -= A[k][j] * temp;
                E[i][j] -= E[k][j] * temp;
            }
        }
    }

    for (int k = N - 1; k > 0; k--)
    {
        for (int i = k - 1; i >= 0; i--)
        {
            temp = A[i][k];

            for (int j = 0; j < N; j++)
            {
                A[i][j] -= A[k][j] * temp;
                E[i][j] -= E[k][j] * temp;
            }
        }
    }

    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            A[i][j] = E[i][j];

    for (int i = 0; i < N; i++)
        delete [] E[i];

    delete [] E;
}
*/
