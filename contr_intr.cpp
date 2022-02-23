#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <vector>
#include "zada21.h"
//211-Ziatdinov-Nail
//File s function
Controller::Controller(){
    if(LogOn==true)
        {log1<<"Contoller sozdan"<<std::endl;}
    Status_Work=true;
}
Controller::~Controller(){
    m_Launch.clear();
}
int Controller::Create_Cloud(int numcl, float ox, float oy,float xdisp, float ydisp, int num){ //sozdaet cloud s centrom ox i oy, kolvo tochek - num
    if(LogOn==true)
        {log1<<"Start create cloud"<<std::endl;}
    if(m_Field.getCl()==(numcl-1)){
        if(m_Field.work()==1){
            m_Field.m1_Cloud.push_back(Cloud(ox,oy,xdisp,ydisp,m_Field.getN(),num));
            //m_Field.m1_Cloud.push_back(Cloud());
            m_Field.m1_Cloud[numcl-1].xdisp=xdisp;
            m_Field.m1_Cloud[numcl-1].ydisp=ydisp;
            m_Field.createcloud(ox,oy,xdisp,ydisp,num);
            m_Field.m1_Cloud[numcl-1].O.newpoint(ox,oy);
            m_Field.SetCloud(1);
        }
        else{
            std::cout << "Error,don't add Cloud, because you saved matrix dist." << std::endl;
            if(LogOn==true)
                {log1<<"Error with adding cloud "<<numcl<<std::endl;}
            return -1;
        }
    }
    else{
        std::cout << "Error,don't understand number of Cloud." << std::endl;
        return -1;
    }
    if(LogOn==true)
        {log1<<"Finish create cloud"<<std::endl;}
	return 0;
}
int Controller::Print_Cloud(int Ncl,std::ofstream& outdata){ //pechataet cloud
    m_Field.m1_Cloud[Ncl-1].Print_Cloud(outdata,Ncl-1);
	return 0;
}
int Controller::Print_Cluster(int Ncl,std::ofstream& outdata){ //pechataet cloud
    m_Launch[Ncl-1].Print_Cluster(outdata);
	return 0;
}
int Controller::Print_Field(std::ofstream& outdata){ //pechataet field
    m_Field.Print_Field(outdata);
	return 0;
}
int Controller::Print_tree(int Ncl,std::ofstream& outdata){
    m_Launch[Ncl-1].Print_tree(outdata);
    return 0;
}
int Controller::Print_hist(int Ncl,std::ofstream& outdata){
    m_Launch[Ncl-1].m1_Histogramm->Print_hist(outdata);
    return 0;
}
int Controller::Wave(float d){ //sozdaet volnovoy algoritm c parametrom d
    int n;
    n = m_Launch.size();
    m_Launch.push_back(Launch(d,1));
    m_Launch[n].m_Wave_alg.m_Field=&m_Field;
    m_Launch[n].Start_alg(1);
	return 0;
}
int Controller::DBScan(float d,int k){ //sozdaet DBscan c parametrami d i k
    int n;
    n = m_Launch.size();
    m_Launch.push_back(Launch(d,k));
    m_Launch[n].m_DBScan.m1_Wave_alg.m_Field=&m_Field;
    m_Launch[n].Start_alg(2);
	return 0;
}
int Controller::Min_tree(){ //sozdaet Min tree
    int n;
    n = m_Launch.size();
    m_Launch.push_back(Launch(0,1));
    m_Launch[n].m_Min_tree.m_Field=&m_Field;
    m_Launch[n].Start_alg(3);
	return 0;
}
int Controller::hist(int num, int k){
    int i,j;
    for(i=0;i<(m_Field.getN()/2);i++){
        for(j=0;j<int(m_Launch[num-1].m_Min_tree.m_Matrix_graf.Mat[i].m.size());j++){
            m_Launch[num-1].m_Min_tree.m_Histogramm.q.push_back(m_Launch[num-1].m_Min_tree.m_Matrix_graf.Mat[i].m[j]);
        }
    }
    m_Launch[num-1].m_Min_tree.m_Histogramm.Start(k);
    m_Launch[num-1].m1_Histogramm=&m_Launch[num-1].m_Min_tree.m_Histogramm;
    m_Launch[num-1].Start_alg(4);
    return 0;
}
int Controller::k_means(int k){ //sozdaet k means c parametrom k
    int n;
    m_Launch.push_back(Launch(0.0,k));
    n = m_Launch.size();
    m_Launch[n-1].m_k_means.m_Field=&m_Field;
    m_Launch[n-1].Start_alg(5);
    return 0;
}
int Controller::EM(int k){ //sozdaet EM c parametrom k
    int n;
    n = m_Launch.size();
    m_Launch.push_back(Launch(0.0,k));
    m_Launch[n].m_EM.m_Field=&m_Field;
    m_Launch[n].Start_alg(6);
    return 0;
}
Interface::Interface(){
}
int Interface::Starts(){ //function interface
    bool k=true;
    int j;
    std::string config;
    std::string logfile;
    std::ifstream conf("config.txt");
    while(conf){
        std::string comm;
        getline(conf, comm);
        config=config+comm;
    }
    std::ofstream ofs;
    std::string file;
    int m=config.find("Cloudtxt");
    m=m+8+1;
    while(config[m]!=';')
    {
        file.push_back(config[m]);
        m=m+1;
    }
    ofs.open(file, std::ofstream::out | std::ofstream::trunc);
    ofs.close();
    file.clear();
    m=config.find("Clustertxt");
    m=m+10+1;
    while(config[m]!=';')
    {
        file.push_back(config[m]);
        m=m+1;
    }
    ofs.open("cluster.txt", std::ofstream::out | std::ofstream::trunc);
    ofs.close();
    m=config.find("LogRec");
    m=m+6+1;
    if(config[m]=='1'){
        if(LogOn==false){
            LogOn=true;
        }
    }
    else{
        LogOn=false;
    }
    std::ofstream log2;
    j=config.find("LogInter");
    j=j+8+1;
    while(config[j]!=';')
    {
        logfile.push_back(config[j]);
        j=j+1;
    }
    log2.open(logfile);
	while (k)
	{
		std::string strcomm;
		getline(std::cin, strcomm);
		k=command(strcomm,config,log2);
	}
	log2.close();
	return 0;
}
bool Interface::command(std::string strcomm, std::string config, std::ofstream& logfile) {//function raspoznavaniya command
        int i,j, k, m;
        std::string comm;
		j=1;
		if(strcomm[strcomm.length()-1]!=';')
        {
            strcomm.push_back(';');
        }
        for(i=0; i<int(strcomm.length());i++)
        {
            if(strcomm[i]==';')
            {
                j=0;
                i=i+1;
                k=config.find(comm);
                k=k+int(comm.length())+2;
                break;
            }
            if(strcomm[i]==':')
            {
                j=1;
                i=i+1;
                break;
            }
            comm.push_back(strcomm[i]);
        }
        if(LogOn==true)
            {logfile<<"Schitana commanda:"<<comm<<std::endl;}
        if((comm=="PRINT_FIELD")||(comm=="print_field"))
        {
            std::ofstream outdata;
            std::string fail;
            if(j==0){
                std::string file;
                int m=config.find("Cloudtxt");
                m=m+8+1;
                while(config[m]!=';')
                {
                    file.push_back(config[m]);
                    m=m+1;
                }
                outdata.open(file, std::ios_base::app);
            }
            else{
                while(i<int(strcomm.length())-1)
                {
                    fail.push_back(strcomm[i]);
                    i=i+1;
                }
                std::cout << fail << std::endl;
                outdata.open(fail, std::ios_base::app);
            }
            m_Controller.Print_Field(outdata);
            outdata.close();
            return true;
        }
        if((comm=="PRINT_CLOUD")||(comm=="print_cloud"))
        {
            std::string ncl;
            std::string fail;
            std::ofstream outdata;
            while(i<int(strcomm.length()))
            {
                if(strcomm[i]=='_')
                {
                    j=2;
                    i=i+1;
                    break;
                }
                ncl.push_back(strcomm[i]);
                i=i+1;
            }
            while(i<int(strcomm.length())-1)
            {
                fail.push_back(strcomm[i]);
                i=i+1;
            }
            if(j==2){
                outdata.open(fail, std::ofstream::out | std::ofstream::trunc);
                outdata.close();
                outdata.open(fail, std::ios_base::app);
            }
            else{
                std::string file;
                int m=config.find("Cloudtxt");
                m=m+8+1;
                while(config[m]!=';')
                {
                    file.push_back(config[m]);
                    m=m+1;
                }
                outdata.open(file, std::ios_base::app);
            }
            m_Controller.Print_Cloud(stoi(ncl), outdata);
            outdata.close();

            return true;
        }
        if((comm=="RUN")||(comm=="run"))
        {
            std::string fail;
            std::ifstream inf;
            if(j==0){
                std::string file;
                int m=config.find("Commandtxt");
                m=m+10+1;
                while(config[m]!=';')
                {
                    file.push_back(config[m]);
                    m=m+1;
                }
                inf.open(file);
            }
            else{
                while(i<int(strcomm.length())-1)
                {
                    fail.push_back(strcomm[i]);
                    i=i+1;
                }
                std::cout << fail << std::endl;
                inf.open(fail);
            }
            while(getline(inf, strcomm)){
                command(strcomm,config,logfile);
            }
            inf.close();

            return true;
        }
        if((comm=="BUFFER_ADD_CLOUD")||(comm=="buffer_add_cloud"))
        {
            std::string ncl;
            while(i<int(strcomm.length()))
            {
                ncl.push_back(strcomm[i]);
                i=i+1;
            }
            if(LogOn==true)
                {logfile<<"Schitani parametri:"<<stoi(ncl)<<std::endl;}
            for(m=0; m<m_Controller.m_Field.getN();m++){
                if(m_Controller.m_Field.m1_Point[m].Cl==(stoi(ncl)-1))
                m_Controller.m_Buffer.m2_Point.push_back(m_Controller.m_Field.m1_Point[m]);
            }
            m_Controller.m_Buffer.copyCl(&m_Controller.m_Field.m1_Cloud[int(stoi(ncl)-1)]);
            return true;
        }
         if((comm=="BUFFER_SHIFT_CLOUD")||(comm=="buffer_shift_cloud"))
        {
            std::string xcor;
            std::string ycor;
            while(i<int(strcomm.length()))
            {
                if(strcomm[i]=='_')
                {
                    i=i+1;
                    break;
                }
                xcor.push_back(strcomm[i]);
                i=i+1;
            }
            while(i<int(strcomm.length()))
            {
                ycor.push_back(strcomm[i]);
                i=i+1;
            }
            if(LogOn==true)
                {logfile<<"Schitani parametri:"<<stod(xcor)<<" "<<stod(ycor)<<std::endl;}
            m_Controller.m_Buffer.shift(std::stod(xcor),std::stod(ycor));
            return true;
        }
        if((comm=="BUFFER_TURN_CLOUD")||(comm=="buffer_turn_cloud"))
        {
            std::string ang;
            while(i<int(strcomm.length()))
            {
                ang.push_back(strcomm[i]);
                i=i+1;
            }
            if(LogOn==true)
                {logfile<<"Schitani parametri:"<<stod(ang)<<std::endl;}
            m_Controller.m_Buffer.rotateCl(std::stod(ang));
            return true;
        }
            if((comm=="BUFFER_COMPRESSION_CLOUD")||(comm=="buffer_compression_cloud"))
        {
            std::string ly;
            while(i<int(strcomm.length()))
            {
                ly.push_back(strcomm[i]);
                i=i+1;
            }
            if(LogOn==true)
                {logfile<<"Schitani parametri:"<<stod(ly)<<std::endl;}
            m_Controller.m_Buffer.zoom(std::stod(ly));
            return true;
        }
        if((comm=="BUFFER_UPLOAD_CLOUD")||(comm=="buffer_upload_cloud"))
        {
            std::string ncl;
            while(i<int(strcomm.length()))
            {
                ncl.push_back(strcomm[i]);
                i=i+1;
            }
            if(LogOn==true)
                {logfile<<"Schitani parametri:"<<stoi(ncl)<<std::endl;}
            if(m_Controller.m_Buffer.m2_Point.size()<1){return false;}
            if(LogOn==true)
                {log1<<"Start bufer upload"<<std::endl;}
            m_Controller.m_Field.m1_Cloud.push_back(Cloud());
            for(i=0;i<m_Controller.m_Buffer.getN();i++){
                m_Controller.m_Buffer.m2_Point[i].setNum(m_Controller.m_Field.getN());
                m_Controller.m_Buffer.m2_Point[i].Cl=stoi(ncl)-1;
                m_Controller.m_Field.m1_Point.push_back(m_Controller.m_Buffer.m2_Point[i]);
                m_Controller.m_Field.SetP(1);
                m_Controller.m_Field.m1_Cloud[stoi(ncl)-1].m_Point.push_back(&m_Controller.m_Field.m1_Point[m_Controller.m_Field.m1_Point.size()-1]);
            }
            m_Controller.m_Field.m1_Cloud[stoi(ncl)-1].SetN(m_Controller.m_Buffer.getN());
            m_Controller.m_Field.m1_Cloud[stoi(ncl)-1].O=m_Controller.m_Buffer.O;
            m_Controller.m_Field.m1_Cloud[stoi(ncl)-1].searchekstr();
            m_Controller.m_Field.SetCloud(1);
            m_Controller.m_Buffer.upload();
            return true;
        }
        if((comm=="GENERATE_CLOUD")||(comm=="generate_cloud"))
        {
            std::string numcl;
            std::string x1cor;
            std::string x2cor;
            std::string x1disp;
            std::string x2disp;
            std::string ncor;
            if(j==1){
                while(i<int(strcomm.length()))
                {
                    if(strcomm[i]=='_')
                    {
                        i=i+1;
                        break;
                    }
                    numcl.push_back(strcomm[i]);
                    i=i+1;
                }
                while(i<int(strcomm.length()))
                {
                    if(strcomm[i]=='_')
                    {
                        i=i+1;
                        break;
                    }
                    x1cor.push_back(strcomm[i]);
                    i=i+1;
                }
                 while(i<int(strcomm.length()))
                {
                    if(strcomm[i]=='_')
                    {
                        i=i+1;
                        break;
                    }
                    x2cor.push_back(strcomm[i]);
                    i=i+1;
                }
                while(i<int(strcomm.length()))
                {
                    if(strcomm[i]=='_')
                    {
                        i=i+1;
                        break;
                    }
                    x1disp.push_back(strcomm[i]);
                    i=i+1;
                }
                 while(i<int(strcomm.length()))
                {
                    if(strcomm[i]=='_')
                    {
                        i=i+1;
                        break;
                    }
                    x2disp.push_back(strcomm[i]);
                    i=i+1;
                }
                 while(i<int(strcomm.length()))
                {
                    ncor.push_back(strcomm[i]);
                    i=i+1;
                }
            }
            else{
                while(config[k]!=';')
                {
                    if(config[k]=='_')
                    {
                        k=k+1;
                        break;
                    }
                    numcl.push_back(config[k]);
                    k=k+1;
                }
                while(config[k]!=';')
                {
                    if(config[k]=='_')
                    {
                        k=k+1;
                        break;
                    }
                    x1cor.push_back(config[k]);
                    k=k+1;
                }
                 while(config[k]!=';')
                {
                    if(config[k]=='_')
                    {
                        k=k+1;
                        break;
                    }
                    x2cor.push_back(config[k]);
                    k=k+1;
                }
                while(config[k]!=';')
                {
                    if(config[k]=='_')
                    {
                        k=k+1;
                        break;
                    }
                    x1disp.push_back(config[k]);
                    k=k+1;
                }
                 while(config[k]!=';')
                {
                    if(config[k]=='_')
                    {
                        k=k+1;
                        break;
                    }
                    x2disp.push_back(config[k]);
                    k=k+1;
                }
                 while(config[k]!=';')
                {
                    ncor.push_back(config[k]);
                    k=k+1;
                }
            }
            if(LogOn==true)
                {logfile<<"Schitani parametri:"<<stoi(numcl)<<" "<<stod(x1cor)<<" "<<stod(x2cor)<<" "<<stod(x1disp)<<" "<<stod(x2disp)<<" "<<stoi(ncor)<<std::endl;}
            m_Controller.Create_Cloud(std::stoi(numcl),std::stod(x1cor),std::stod(x2cor),std::stod(x1disp),std::stod(x2disp),std::stoi(ncor));
            return true;

        }
        if((comm=="WAVE")||(comm=="wave"))
        {
            std::string d;
            while(i<int(strcomm.length()))
            {
                d.push_back(strcomm[i]);
                i=i+1;
            }
            if(LogOn==true)
                {logfile<<"Schitani parametri:"<<stod(d)<<std::endl;}
            m_Controller.m_Field.savematrix_dist();
            m_Controller.Wave(std::stod(d));
            return true;
        }
        if((comm=="KMEANS")||(comm=="kmeans"))
        {
            std::string k;
            while(i<int(strcomm.length()))
            {
                k.push_back(strcomm[i]);
                i=i+1;
            }
            if(LogOn==true)
                {logfile<<"Schitani parametri:"<<stoi(k)<<std::endl;}
            m_Controller.m_Field.savematrix_dist();
            m_Controller.k_means(std::stoi(k));
            return true;
        }
        if((comm=="EM")||(comm=="em"))
        {
            std::string k;
            while(i<int(strcomm.length()))
            {
                k.push_back(strcomm[i]);
                i=i+1;
            }
            if(LogOn==true)
                {logfile<<"Schitani parametri:"<<stoi(k)<<std::endl;}
            m_Controller.m_Field.savematrix_dist();
            m_Controller.EM(std::stoi(k));
            return true;
        }
        if((comm=="HISTOGRAMM")||(comm=="histogramm"))
        {
            std::string d;
            std::string k;
            while(i<int(strcomm.length()))
            {
                if(strcomm[i]=='_')
                {
                    i=i+1;
                    break;
                }
                d.push_back(strcomm[i]);
                i=i+1;
            }
            while(i<int(strcomm.length()))
            {
                k.push_back(strcomm[i]);
                i=i+1;
            }
            if(LogOn==true)
                {logfile<<"Schitani parametri:"<<stoi(d)<<" "<<stoi(k)<<std::endl;}
            m_Controller.hist(stoi(d),stoi(k));
            return true;
        }
        if((comm=="PRINT_HISTOGRAMM")||(comm=="print_histogramm"))
        {
            std::string ncl;
            std::string fail;
            std::ofstream outdata;
            while(i<int(strcomm.length()))
            {
                if(strcomm[i]=='_')
                {
                    j=2;
                    i=i+1;
                    break;
                }
                ncl.push_back(strcomm[i]);
                i=i+1;
            }
            while(i<int(strcomm.length())-1)
            {
                fail.push_back(strcomm[i]);
                i=i+1;
            }
            if(j==2){
                outdata.open(fail, std::ofstream::out | std::ofstream::trunc);
                outdata.close();
                outdata.open(fail, std::ios_base::app);
            }
            else{
                std::string file;
                int m=config.find("Histtxt");
                m=m+8+1;
                while(config[m]!=';')
                {
                    file.push_back(config[m]);
                    m=m+1;
                }
                outdata.open(file, std::ios_base::app);
            }
            m_Controller.Print_hist(stoi(ncl), outdata);
            outdata.close();

            return true;
        }
        if((comm=="MIN_TREE")||(comm=="min_tree"))
        {

            m_Controller.m_Field.savematrix_dist();
            m_Controller.Min_tree();
            return true;
        }
        if((comm=="PRINT_TREE")||(comm=="print_tree"))
        {
            std::ofstream outdata;
            std::string ncl,fail;
            if(j==0){
                std::string file;
                int m=config.find("treetxt");
                m=m+7+1;
                while(config[m]!=';')
                {
                    if(config[m]=='_')
                    {
                        m=m+1;
                        break;
                    }
                    ncl.push_back(config[m]);
                    m=m+1;
                }
                while(config[m]!=';')
                {
                    file.push_back(config[m]);
                    m=m+1;
                }
                outdata.open(file, std::ofstream::out | std::ofstream::trunc);
                outdata.close();
                outdata.open(file, std::ios_base::app);
            }
            else{
                while(i<int(strcomm.length()))
                {
                    if(strcomm[i]=='_')
                    {
                        i=i+1;
                        break;
                    }
                    ncl.push_back(strcomm[i]);
                    i=i+1;
                }
                while(i<int(strcomm.length())-1)
                {
                    fail.push_back(strcomm[i]);
                    i=i+1;
                }
                outdata.open(fail, std::ofstream::out | std::ofstream::trunc);
                outdata.close();
                outdata.open(fail, std::ios_base::app);
            }
            m_Controller.Print_tree(stoi(ncl),outdata);
            outdata.close();
            return true;
        }
        if((comm=="DBSCAN")||(comm=="dbscan"))
        {
            std::string d;
            std::string k;
            while(i<int(strcomm.length()))
            {
                if(strcomm[i]=='_')
                {
                    i=i+1;
                    break;
                }
                d.push_back(strcomm[i]);
                i=i+1;
            }
            while(i<int(strcomm.length()))
            {
                k.push_back(strcomm[i]);
                i=i+1;
            }
            if(LogOn==true)
                {logfile<<"Schitani parametri:"<<stod(d)<<" "<<stoi(k)<<std::endl;}
            m_Controller.m_Field.savematrix_dist();
            m_Controller.DBScan(stod(d),stoi(k));
            return true;
        }
        if((comm=="PRINT_CLUSTER")||(comm=="print_cluster"))
        {
            std::string ncl;
            std::string fail;
            std::ofstream outdata;
            while(i<int(strcomm.length()))
            {
                if(strcomm[i]=='_')
                {
                    j=2;
                    i=i+1;
                    break;
                }
                ncl.push_back(strcomm[i]);
                i=i+1;
            }
            while(i<int(strcomm.length())-1)
            {
                fail.push_back(strcomm[i]);
                i=i+1;
            }
            if(j==2){
                outdata.open(fail, std::ofstream::out | std::ofstream::trunc);
                outdata.close();
                outdata.open(fail, std::ios_base::app);
            }
            else{
                std::string file;
                int m=config.find("Clustertxt");
                m=m+10+1;
                while(config[m]!=';')
                {
                    file.push_back(config[m]);
                    m=m+1;
                }
                outdata.open(file, std::ios_base::app);
            }
            m_Controller.Print_Cluster(stoi(ncl), outdata);
            outdata.close();

            return true;
        }
        if((comm=="help")||(comm=="HELP"))
        {
            std::string line;
            std::ifstream in("help.txt");
            if (in)
            {
                while (getline(in, line))
                {
                    std::cout << line << std::endl;
                }
            }
            in.close();
            return true;
        }
        if((comm=="info")||(comm=="INFO"))
        {
            std::cout <<"Kol-vo tochek: "<< m_Controller.m_Field.getN() << std::endl;
            std::cout <<"Kol-vo cloud: "<< m_Controller.m_Field.getCl() << std::endl;
            return true;
        }
        if((comm=="info_all")||(comm=="INFO_ALL"))
        {
            float min1,min2,max1,max2;
            for(int i=0;i<m_Controller.m_Field.getCl();i++){
                std::cout <<"Number cloud: "<< i+1 << std::endl;
                std::cout <<"Kol-vo tochek: "<< m_Controller.m_Field.m1_Cloud[i].getN() << std::endl;
                std::cout <<"Centr: "<< m_Controller.m_Field.m1_Cloud[i].O.getx1()<<" " <<m_Controller.m_Field.m1_Cloud[i].O.getx2()<< std::endl;
                m_Controller.m_Field.m1_Cloud[i].getekstr(max1,max2,min1,min2);
                std::cout <<"Max po x1: "<< max1<<" Max po x2: "<< max2 <<" Min po x1: "<< min1<<" Min po x2: "<< min2<< std::endl;
            }
            return true;
        }
        if((comm=="info_cluster")||(comm=="INFO_CLUSTER"))
        {
            std::string numcl;
            while(i<int(strcomm.length()))
            {
                numcl.push_back(strcomm[i]);
                i=i+1;
            }
            float min1,min2,max1,max2;
            for(int i=0;i<m_Controller.m_Launch[stoi(numcl)-1].getCl();i++){
                std::cout <<"Number cluster: "<< i+1 << std::endl;
                std::cout <<"Kol-vo tochek: "<< m_Controller.m_Launch[stoi(numcl)-1].m1_Cluster[i]->getN() << std::endl;
                std::cout <<"Centr: "<< m_Controller.m_Launch[stoi(numcl)-1].m1_Cluster[i]->O.getx1()<<" " <<m_Controller.m_Launch[stoi(numcl)-1].m1_Cluster[i]->O.getx2()<< std::endl;
                m_Controller.m_Launch[stoi(numcl)-1].m1_Cluster[i]->getekstr(max1,max2,min1,min2);
                std::cout <<"Max po x1: "<< max1<<" Max po x2: "<< max2 <<" Min po x1: "<< min1<<" Min po x2: "<< min2<< std::endl;
            }
            return true;
        }
        if((comm=="exit")||(comm=="EXIT"))
        {
            return false;
        }

        else
        {
            std::cout << comm << std::endl;
            std::cout << "Not identified command" << std::endl;
            return true;
        }
}
