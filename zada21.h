//211-Ziatdinov-Nail
const int LIMIT = 1200;
extern std::ofstream log1;
extern bool LogOn;
class Field;
class Controller;
class arr
{
private:
    int N;
public:
    std::vector<int> mass;
    std::vector<float> m;
    arr();
    int getN();
    int setN(int N1);
    int operator[](const int k);
    void operator+=(const int k);
};
class Matrix
{

public:
	Matrix();

	int getN();
	void setN(int n);
    std::vector<arr> Mat;
private:
    int N1;

};
class Matrix_dist
{

public:
	Matrix_dist();
	int getN();
	void setN(int n);
    std::vector<float> Mat_dist;
private:
    int N1;
};
class Point
{

public:
	Point();
	Point(float x1, float x2, float y, int N, int Cl1);
	float getx1();
	float getx2();
	float gety();
	int getN();
	void setNum(int N);
	int newpoint(float ox, float oy);
	Point & operator=(Point arg);
	int Cl;
private:
	float X1;
	float X2;
	float Y;
	int Num;
};
class Cloud
{
public:
	std::vector<Point*> m_Point;
    Cloud();
	Cloud(float ox, float oy,float xdisper, float ydisper,int from, int num);
	Point createpoint(float ox, float oy,float xdisper, float ydisper,int from, int Cl);
	int Print_Cloud(std::ofstream& outdata,int CL);
	int SetN(int N1);
	int getN();
	void getekstr(float &max1,float &max2,float &min1,float &min2);
	int searchekstr();
    Point O;
    float xdisp;
    float ydisp;
private:
    float maxx1;
    float maxx2;
    float minx1;
    float minx2;
	int N_points;
};
class Buffer
{
public:
	Buffer();
	Point O;
	int getN();
	std::vector<Point> m2_Point;
	int rotateCl(float ang);
	int shift(float ox, float oy);
	int zoom(float ly);
	int copyCl(Cloud* n_Cloud);
	int upload();
private:
    int N;
    bool Status_Work;
};
class Field
{
public:
    std::vector<Cloud> m1_Cloud;
    std::vector<Point> m1_Point;
    Matrix_dist m_Matrix_dist;
    Field();
    int createcloud(float ox, float oy,float xdisper, float ydisper, int num);
	int Print_Field(std::ofstream& outdata);
	int SetP(int N1);
	int SetCloud(int N1);
	int getCl();
	int getN();
	int work();
    int savematrix_dist();
    float rastpoint(int n1, int n2);
private:
	int N_clouds;
	int N_points;
	bool Status_Work;
};
class Cluster
{

public:
	Cluster();
	std::vector<Point*> m_Point;
	int Print_Cluster(std::ofstream& outdata);
    int SetP(int N1);
    int getN();
    void getekstr(float &max1,float &max2,float &min1,float &min2);
	int searchekstr();
	Point O;
private:
    float maxx1;
    float maxx2;
    float minx1;
    float minx2;
	int N_points;

};

class Wave_alg
{

public:
    Field* m_Field;
    std::vector<Cluster> m_Cluster;
    Matrix m_Matrix_graf;
    std::vector<int> param;
    int savematrix_graf(float d);
	Wave_alg();
	int Start_alg(float d);
    int getN();
    int getCl();
    int SetCluster(int N1);
    int SetP(int N1);
private:
	int N_cluster;
	int N_points;

};
class k_means
{
public:
    Field* m_Field;
    std::vector<Cluster> m_Cluster;
    Matrix m_Matrix;
    std::vector<float> centrs;
    int creat_matr(int k);
	k_means();
	int Start_alg(int k);
    int getN();
    int getCl();
    int SetCluster(int N1);
    int SetP(int N1);
private:
	int N_cluster;
	int N_points;

};
class EM
{
public:
    Field* m_Field;
    std::vector<Cluster> m_Cluster;
    Matrix m_Matrix;
    std::vector<float> centrs;
    std::vector<float> sigma;
    std::vector<float> massa;
    int creat_matr(int k);
	EM();
	int Start_alg(int k);
	float deter(int j);
	float delt(int i, int j);
	float fi(int i, int j);
	float prov();
	int razl();
    int getN();
    int getCl();
    int SetCluster(int N1);
    int SetP(int N1);
private:
    int k1;
	int N_cluster;
	int N_points;

};
class DBScan
{

public:
    Wave_alg m1_Wave_alg;
	DBScan();
	int Start_alg(float d, int k);
    int getN();
    int getCl();
    int SetCluster(int N1);
    int SetP(int N1);
private:
	int N_cluster;
	int N_points;

};
class Histogramm
{

public:
	Histogramm();
	std::vector<float> q;
	int Start(int Nq);
	int Print_hist(std::ofstream& outdata);
	float getQ();
private:
    std::vector<float> otr;
    std::vector<int> num;
	float qmax;
	float qmin;
	int NQ;
	float recQ;
};
class Min_tree
{

public:
    Field* m_Field;
    Histogramm m_Histogramm;
    Matrix m_Matrix_graf;
    std::vector<Cluster> m_Cluster;
	Min_tree();
	int Print_tree(std::ofstream& outdata);
	int Start_alg();
	int ClusterS();
    int getN();
    int getCl();
    int SetP(int N1);
private:
	int N_points;
	int N_cluster;

};
class Launch
{

public:
    Histogramm* m1_Histogramm;
    Launch();
	Launch(float delt, int k1);
	int Start_alg(int n);
    std::vector<Cluster*> m1_Cluster;
    Wave_alg m_Wave_alg;
    DBScan m_DBScan;
    Min_tree m_Min_tree;
    k_means m_k_means;
    EM m_EM;
    int Print_Cluster(std::ofstream& outdata);
    int Print_tree(std::ofstream& outdata);
    int getN();
    int getCl();
    int SetCluster(int N1);
    int SetP(int N1);
    float d;
	int k;
private:
	int N_cluster;
	int N_points;

};
class Controller
{

public:
	Field m_Field;
	Buffer m_Buffer;
	std::vector<Launch> m_Launch;
	Controller();
	~Controller();
	int Create_Cloud(int numcl, float ox, float oy, float xdisp, float ydisp, int num);
	int Print_Cloud(int Ncl,std::ofstream& outdata);
	int Print_Cluster(int Ncl,std::ofstream& outdata);
	int Print_Field(std::ofstream& outdata);
	int Print_tree(int Ncl,std::ofstream& outdata);
	int Wave(float d);
	int k_means(int k);
	int EM(int k);
	int DBScan(float d,int k);
	int Min_tree();
	int hist(int num, int k);
	int Print_hist(int Ncl,std::ofstream& outdata);

private:
	bool Status_Work;
};
class Interface
{
public:
	Interface();
	~Interface(){}
	Controller m_Controller;
	int Starts();
	bool command(std::string strcomm, std::string config, std::ofstream& logfile);
};
//void inversion(float **A, int N);
//void jacobi ( const unsigned int n,float **a, std::vector<float> d,float **v );