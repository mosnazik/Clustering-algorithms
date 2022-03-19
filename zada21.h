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
class Matrix //klass khraneniya matritsy lyubogo razmera
{

public:
	Matrix();

	int getN();
	void setN(int n);
    std::vector<arr> Mat;
private:
    int N1;

};
class Matrix_dist  //klass khraneniya matritsy rasstoyaniy
{

public:
	Matrix_dist();
	int getN();
	void setN(int n);
    std::vector<float> Mat_dist;
private:
    int N1;
};
class Point //klass tochki khranit dve koordinaty i znacheniye. nomer tochki. mozhet ikh izmenyat
{

public:
	Point();
	Point(float x1, float x2, float y, int N, int Cl1);
	float getx1();
	float getx2();
	float gety();
	int getN();
	void setNum(int N);
	int newpoint(float ox, float oy); //sozdanie nowoy tochki
	Point & operator=(Point arg);
	int Cl;
private:
	float X1;
	float X2;
	float Y;
	int Num;
};
class Cloud
//klass khraneniya mnozhestva ssylok tochek v pole. mozhet pechatat oblako iskat ekstremumy. sozdavat tochku. khranit
//ekstremumy tsentr i dispersii. esli ne iz bufera
{
public:
	std::vector<Point*> m_Point;
    Cloud();
	Cloud(float ox, float oy,float xdisper, float ydisper,int from, int num);
	Point createpoint(float ox, float oy,float xdisper, float ydisper,int from, int Cl);
	int Print_Cloud(std::ofstream& outdata,int CL); //pechat' cloud
	int SetN(int N1);
	int getN();
	void getekstr(float &max1,float &max2,float &min1,float &min2);  //poisk ekstremuma
	int searchekstr();  //poisk ekstremuma
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
//klass dlya izmeneniya oblaka. zagruzhayet oblako v vide mnozhestva tochek. izmenyayet ego cherez sdvig. povorot.
//masshtabirovaniye. vygruzhayet v novoye oblako
{
public:
	Buffer();
	Point O;
	int getN();
	std::vector<Point> m2_Point;
	int rotateCl(float ang);         //povorot
	int shift(float ox, float oy);   //sdvig
	int zoom(float ly);              //mashtab
	int copyCl(Cloud* n_Cloud);
	int upload();
private:
    int N;
    bool Status_Work;
};
class Field
//klass v kotorom khranitsya mnozhestvo tochek. gruppirovka ikh po oblakam. mozhet raspechatyvat oblako ili pole.
//dobavlyat oblako. sokhranyat matritsu rasstoyaniy. nakhodit rasstoyaniye mezhdu tochkami
{
public:
    std::vector<Cloud> m1_Cloud;
    std::vector<Point> m1_Point;
    Matrix_dist m_Matrix_dist;
    Field();
    int createcloud(float ox, float oy,float xdisper, float ydisper, int num); //sozdanie oblaka
	int Print_Field(std::ofstream& outdata); //pechat' polya
	int SetP(int N1);
	int SetCloud(int N1);
	int getCl();
	int getN();
	int work();
    int savematrix_dist(); //poisk matrix rastoyaniy
    float rastpoint(int n1, int n2); //poisk rastoyanitya ot tochki n1 do n2
private:
	int N_clouds;
	int N_points;
	bool Status_Work;
};
class Cluster // klass khranyashchiy ssylki na mnozhestvo tochek klastera. takzhe ego tsentr. ekstremumy. kol-vo tochek. mozhet raspechatat tochki
{

public:
	Cluster();
	std::vector<Point*> m_Point;
	int Print_Cluster(std::ofstream& outdata); //pechat' clustera
    int SetP(int N1);
    int getN();
    void getekstr(float &max1,float &max2,float &min1,float &min2); //vozvrashenie ekstremuma
	int searchekstr(); //poisk ekstremuma
	Point O; //centr cluster
private:
    float maxx1;
    float maxx2;
    float minx1;
    float minx2;
	int N_points;

};

class Wave_alg //vypolnyayet volnovoy algoritm
{

public:
    Field* m_Field;
    std::vector<Cluster> m_Cluster;
    Matrix m_Matrix_graf;
    std::vector<int> param;
    int savematrix_graf(float d); //poisk matrix grafa
	Wave_alg();
	int Start_alg(float d); //start volnovogo algoritma
    int getN();
    int getCl();
    int SetCluster(int N1);
    int SetP(int N1);
private:
	int N_cluster;
	int N_points;

};
class k_means  // vypolnyayet algorit k srednikh
{
public:
    Field* m_Field;
    std::vector<Cluster> m_Cluster;
    Matrix m_Matrix;
    std::vector<float> centrs;
    int creat_matr(int k); //sozdanie matrix indicatorov
	k_means();
	int Start_alg(int k);  //start k srednih
    int getN();
    int getCl();
    int SetCluster(int N1);
    int SetP(int N1);
private:
	int N_cluster;
	int N_points;

};
class EM  //vypolnyayet EM algoritm
{
public:
    Field* m_Field;
    std::vector<Cluster> m_Cluster;
    Matrix m_Matrix;
    std::vector<float> centrs;
    std::vector<float> sigma;
    std::vector<float> massa;
    int creat_matr(int k);  //sozdanie matrix indicatorov
	EM();
	int Start_alg(int k);  //start em algoritma
	float deter(int j); //poisk opredelitelya
	float delt(int i, int j); //poisk delta
	float fi(int i, int j); //poisk fi
	float prov(); //proverka algoritma
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
class DBScan //vypolnyayet algoritm DBScan. to est nakhoditsya nuzhnoye mnozhestvo. zatem s nim vypolnyayetsya volnovoy algoritm
{

public:
    Wave_alg m1_Wave_alg;
	DBScan();
	int Start_alg(float d, int k);  //start dbscan
    int getN();
    int getCl();
    int SetCluster(int N1);
    int SetP(int N1);
private:
	int N_cluster;
	int N_points;

};
class Histogramm // klass. stroyashchiy gistogrammu
{

public:
	Histogramm();
	std::vector<float> q;
	int Start(int Nq); //start poisk histogramm
	int Print_hist(std::ofstream& outdata); //pechat' histogramm
	float getQ();
private:
    std::vector<float> otr;
    std::vector<int> num;
	float qmax;
	float qmin;
	int NQ;
	float recQ;
};
class Min_tree  //vypolnyayet algoritm minimalnogo pokryvayushchego dereva. zatem po nemu stroit Histogramm
{

public:
    Field* m_Field;
    Histogramm m_Histogramm;
    Matrix m_Matrix_graf;
    std::vector<Cluster> m_Cluster;
	Min_tree();
	int Print_tree(std::ofstream& outdata);  //pechat' dereva
	int Start_alg();  //start minimalnogo dereva
	int ClusterS();
    int getN();
    int getCl();
    int SetP(int N1);
private:
	int N_points;
	int N_cluster;

};
class Hierarchical  // vypolnyayet algorit Hierarchical
{
public:
    Field* m_Field;
    std::vector<Cluster> m_Cluster;
    Matrix m_Matrix;
    std::vector<int> number;
    std::vector<float> coord;
    int copy_matr(); //copy matrix rast
	Hierarchical();
	int Start_alg(int k);  //start Hierarchical
	float poisk_rast(int k, int n, int m);
	float poisk_cluster(int numk);
	int addpoint(int clus, int num);
    int getN();
    int getCl();
    int SetCluster(int N1);
    int SetP(int N1);
private:
	int N_cluster;
	int N_points;

};
class Launch  //klass zapuskov metodov klassternogo analiza EM.k_means.Min_tree.DBScan.Wave_alg. pechatayet klastery i vydayet po nim informatsiyu
{

public:
    Histogramm* m1_Histogramm;
    Launch();
	Launch(float delt, int k1);
	int Start_alg(int n); //start odnogo iz algoritmov 1 - wave, 2 - dbscan, 3,4 - min derevo, 5 - k means, 6 - EM
    std::vector<Cluster*> m1_Cluster;
    Wave_alg m_Wave_alg;
    DBScan m_DBScan;
    Min_tree m_Min_tree;
    k_means m_k_means;
    Hierarchical m_Hierarchical;
    EM m_EM;
    int Print_Cluster(std::ofstream& outdata);  //pechat' cluster
    int Print_tree(std::ofstream& outdata);  //pechat' dereva
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
class Controller //poluchennyye komandy otpravlyayet libo v Buffer. Field. Launch dlya ikh realizatsii
{

public:
	Field m_Field;
	Buffer m_Buffer;
	std::vector<Launch> m_Launch;
	Controller();
	~Controller();
	int Create_Cloud(int numcl, float ox, float oy, float xdisp, float ydisp, int num); //sozdanie cloud
	int Print_Cloud(int Ncl,std::ofstream& outdata);  //pechat' cloud
	int Print_Cluster(int Ncl,std::ofstream& outdata);  //pechat' cluster
	int Print_Field(std::ofstream& outdata);  //pechat' polya
	int Print_tree(int Ncl,std::ofstream& outdata);  //pechat' dereva
	int Wave(float d);
	int k_means(int k);
	int EM(int k);
	int DBScan(float d,int k);
	int Hierarchical(int k);
	int Min_tree();
	int hist(int num, int k);
	int Print_hist(int Ncl,std::ofstream& outdata);  ////pechat' histogramm

private:
	bool Status_Work;
};
class Interface //klass otvechayet za raspaznovaniye komand i otpravleniye ikh v Controller. takzhe chitayet config help fayly
{
public:
	Interface();
	~Interface(){}
	Controller m_Controller;
	int Starts(); //start
	bool command(std::string strcomm, std::string config, std::ofstream& logfile); //chtenie command
};
//void inversion(float **A, int N);
//void jacobi ( const unsigned int n,float **a, std::vector<float> d,float **v );
