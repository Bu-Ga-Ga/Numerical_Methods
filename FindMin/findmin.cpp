#include <cmath>
#include <iostream>
#include <vector>

static const int N = 3;
static const double EPS = 0.000000001;
static const double EPS1 = sqrt(EPS);
static const double TAU = sqrt(EPS) * 0.1;
	
using std::cout;
using std::endl;
using std::vector;

double f1(vector<double> x);
double f(vector<double> x, vector<double> grad, double alpha);
vector<double> epsmassplus(vector<double> mass, int i, int N);
vector<double> epsmassminus(vector<double> mass, int i, int N);
vector<double> grad(vector<double> mass,int N);
vector<vector<double>> Gess(vector<double> mass,int N);
double Determinant(vector<vector<double>> mas, int m);
vector<vector<double>> GetMatr(vector<vector<double>> mas, int row, int col,int N);
vector<vector<double>> Transpone(vector<vector<double>> mas, int N);
vector<vector<double>> Mreverse(vector<vector< double>> mas, int N);
vector<double> smin(vector<double> sp, int N);
vector<double> gmin(vector<double> sp, int N);




int main(){
	vector<double> x = {1,1,1};
	vector<double> t(N);
	vector<double> t1(N);
	vector<double> t2(N);

    t = smin(x,N);
    for(int i = 0 ; i < N ; i++) cout << t[i] << "  ";
    cout << endl << f(t,t,0) << endl << endl << endl;
    t2 = gmin(t,N);
    for(int i = 0 ; i < N ; i++) cout << t2[i] << "  ";
    cout << endl << f(t2,t2,0) << endl << endl << endl;

	}


double f(vector<double> x, vector<double> grad, double alpha){
	vector<double> tmp(N);
	for(int i = 0; i < N ; i++){tmp[i] = (x[i] - (alpha * grad[i]));
	}
	return (2 * pow(tmp[0]+5,2) + pow(tmp[1]-10, 2) + 10 * pow(x[2],2));
}







//-------------
vector<double> epsmassplus(vector<double> mass, int i, int N){
    vector<double> mass1(N);
	mass1[i - 1] = mass[i - 1] + TAU;
	for(int j = 0; j < N; j++){
		if(j != (i - 1)) mass1[j] = mass[j];
	}
    return mass1;
}
//------------
vector<double> epsmassminus(vector<double> mass, int i, int N){
    vector<double> mass1(N);
	mass1[i - 1] = mass[i - 1] - TAU;
	for(int j = 0; j < N; j++){
		if(j != (i - 1)) mass1[j] = mass[j];
	}
    return mass1;
}
//------------
vector<double> grad(vector<double> mass,int N){
	
    vector<double> grad1(N);
	for(int i = 0; i < N; i++){
		vector<double> ap(N);
		vector<double> am(N);
		ap = epsmassplus(mass, i + 1, N);
		am = epsmassminus(mass, i + 1, N);
		grad1[i] = (f(ap,ap,0) - f(am,am,0))/(2 * TAU);
	}
    return grad1;
}
//-----------
vector<vector<double>> Gess(vector<double> mass,int N){
     vector<vector<double>> gess(N, vector<double> (N));

	for(int i = 0; i < N; i++){
		for(int j = 0; j < N; j++){
			vector<double> ap(N);
			vector<double> am(N);
           			ap = epsmassplus(mass, i + 1, N);
			am = epsmassminus(mass, i + 1, N);
			if(i == j){
				gess[i][j] =( f(ap,ap,0) - 2 * f(mass,mass,0) + f(am,am,0) ) / pow(TAU, 2);
		//		cout << f(ap) <<" "<<f(mass) <<" "<<f(am) << endl;
			}

			if(i != j){
			//	double sap[N];
			//	double sam[N];
				vector<double> comp(N);
				vector<double> comp1(N);
				vector<double> amj(N);
				comp = epsmassplus(mass, i + 1, N);
				comp1 = epsmassminus(comp, j + 1, N);
				amj = epsmassminus(mass, j + 1, N);	
			//	epsmassplus(mass, sap, j + 1, 2);
			//	epsmassminus(mass, sam, j + 1, 2);
	
				gess[i][j] = (( f(ap,ap,0) - f(mass,mass,0) ) - (f(comp1,comp1,0) - f(amj,amj,0))) /(TAU * TAU);

		//		cout << f(ap) <<" "<<f(mass) <<" "<<f(comp1) <<" "<<f(amj) << endl;
			}	
		}
	}
    return gess;
}
//------------
double Determinant(vector<vector<double>> mas, int m) {
  int k;
  vector<vector<double>> p(N, vector<double> (N));
  double d = 0;
  k = 1; //(-1) в степени i
  if (m < 1) {cout <<" ERROR "; return 0; }
  if (m == 1) { d = mas[0][0]; return(d); }
  if (m == 2) { d = mas[0][0] * mas[1][1] - (mas[1][0] * mas[0][1]); return(d); }
  if (m > 2) {
    for (int i = 0; i < m; i++) {
      p = GetMatr(mas, i, 0,m);
      d = d + k * mas[i][0] * Determinant(p, m - 1);
      k = -k;
    }
  }
  return(d);
}
//------------
vector<vector<double>> GetMatr(vector<vector<double>> mas, int row, int col,int N) {
  int di, dj;
  di = 0;
  vector<vector<double>> mas1(N, vector<double> (N));
  for (int i = 0; i < N - 1; i++) { 
    if (i == row)  
      di = 1;   
    dj = 0;
    for (int j = 0; j < N - 1; j++) { 
      if (j == col) 
        dj = 1; 
      mas1[i][j] = mas[i + di][j + dj];
      
    }
  }
  return mas1;
}
//-------------
vector<vector<double>> Transpone(vector<vector<double>> mas, int N) {
  vector<vector<double>> mas1(N, vector<double>(N));
  for (int i = 0; i < N; i++){
    for (int j = 0; j < N; j++)
      mas1[i][j] = mas[j][i];
  }
  return mas1;
}
//-------------
vector<vector<double>> Mreverse(vector<vector< double>> mas, int N) {  
  vector<vector<double>> rez(N, vector<double>(N));
  double det = Determinant(mas, N); // находим определитель исходной матрицы
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      rez[i][j] = Determinant(GetMatr(mas,i, j, N), N - 1);
      if ((i + j) % 2 == 1)       // если сумма индексов строки и столбца нечетная
        rez[i][j] = -rez[i][j];    // меняем знак минора
      rez[i][j] = rez[i][j] / det;
    }
  }
  return Transpone(rez,N);   
}
//------------
vector<double> smin(vector<double> sp, int N){
    vector<double> res(N);
    vector<double> Grad(N);
    Grad = grad(sp,N);
    vector<double> pointa(N);
    vector<double> pointb(N);
    vector<double> point1(N);
    double lambda = 0.5;
    double mu = 3;
    double alpha = 0.01;
    double beta = 1;
    int flag = 0;
    int c = 0;
    int d = 0;
    while(1){
        double vec[N];
        for(int i = 0; i < N; i++){ pointa[i] = sp[i] - alpha * Grad[i];}
        if(f(pointa,pointa,0) <= f(sp,sp,0) ) { break; } 
        if(f(pointa,pointa,0) > f(sp,sp,0) ) { alpha = alpha * lambda; }
    }
    for(int i = 0; i < N; i++){pointa[i] = sp[i];}
    for(int i = 0; i < N; i++){pointb[i] = sp[i] - alpha * Grad[i]; }
    for(int i = 0; i < N; i++) {res[i] = pointb[i];}
    for(int i = 0; i < N; i++) point1[i] = pointa[i] - pointb[i];
    double a = 0;
    for(int i = 0; i < N; i++) a += pow(point1[i],2);
    double sqr = sqrt(a);
    Grad = grad(pointa, N);
    int gr = 0;
    for(int i = 0; i < N; i++) gr += pow(Grad[i],2);
    double grsqr = sqrt(gr);

    if(std::abs( f(pointb,pointb,0) - f(pointa,pointa,0)) > EPS1 || sqr > EPS1 || gr > EPS1) pointb = smin(pointb,N);
    return pointb;
}
//-----------
vector<double> gmin(vector<double> sp, int N){
	vector<double> res(N);
    vector<double> Grad(N);
	vector<vector<double>> gess(N, vector<double>(N));
	vector<vector<double>> rgess(N, vector<double>(N));
	vector<double> h(N);
	double cn;
	double cnn = 0;
	double x1 = 2;
	double x2 = 0;
	int flag = 0;
	gess = Gess(sp, N);
	Grad = grad(sp, N);

	rgess = Mreverse(gess,N);
	for(int i = 0; i < N; i++){
		for(int j = 0; j < N; j++){
			h[i] +=( rgess[i][j] * Grad[i]);
		}
	}
	while(1){
		cn = ( f(sp,h,x2) - f(sp, h,x1) +  (((f(sp,h,x1+TAU)-f(sp,h,x1-TAU) )/(2*TAU)) * x1) -  (((f(sp,h,x2+TAU)-f(sp,h,x2-TAU) )/(2*TAU)) * (x2))  ) /( (((f(sp,h,x1+TAU)-f(sp,h,x1-TAU) )/(2*TAU))) - ((f(sp,h,x2+TAU)-f(sp,h,x2-TAU) )/(2*TAU)));
    //    cout << cn << " "<< cnn <<endl;

		if(((f(sp,h,cn+TAU)-f(sp,h,cn-TAU) )/(2*TAU)) <= 0) {cnn = x2; x2 = cn;flag++;}
		if(((f(sp,h,cn+TAU)-f(sp,h,cn-TAU) )/(2*TAU)) > 0) {cnn = x1; x1 = cn;flag++;}
      //  cout << "^"<<std::abs(((f(sp,h,cn+TAU)-f(sp,h,cn-TAU) )/(2*TAU))) << "^"<<endl;
	if(std::abs(((f(sp,h,cn+TAU)-f(sp,h,cn-TAU) )/(2*TAU))) < EPS && std::abs(cn-cnn)<EPS1) { break;}
	
	}
//	cout << "!!!!!!!!!!!!"<<endl<<cn<<endl<<cnn<<endl;
    vector<double> tres(N);
	vector<double> res1(N);
	for(int i = 0; i < N; i++) res[i] = sp[i] - cn * h[i];

	for(int i = 0; i < N; i++) {res1[i] = sp[i] - cn * h[i];/*cout <<"<" <<res1[i]<<">  ";*/}

	for(int i = 0; i < N; i++) tres[i] = sp[i] - cnn * h[i];

	vector<double> point1(N);
	for(int i = 0; i < N; i++) point1[i] = res1[i] - tres[i];
	double a = 0;
	for(int i = 0; i < N; i++) a += pow(point1[i],2);
	double sqr = sqrt(a);

	Grad = grad(res1,N);
	int gr = 0;
	for(int i = 0; i < N; i++) gr += pow(Grad[i],2);
	double grsqr = sqrt(gr);


	if(std::abs( f(sp,h,cn) - f(sp,h,cnn)) > EPS || sqr > EPS || gr > EPS) res = gmin(res1, N);

	return res;
}


//--------------
#if 0
double f1(vector<double> x){
	return (3 * pow(x[0]+1,2) + pow(x[1],2) + 10 * pow(x[2],2));
}
#endif
//--------------























































































































