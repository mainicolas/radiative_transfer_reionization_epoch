#include<algorithm> //find function
#include<stdlib.h>
#include<iostream>
#include<fstream> //ofstream
#include<stdio.h>
#include<math.h> //math things
#include<vector> //vector
#include<functional> //trnasform function

#include<omp.h> //parallelism

using namespace std;


double alpha(double T){
	double lambda = (2*157807)/T;
	double res = 1.269e-13;
	res *= pow(lambda, 1.503);
	res /= pow(1 + pow((lambda/0.522),0.47), 1.923);
	return res;
}

double alpha_b(double T){
    double lambda = (2*157807)/T;
    double res    = 2.753e-14;
    res   *= pow(lambda, 1.5);
    res   /= pow(1 + pow((lambda/2.74),0.407), 2.242);
    return res;
}

double beta(double T){
	double lambda = (2*157807)/T;
    double res 	= 21.11*pow(T, (-3/2))*exp(-lambda/2)*pow(lambda, -1.089);
    res /= pow(1 + pow(lambda/0.354, 0.874), 1.01);
    return res;
}

float polynom(double &x, double m, double n, double p, double q){
	return m*x*x*x + n*x*x + p*x + q;
}


float deriv_polynom(double &x, double m, double n, double p){
	return 3*m*x*x + 2*n*x + p;
}


float newton_raphson(double m, double n, double p, double q){
	float eps = 0.0001;
	double x = 0.5;
	double h = polynom(x,m,n,p,q) / deriv_polynom(x,m,n,p);
	while(abs(h) >= eps) 
	{
		h = polynom(x,m,n,p,q) / deriv_polynom(x,m,n,p);
		x -= h;
	}
	return x;
}


float third_order_polynom(double T, double N, int rho, double sigma, double X, double dt){
	int c = 3e8;
	double m = (alpha_b(T) + beta(T)) * (rho*rho) * dt;
	double n = rho - ((alpha(T) + beta(T)) * rho)/(sigma*c) - alpha_b(T) * (rho*rho) * dt - 2 * beta(T) * (rho*rho) * dt;
	double p = -rho * (1 + X) - N - 1/(sigma * c * dt) + (beta(T) * rho)/(sigma * c) + beta(T) * (rho*rho) * dt;
    double q = N + rho * X + X/(sigma * c * dt);

	return newton_raphson(m,n,p,q);
}

void new_advection_equation(double &N, double &F, double T, int rho, double sigma, double X, double x, double dt){
	double c = 3e8;
	N += beta(T)*rho*rho*(1-X)*X*dt - alpha_b(T)*rho*rho*X*X*dt - rho*(X-x);
	F /=  (1 + c*sigma*rho*dt*(1 - X));
}


void glf(int i, int n_cell, vector<vector<double>> &N, vector<vector<double>> &F, vector<vector<double>> &P, vector<double>&glf_f, vector<double>&glf_p){
	int j;
	double c = 3e8;
	
	for(j=0; j<n_cell; j++)
	{
		if (i == 0 && j == 0)
		{
			glf_f.push_back(0);
			glf_p.push_back(0);
		}
		if (i != 0 && j == 0)
		{
			glf_f.push_back( (F[i - 1][n_cell - 1] + F[i][j]) * 0.5 - c/2 * (N[i][j] - N[i - 1][n_cell - 1]) );
			glf_p.push_back( (P[i - 1][n_cell - 1] + P[i][j]) * 0.5 - c/2 * (F[i][j] - F[i - 1][n_cell - 1]) );
		}
		if (j != 0)
		{
			glf_f.push_back((F[i][j - 1] + F[i][j]) * 0.5 - c/2 * (N[i][j] - N[i][j - 1]) );
			glf_p.push_back((P[i][j - 1] + P[i][j]) * 0.5 - c/2 * (F[i][j] - F[i][j - 1]) );
		}
	}
}



void numerical_scheme(double n_photon, vector<vector<double>> &N, vector<vector<double>> &F, vector<vector<double>> &P, vector<vector<double>> &X, vector<vector<double>> &T, vector<int> &source_pos, double dt, double tf, int n_cell, double dx, int rho, double sigma, bool continuous, bool chemistery){
	double c = 3e8;
	vector<double> glf_f, glf_p;
	for(int i = source_pos.front(); i < source_pos.back(); i++){
		N[0][i] = n_photon;
	}
	for(int i = 0; i < (tf/dt) - 1; i++){
		glf_f.clear();
		glf_p.clear();
		glf(i,n_cell,N,F,P,glf_f,glf_p);

		for(int k = 0; k < n_cell; k++){
			if(continuous == true && find(source_pos.begin(), source_pos.end(), k) != source_pos.end()){
				if (k < n_cell -1){
					N[i + 1][k] = n_photon + N[i][k] - dt/dx * (glf_f[k + 1] - glf_f[k]);
					F[i + 1][k] = F[i][k] - dt/dx * (glf_p[k + 1] - glf_p[k]);
				}
				else{
					N[i + 1][k] = n_photon + N[i][k] - dt/dx * (glf_f[0] - glf_f[k]);
					F[i + 1][k] = F[i][k] - dt/dx * (glf_p[0] - glf_p[k]);
				}
			}
			else{
				if (k < n_cell - 1){
                    N[i + 1][k] = N[i][k] - dt/dx * (glf_f[k + 1] - glf_f[k]);
                    F[i + 1][k] = F[i][k] - dt/dx * (glf_p[k + 1] - glf_p[k]);
				}
				else{
                    N[i + 1][k] = N[i][k] - dt/dx * (glf_f[0] - glf_f[k]);
                    F[i + 1][k] = F[i][k] - dt/dx * (glf_p[0] - glf_p[k]);
				}		
			}
			if(N[i + 1][k] < 1e-25){
				N[i + 1][k] = 0.0001;
			}
			if(chemistery == true){
				X[i + 1][k] = third_order_polynom(T[i + 1][k],N[i + 1][k],rho,sigma,X[i][k],dt);
				new_advection_equation(N[i + 1][k],F[i + 1][k],T[i + 1][k],rho,sigma,X[i + 1][k],X[i][k],dt);
			}
			double f = F[i + 1][k] / (c * N[i + 1][k]);
			double chi = (3 + 4 * (f*f)) / (5 + 2 * sqrt(4 - 3 * (f*f)));
			P[i + 1][k] = chi * N[i + 1][k] * (c*c);
		}
	}
}



int main(){
	int n_cell = 61;
	double dt = 1e9;
	double tf = 5e11;
	double x = 0.0012;
	double t = 2e3;

	vector<int> source_pos{30};

	vector<vector<double>> N;
	vector<vector<double>> F;
	vector<vector<double>> P;
	vector<vector<double>> X;
	vector<vector<double>> T;

	vector<double> zeros(n_cell, 0.);
	vector<double> ones(n_cell, 1.);

	for(int i = 0; i < tf/dt; i++){
		N.push_back(zeros);
		F.push_back(zeros);
		P.push_back(zeros);
		X.push_back(ones);
		T.push_back(ones);

		transform(X[i].begin(), X[i].end(), X[i].begin(), bind1st(multiplies<double>(),x));
		transform(T[i].begin(), T[i].end(), T[i].begin(), bind1st(multiplies<double>(),t));
		/*cout << T[i][0] << T[i][15] << T[i][30] << T[i][45] << T[i][60] << endl;*/
	}

	double n_photon = 5e48;
	double sigma = 1.63e-18;
	int rho = 3;
	double dx = 19285e15;
	int c = 3e8;

	bool continuous = true;
	bool chemistery = true;

	cout << "Current condition: " << c*dt / dx << "\n" << endl;

	numerical_scheme(n_photon,N,F,P,X,T,source_pos,dt,tf,n_cell,dx,rho,sigma,continuous,chemistery);
	for(int i = 0; i< tf/dt;i++){
		cout << X[i][0] << "\t" << X[i][10] << "\t"  << X[i][20] << "\t"  << X[i][30] << "\t"  << X[i][40] << "\t"  << X[i][50] << "\t"  << X[i][60] << endl;
	}
	
	/*
	ofstream plot("n.txt");
	for(int i = 0; i < N.size(); i++){
		if(i != 0 && i%6000 == 0){
			for(int j = 0; j < N[i].size(); j++){
				plot<<N[i][j]<<" "<<j<<endl;
			}
			plot<<"\n"<<endl;
		}
	}

	ofstream lex("x.txt");
	for(int i = 0; i < X.size(); i++){
		if(i != 0 && i%6000 == 0){
			for(int j = 0; j < X[i].size(); j++){
				lex<<X[i][j]<<" "<<j<<endl;
			}
			lex<<"\n"<<endl;
		}
	}
	*/
	
	return 0;
}


/*
gnuplot
set title "photon density"
set grid
plot "plot.txt" using 2:1 with lines
*/