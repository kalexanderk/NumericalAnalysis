#include <math.h>
#include <fstream>
#include <iostream>
#include <string>
using namespace std;

long double task (double t) {
	return log(t*t/(sqrt(1-8*t*t*t*t)));
}

double *GausssFunc(int dim, double **A, double *p) { //Ax=p
			int i, j, k, l;
			double koef;
      		double max_el;
      		double **matr;
      		double *q=new double [dim];
      		matr=new double *[dim];
			for(i=0; i<dim; i++) {
				matr[i]=new double[dim+1];	
			}
			for(i=0; i<dim; i++)
          		for(j=0; j<dim+1; j++) {
            		if (j!=dim) matr[i][j]=A[i][j];
            		else matr[i][j]=p[i];
            }
			for(k=0; k<dim; k++) {
          		max_el=matr[k][k];
	      		for(i=k; i<dim; i++) if (fabs(matr[i][k])>max_el) max_el=fabs(matr[i][k]);
	      		if (max_el==0) return 0;
 				for(l=k+1; l<dim; l++) {
		      		koef=matr[l][k]/matr[k][k];
	          		for(j=0; j<dim+1; j++) 
		        		matr[l][j]+=-matr[k][j]*koef;
          		} 
      		}
      		q[dim-1]=matr[dim-1][dim]/matr[dim-1][dim-1];
			for(i=dim-2; i>=0; i--) {
	      		q[i]=matr[i][dim]/matr[i][i];
	      		for(j=i+1; j<dim; j++) q[i]+=-matr[i][j]/matr[i][i]*q[j];
	  		}
	  		return q;
		}

class Interpolator {
	public:
		double begin;
		double end;
		int n; //n+1 nods
		double *x; //nod's x
		double *y; //nod's y
		double *f; //coefs of Newton's polymon
		double **cub; //for cubic spline
		Interpolator(int N) {
			int i;
			cout<<"Input the beginning of interpolation: ";
			do {
				i=scanf("%lf", &begin);
				fflush(stdin);
				if (i!=1) cout<<"Error, enter value again: ";
			} while(i!=1);
			cout<<"Input the ending if interpolation: ";
			do {
				i=scanf("%lf", &end);
				fflush(stdin);
				if (i!=1) cout<<"Error, enter value again: ";
			} while(i!=1);
			n=N;
			x=new double[n+1];
			y=new double[n+1];
			f=new double[n+1];
			cub=new double*[n];
			for (i=0; i<n; i++) cub[i]=new double[3];
		}
		void Transfer() {
			double a;
			a=begin; begin=end; end=a;
		}
		void GetNet(long double q(double)) {
			int i;
			double step=(end-begin)/n;
			double tx=begin;
			for (i=0; i<n+1; i++) {
				x[i]=tx;
				y[i]=q(tx);
				tx+=step;
			}
			x[n]=end; y[n]=q(end);
		}
		void NewtonInterpolator() {//creates Newt Intpltr ans solves the system
			int i,j=0,k;
			double **fcc;
			fcc=new double *[n+1];
			for(i=0; i<n+1; i++) fcc[i]=new double [n+1];
			for (i=0; i<n+1; i++) {
				fcc[i][0]=1;
				for (j=0; j<n+1; j++) {
					if (i<j) fcc[i][j]=0;
					if (i>0) fcc[i][1]=x[i]-x[0];
					if (i>1) for (k=2; k<i+1; k++) {
						fcc[i][k]=(x[i]-x[k-1])*fcc[i][k-1];
					}
				}
			}
			f=GausssFunc(n+1,fcc,y);
		}
		double Newton(double t) {
			int i, j;
			double S=f[0],P;
			for (i=1; i<n+1; i++) {
				for (j=0, P=f[i]; i>j; j++) P*=(t-x[j]);
				S+=P;
			}
			return S;
		}
		double Lagrange(double t) {
			int i,j;
			double S, P1, P2;
			for (i=0; i<n+1; i++) {
				for (j=0, P1=y[i], P2=1; j<n+1; j++) 
					if (i!=j) {
						P1*=(t-x[j]);
						P2*=(x[i]-x[j]);
					}
				S+=P1/P2;
			}
			return S;
		}
		void CubicSplineInterpolator() {//finds coeff of eb cub functions
			int i,j;
			double *h=new double[n];
			double *r=new double[n+1];
			double *p=new double[n+1];
			double **s;
			s=new double *[n+1];
			for(i=0; i<n+1; i++) s[i]=new double[n+1];
			for (i=0; i<n; i++) h[i]=x[i+1]-x[i];
			for (i=0; i<n+1; i++) {//formula(to solve equ)
				if ((i==0) or (i==n)) r[i]=0;
				else r[i]=((y[i+1]-y[i])/h[i]-(y[i]-y[i-1])/h[i-1])*6;
			}
			for (i=0; i<n+1; i++) {//3diagonal matrix
				for (j=0; j<n+1; j++) {
					if (i==j) {
						if ((i!=0)&&(i!=n)) s[i][j]=2*(h[i-1]+h[i]);
						else s[i][j]=1;
					}
				else if (j==i-1) {
					if (i!=n) s[i][j]=h[i-1];
					else s[i][j]=0;
				}
				else if (j==i+1) {
					if (i!=0) s[i][j]=h[i];
					else s[i][j]=0;
				}
				else s[i][j]=0;
				}
			}///////////////////////
			r=GausssFunc(n+1,s,r);
			for(i=0; i<n; i++) {
				cub[i][1]=r[i+1]*0.5;
				cub[i][2]=(r[i+1]-r[i])/(h[i]*6);
				cub[i][0]=(h[i]*cub[i][1]-h[i]*h[i]*cub[i][2]+(y[i+1]-y[i])/h[i]);
			}
		}
		long double CubicSpline(double t) {//finds f vls in points t
			int i;
			if (x[0]<x[n]) {
				if (t<x[0]) i=0;
				else if (t>=x[n]) i=n-1;
				else i=(t-x[0])*n/(x[n]-x[0]);	
			}
			else {
				if (t>x[0]) i=0;
				else if (t<=x[n]) i=n-1;
				else i=(t-x[0])*n/(x[n]-x[0]);
	
			}
			double dx=t-x[i+1];				
			return cub[i][2]*dx*dx*dx+cub[i][1]*dx*dx+cub[i][0]*dx+y[i+1];
		}
		void Out(long double q(double), string s) {
			double tx=x[0], k; 
			double step=(end-begin)/(n);
			ofstream fout(s,ios::out);
			fout << "|  x  |  F(x)   |(Newton-F)(x)|(Lagrange-F)(x)|(CubSpline-F)(x)|" << endl;
			do {
				fout <<"|";
				fout.precision(3);
				fout.width(5);
				fout << tx << "|";
				fout.precision(6);
				fout.width(9);
				fout << q(tx) << "|";
				fout.precision(6);
				fout.width(13);
				fout << fixed << fabs(Newton(tx)-q(tx)) << "|";
				fout.precision(6);
				fout.width(15);
				fout << fixed << fabs(Lagrange(tx)-q(tx)) << "|";
				fout.precision(6);
				fout.width(16);
				fout << fixed << fabs(CubicSpline(tx)-q(tx)) << "|";
				fout << endl;
				tx+=step;
			} while (((tx<x[n]) && (step>0))||((step<0)&&(tx>x[n])));
			tx=x[n];
			fout <<"|";
				fout.precision(3);
				fout.width(5);
				fout << tx << "|";
				fout.precision(6);
				fout.width(9);
				fout << q(tx) << "|";
				fout.precision(6);
				fout.width(13);
				fout << fixed << fabs(Newton(tx)-q(tx)) << "|";
				fout.precision(6);
				fout.width(15);
				fout << fixed << fabs(Lagrange(tx)-q(tx)) << "|";
				fout.precision(6);
				fout.width(16);
				fout << fixed << fabs(CubicSpline(tx)-q(tx)) << "|";
				fout << endl;
			fout << "#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@" << endl;
			tx = x[0]; step=step/5;  int i=0;
			do {
				if (i==5*n) tx=x[n];
				fout <<"|";
				fout.precision(3);
				fout.width(5);
				fout << tx << "|";
				fout.precision(6);
				fout.width(9);
				fout << q(tx) << "|";
				fout.precision(6);
				fout.width(13);
				fout << fixed << fabs(Newton(tx)-q(tx)) << "|";
				fout.precision(6);
				fout.width(15);
				fout << fixed << fabs(Lagrange(tx)-q(tx)) << "|";
				fout.precision(6);
				fout.width(16);
				fout << fixed << fabs(CubicSpline(tx)-q(tx)) << "|";
				fout << endl;
				tx+=step;
				i++;
			} while (((tx<x[n]) && (step>0))||((step<0) && (tx>x[n])));
			tx=x[n];
			fout <<"|";
			fout.precision(3);
			fout.width(5);
			fout << tx << "|";
			fout.precision(6);
			fout.width(9);
			fout << q(tx) << "|";
			fout.precision(6);
			fout.width(13);
			fout << fixed << fabs(Newton(tx)-q(tx)) << "|";
			fout.precision(6);
			fout.width(15);
			fout << fixed << fabs(Lagrange(tx)-q(tx)) << "|";
			fout.precision(6);
			fout.width(16);
			fout << fixed << fabs(CubicSpline(tx)-q(tx)) << "|";
			fout << endl;
		}
};

int main() {
	Interpolator B(5);
	B.GetNet(task);
	B.NewtonInterpolator();
	B.CubicSplineInterpolator();
	B.Out(task,"output_up.txt");
	cout<<endl;
	B.Transfer();
	B.GetNet(task);
	B.NewtonInterpolator();
	B.CubicSplineInterpolator();
	B.Out(task,"output_down.txt");
}
