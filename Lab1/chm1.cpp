#include <iostream>
#include <string>
#include <cmath>
#include <fstream>
using namespace std;

typedef float (*func) (float);

float min(float x, float y) {
	if (x>=y) return y;
	else return x;
}

float max(float x, float y) {
	if (x>=y) return x;
	else return y;
}

float bisections(float a, float b, float eps, func f) {
	if (f(a)==0) return a;
	if (f(b)==0) return b;
	float t1=a,t2=b,t3, p1, p2;
	int i=0;
	while(abs(t2-t1)>eps || abs(f(t1))>eps || abs(f(t2))>eps) {
		i++;
		t3=(t2+t1)/2;
		if (f(t3)==0) return t3;
		if (f(t3)*f(t2)<0) t1=t3;
		else t2=t3;
		p1=f(t1); p2=f(t2);
		cout<<i<<" iteration on ["<<t1<<","<<t2<<"];"<<endl;
		cout<<"f("<<t1<<")="<<p1<<", f("<<t2<<")="<<p2<<endl;
	}
	return t3;
}

float hord(float a, float b, float eps, func f) {
	float t1=a, t2=b, t3, p1, p2, d;
	int i=0;
	while(abs(t3-d)>eps || abs(f(t3))>eps) {
		i++;
		d=t3;
		t3=t1-f(t1)*(t2-t1)/(f(t2)-f(t1));
		if (f(t3)==0) return t3;
		if (f(t2)*f(t3)<=0) t1=t3;
		else t2=t3;
		p1=f(min(t2,t1)); p2=f(max(t1,t2));
		cout<<i<<" iteration on ["<<min(t1,t2)<<","<<max(t1,t2)<<"];"<<endl;
		cout<<"f("<<min(t1,t2)<<")="<<p1<<", f("<<max(t1,t2)<<")="<<p2<<endl;
	}
	return t3;
}

float N(float a, float eps, func f, func df) {
	float t1=a, t2=a-f(a)/df(a), rab;
	int i=0;
	while(abs(t2-t1)>eps || abs(f(t2))>eps) {
		i++;
		cout<<i<<" iteration: "<<t2<<","<<" f("<<t2<<")="<<f(t2)<<endl;
	    t1=t2;
		t2=t1-f(t1)/df(t1);
	}
	return t2;
}

float f(float x) {
	return -x*x*x*x + 3*x*x*x - 2*x +3;
}

float df(float x){
	return -4*x*x*x + 9*x*x -2;
}

int main() {
	float a,b,eps=0.00001,res;
	int n,l;
	cout<<"1 - Bisections: "<<endl;
	cout<<"2 - Hords: "<<endl;
	cout<<"3 - Newton: "<<endl;
	do {
    cout<<"Input a=";
			do{
			   l=scanf("%f",&a);
		   	   fflush(stdin);
			   if (l!=1) cout<<"Error, input again: ";
			} while (l!=1);
		cout<<"Input b=";
			do{
			   l=scanf("%f",&b);
		   	   fflush(stdin);
			   if (l!=1) cout<<"Error, input again: ";
			} while (l!=1);
	if (f(a)*f(b)>0 && f(a)*f((a+b)/2)>0 && f((a+b)/2)*f(b)>0) cout<<"There are maybe no roots, input again more correct:"<<endl;		
	} while(f(a)*f(b)>0 && f(a)*f((a+b)/2)>0 && f((a+b)/2)*f(b)>0);
		cout<<"1:"<<endl;	
		res=bisections(a,b,eps,&f);
		cout<<"Result="<<res<<endl;
		cout<<"2:"<<endl;
		res=hord(a,b,eps,&f);
		cout<<"Result="<<res<<endl;
		cout<<"3:"<<endl;
		cout<<"To begin iteration took x0=0.5*(a+b)"<<endl;
		a=(a+b)/2;
		res=N(a,eps,&f,&df);
		cout<<"Result="<<res;
}

