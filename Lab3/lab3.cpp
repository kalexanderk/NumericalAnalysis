#include <iostream>
#include <fstream>
#include <cmath>

#define n 5
#define eps 0.00001

using namespace std;


double Residual(double** A, double* x, int k){
	double *r=new double[k], tmp, norm;
	for(int i=0; i<k; i++){
		tmp=0;
		for(int j=0; j<k; j++)
			tmp+=A[i][j]*x[j];
		r[i]=tmp-A[i][k];
	}
	cout << "Vector of residual: "<< endl;
	for (int i=0; i<k; i++){
		cout.width(15);
		cout << r[i];
	}

	norm = 0;
	for(int i = 0; i < k; i++) norm = norm + r[i]*r[i];
	norm = sqrt(norm);
	cout<<"\nNorma: "<<norm<<"\n\n";
	
	return norm;
}

int simple_iteration(double **mat, double *res, int m){
    
    int i, j;
    
    //check for zeros on the diagonal
	for(i = 0; i < m; i++) {  
	    if(mat[i][i] == 0) return 1;
    }
	
	//install values for free values (the last one column in our matrix)
	double *free_vec = new double[m];
	for(i = 0; i < m; i++)
	    free_vec[i] = mat[i][m]/mat[i][i];
	
	//directly our cycle
	double *temp = new double[m];
	double sum;
	double norm;
	int count = 0;
	
	do {	
		//resend values of results to temporary array to use 'em later/again&again
	    for(i = 0; i < m; i++)
	        temp[i] = res[i];
	    
	    //calculation of the values of x^i
		for(i = 0; i < m; i++){
            sum = 0;
			for(j = 0; j < m; j++) if(i != j) sum -= res[j]*mat[i][j]/mat[i][i];
			res[i] = sum + free_vec[i];  
        }
	    
	    //restriction on the number of iterations
	    count++;
	    if(count > 10000) return 0;
	
		//calculation of the norm of residual
		norm = 0;
		for(i = 0; i < m; i++){
        	norm = norm + (res[i] - temp[i])*(res[i] - temp[i]);
		}
		norm = sqrt(norm);
		
		cout<<"Iteration "<<count<<".\n";
		cout<<"Roots: { ";
		for (i=0;i<m-1;i++) cout<<res[i]<<", ";
		cout<<res[m-1]<<" }\n";
		norm=Residual(mat, res, m);
		
	} while(norm > eps);
	
	return 0;
}

int main() {
	
	int i, j, k;
	
	//results initialisation
	double *result = new double[n];
	
	//matrix initialization
	double **A = new double*[n];
	for(i = 0; i < n; i++)
	    A[i] = new double[n + 1];
	
	//read matrix from "matrix.dat"
	ifstream input_file("matrix.dat");   
	char temp[10];	
	for(i = 0; i < n; i++)
		for(j = 0; j < n + 1; j++){
				input_file>>A[i][j];
				if (input_file.fail()) return 0;					
		}
	input_file.close();
	
	//matrix output
	cout << "Matrix A|b: \n";
	for(i = 0; i < n; i++){
	    for(j = 0; j < n + 1; j++)
	        cout << A[i][j] << " ";
	    cout << "\n";
    }
    cout <<"\n";
	 
	//initial values of the results    
	for(i = 0; i < n; i++)
	    result[i] = A[i][n]/A[i][i];    
	
	//execute iteration method
	if(simple_iteration(A, result, n) != 0) cout << "\nError!";
	else {
		  cout << "\nRoots: ";
		  for(i = 0; i < n; i++)
	          cout << result[i] << " ";
    }
	
	return 0;
}

