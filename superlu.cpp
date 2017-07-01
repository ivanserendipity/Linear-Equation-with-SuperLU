//
//  superlu.cpp
//
//  Created by Ivan zhang on 10/11/16.
//  Copyright © 2016 Ivan zhang. All rights reserved.
//
#include "slu_ddefs.h"
#include <iostream>
#include "Eigen/Core"
#include "Eigen/Sparse"
#include <map>
#include <time.h>
#include <cmath>
#include <complex>
#include <boost/numeric/odeint.hpp>
#include <thread>

using namespace boost::numeric::odeint;

void readCCS(Eigen::SparseMatrix<double>* smt,int r,int c,int nnz, int* ptr,int* ptc,double* value);
Eigen::SparseMatrix<double> eliminate(Eigen::SparseMatrix<double>* smt,int n,int m);
void reform(int type, Eigen::MatrixXd &generator, Eigen::MatrixXd &bus, Eigen::MatrixXd &line,int options);

void solveLinear(std::vector<double> &solve);
Eigen::SparseMatrix<std::complex<double> > admittanceNetwork(Eigen::MatrixXd &bus, Eigen::MatrixXd &line);

std::vector<double> solution(5000);
int i = 0;

void readCCS(Eigen::SparseMatrix<double>* smt,int r,int c,int nnz, int* ptr,int* ptc,double* value){
	int i,j,k;

	if(ptc[0] == 0){
		j = 0;
		for(k = 0;k<nnz;k++){
			i = ptr[k];
			while(ptc[j+1] <= k){
				j = j + 1;			
			}
			smt->coeffRef(i,j) = value[k];		
		}	
	}
	else{
		j = 1;
		for(k = 0;k<nnz;k++){
			i = ptr[k];
			while(ptc[j] <= k+1){
				j = j + 1;			
			}
			smt->coeffRef(i,j) = value[k];
		}	
	}
	return;
}

void eliminate(Eigen::SparseMatrix<double>* smt,double* rhs,int n,int m){

	Eigen::SparseMatrix<double> a = *smt;


	int size = (int)a.rows();
	Eigen::VectorXd rhv(size);
	for(int i = 0;i<size;i++) rhv(i) = rhs[i];

	

	Eigen::SparseMatrix<double> anew;
Eigen::SparseMatrix<double> a11;
Eigen::SparseMatrix<double> a12;
Eigen::SparseMatrix<double> a21;
Eigen::SparseMatrix<double> a22;
	Eigen::VectorXd bnew;
Eigen::VectorXd b1;
Eigen::VectorXd b2;
	clock_t s = clock();

		while(n > m){
			

			a11 = a.block(0,0,n-1,n-1);
			a12 = a.block(0,n-1,n-1,1);
			a21 = a.block(n-1,0,1,n-1);
			a22 = a.block(n-1,n-1,1,1);	
				

			b1 = rhv.head(n-1);
			b2 = rhv.tail(1);

			anew = a11 - (a12 * (1.0 / a22.coeff(0,0)))*a21;
			bnew = b1 - (a12* (1.0 / a22.coeff(0,0)))*b2;
			
	std::cout << a22.coeff(0,0) << std::endl;

			a = anew;
			rhv = bnew;
			n = (int)anew.rows();
		
			a11.data().squeeze();
			a12.data().squeeze();
			a21.data().squeeze();
			a22.data().squeeze();
			anew.data().squeeze();			

	}
	std::cout << "the density time is " << ((double)(clock()-s)/CLOCKS_PER_SEC)*1000 <<" ms" <<std::endl;

	std::cout << anew << std::endl;

	SuperMatrix A,L,U,B;
	SCformat *Lstore;
	NCformat *Ustore;
	mem_usage_t   mem_usage;

	dCreate_CompCol_Matrix(&A, (int)anew.rows(), (int)anew.cols(), (int)anew.nonZeros(), (double*)anew.valuePtr(), (int*)anew.innerIndexPtr(), (int*)anew.outerIndexPtr(), SLU_NC, SLU_D, SLU_GE);

	double *newrhs;
	if ( !(newrhs = doubleMalloc(n)) ) ABORT("Malloc fails for perm_c[].");

	

	for(int i = 0;i<n;i++) newrhs[i] = bnew(i);

	

	dCreate_Dense_Matrix(&B, n, 1, newrhs, n, SLU_DN, SLU_D, SLU_GE);
	//dPrint_Dense_Matrix("B",&B);
	superlu_options_t options;
  	SuperLUStat_t stat;
    
    	set_default_options(&options);
    	options.ColPerm = NATURAL;
    
    	StatInit(&stat);

	int *perm_c; /* column permutation vector */
    	int *perm_r; /* row permutations from partial pivoting */
	int info;
	if ( !(perm_c = intMalloc(n)) ) ABORT("Malloc fails for perm_c[].");
   	if ( !(perm_r = intMalloc(m)) ) ABORT("Malloc fails for perm_r[].");
	
	dgssv(&options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info);

	if ( info == 0 ) {

	/* This is how you could access the solution matrix. */
        double *sol = (double*) ((DNformat*) B.Store)->nzval; 

	 /* Compute the infinity norm of the error. */
	

	Lstore = (SCformat *) L.Store;
	Ustore = (NCformat *) U.Store;
    	std::cout << "No of nonzeros in factor L = " <<Lstore->nnz << std::endl;
    	std::cout << "No of nonzeros in factor U = " <<Ustore->nnz << std::endl;
    	std::cout << "No of nonzeros in L+U = " << Lstore->nnz + Ustore->nnz - n << std::endl;
    	
	
	dQuerySpace(&L, &U, &mem_usage);
	std::cout << "L\\U MB " << mem_usage.for_lu/1e6 << " total MB needed "<< mem_usage.total_needed/1e6 << std::endl;
	
    } else {
	std::cout << "dgssv() error returns INFO= "<<info << std::cout;
	if ( info <= n ) { /* factorization completes */
	    dQuerySpace(&L, &U, &mem_usage);
	    std::cout << "L\\U MB " << mem_usage.for_lu/1e6 << " total MB needed " << mem_usage.total_needed/1e6 << std::endl;
	}
    }

	if ( options.PrintStat ) StatPrint(&stat);
	
	SUPERLU_FREE(newrhs);
    	SUPERLU_FREE(perm_r);
    	SUPERLU_FREE(perm_c);
    
   	Destroy_SuperMatrix_Store(&A);
    	Destroy_SuperMatrix_Store(&B);
    	Destroy_SuperNode_Matrix(&L);
    	Destroy_CompCol_Matrix(&U);
    	StatFree(&stat);
}

void reform(int type, Eigen::MatrixXd &generator, Eigen::MatrixXd &bus, Eigen::MatrixXd &line,int options){

	int i = 0;
	int count = 0;

	while(i < bus.rows()){
		double min = bus(i,4);
		int minIndex = i;
		for(int j = i+1;j < bus.rows();j++){
			if(bus(j,4) < min){
				min = bus(j,4);
				minIndex = j;
			}	
		}
		
		bus.row(i).swap(bus.row(minIndex));
		i++;
	
	}
	
	std::map<double,int> mapping;

		
	for(int i = 0;i<bus.rows();i++){
		mapping[bus(i,0)] = i+1;
	}

	
	for(int i = 0;i<line.rows();i++){
		line(i,0) = mapping[line(i,0)];
		line(i,1) = mapping[line(i,1)];
	}

}

void solveLinear(std::vector<double> &solution){

    SuperMatrix A;
    NCformat *Astore;
    DNformat *Bstore;
    double   *a,*b;
    int      *asub, *xa;
    int      *perm_c; /* column permutation vector */
    int      *perm_r; /* row permutations from partial pivoting */
    SuperMatrix L;      /* factor L */
    SCformat *Lstore;
    SuperMatrix U;      /* factor U */
    NCformat *Ustore;
    SuperMatrix B;
    int      nrhs, ldx, info, m, n, nnz;
    mem_usage_t   mem_usage;
    superlu_options_t options;
    SuperLUStat_t stat;
    FILE      *fp;

    
#if ( DEBUGlevel>=1 )
    CHECK_MALLOC("Enter main()");
#endif

    /* Set the default input options:
	options.Fact = DOFACT;
        options.Equil = YES;
    	options.ColPerm = COLAMD;
	options.DiagPivotThresh = 1.0;
    	options.Trans = NOTRANS;
    	options.IterRefine = NOREFINE;
    	options.SymmetricMode = NO;
    	options.PivotGrowth = NO;
    	options.ConditionNumber = NO;
    	options.PrintStat = YES;
     */
	if(!(fp = fopen("TestData/5000(7).rua","r"))){
		std::cout << "file does not exist" << std::endl;
}
    set_default_options(&options);

    /* Read the matrix in Harwell-Boeing format. */
    dreadhb(fp, &m, &n, &nnz, &a, &b, &asub, &xa);

    dCreate_CompCol_Matrix(&A, m, n, nnz, a, asub, xa, SLU_NC, SLU_D, SLU_GE);
    Astore = (NCformat*)A.Store;
    std::cout << "Dimension " << A.nrow << "x" << A.ncol << " # nonzeros " << Astore->nnz <<std::endl;

    nrhs  = 1;

	dCreate_Dense_Matrix(&B, m, nrhs, b, m, SLU_DN, SLU_D, SLU_GE);
    ldx = n;

    if ( !(perm_c = intMalloc(n)) ) ABORT("Malloc fails for perm_c[].");
    if ( !(perm_r = intMalloc(m)) ) ABORT("Malloc fails for perm_r[].");

    /* Initialize the statistics variables. */
    StatInit(&stat);
    
    dgssv(&options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info);    
 Bstore = (DNformat*)B.Store;
	//dPrint_Dense_Matrix("B",&B);
	
    register int i, lda = Bstore->lda;
    double       *dp;
    
    printf("Stype %d, Dtype %d, Mtype %d\n", B.Stype,B.Dtype,B.Mtype);
    dp = (double *) Bstore->nzval;
    printf("nrow %d, ncol %d, lda %d\n", B.nrow,B.ncol,lda);
    printf("\nnzval: ");
    
    for (i = 0; i < B.nrow; ++i) solution[i] = dp[i];
    
    printf("\n");
    fflush(stdout);

    if ( info == 0 ) {

	/* This is how you could access the solution matrix. */
        double *sol = (double*) ((DNformat*) B.Store)->nzval; 

	 /* Compute the infinity norm of the error. */
	

	Lstore = (SCformat *) L.Store;
	Ustore = (NCformat *) U.Store;
    	std::cout << "No of nonzeros in factor L = " <<Lstore->nnz << std::endl;
    	std::cout << "No of nonzeros in factor U = " <<Ustore->nnz << std::endl;
    	std::cout << "No of nonzeros in L+U = " << Lstore->nnz + Ustore->nnz - n << std::endl;
    	std::cout << "FILL ratio = " << (float)(Lstore->nnz + Ustore->nnz - n)/nnz << std::endl;
	
	dQuerySpace(&L, &U, &mem_usage);
	std::cout << "L\\U MB " << mem_usage.for_lu/1e6 << " total MB needed "<< mem_usage.total_needed/1e6 << std::endl;
	
    } else {
	std::cout << "dgssv() error returns INFO= "<<info << std::cout;
	if ( info <= n ) { /* factorization completes */
	    dQuerySpace(&L, &U, &mem_usage);
	    std::cout << "L\\U MB " << mem_usage.for_lu/1e6 << " total MB needed " << mem_usage.total_needed/1e6 << std::endl;
	}
    }
	

    if ( options.PrintStat ) StatPrint(&stat); 
    StatFree(&stat);

    SUPERLU_FREE (perm_r);
    SUPERLU_FREE (perm_c);
    Destroy_CompCol_Matrix(&A);
    Destroy_SuperMatrix_Store(&B);
    Destroy_SuperNode_Matrix(&L);
    Destroy_CompCol_Matrix(&U);

#if ( DEBUGlevel>=1 )
    CHECK_MALLOC("Exit main()");
#endif

}

Eigen::SparseMatrix<std::complex<double> > admittanceNetwork(Eigen::MatrixXd &bus, Eigen::MatrixXd &line){
	int nSW = 0;
	int nPV = 0;
	int nPQ = 0;
	int swing_bus = 1;
	int gen_bus = 2;
	int load_bus = 3;
	int SB = -1;

	for(int i = 0; i < bus.rows();i++){
		int bus_type = bus(i,9);
		if(bus_type == swing_bus){
			try{
				if(nSW>1) throw;
			}catch(...){
				std::cout<< "Error:nsw cant be greater than 1" << std::endl;
			}
				
			nSW++;		
		}
		else if(bus_type == gen_bus){
			nPV++;		
		}
		else if(bus_type == load_bus){
			nPQ++;
		}
		
	}

	
	Eigen::SparseMatrix<std::complex<double> > yBus(bus.rows(),bus.rows());

	Eigen::MatrixXd from_bus = line.block(0,0,line.rows(),1);
	Eigen::MatrixXd to_bus = line.block(0,1,line.rows(),1);
	Eigen::MatrixXd r = line.block(0,2,line.rows(),1);
	Eigen::MatrixXd rx = line.block(0,3,line.rows(),1); 
	Eigen::MatrixXd chrg = line.block(0,4,line.rows(),1);
	Eigen::MatrixXd tap_ratio = line.block(0,5,line.rows(),1);	

	std::cout << "fb is : " << from_bus << std::endl;
	std::cout << "tb is : " << to_bus << std::endl;
	std::cout << "r is : " << r << std::endl;
	std::cout << "x is : " << rx << std::endl;
	std::cout << "b is : " << chrg << std::endl;
	std::cout << "a is : " << tap_ratio << std::endl;


	Eigen::VectorXcd y = Eigen::VectorXcd::Zero(r.rows());
	Eigen::VectorXcd newb = Eigen::VectorXcd::Zero(r.rows());

	for(int i = 0;i<line.rows();i++){
		y(i) = 1.0 / (r(i,0) + 1i * rx(i,0)); // 1i is jay which equals sqrt(-1)
		newb(i) = 1i * chrg(i,0);
	}

	double nbus = from_bus.maxCoeff();
	double nbranch = from_bus.rows();

	for(int i = 0;i<nbranch;i++){
		yBus.coeffRef(from_bus(i,0) - 1,to_bus(i,0) - 1) = yBus.coeffRef(from_bus(i,0) - 1,to_bus(i,0) - 1) - y(i)/chrg(i,0);
		yBus.coeffRef(to_bus(i,0) - 1,from_bus(i,0) - 1) = yBus.coeffRef(from_bus(i,0) - 1,to_bus(i,0) - 1);
	}


	for(int i = 0;i<nbus;i++){
		for(int j = 0;j<nbranch;j++){
			if(from_bus(j,0) == i+1){
				yBus.coeffRef(i,i) = yBus.coeffRef(i,i) + y(j)/(pow(chrg(j,0),2)) + newb(j);		
			}
			else if(to_bus(j,0) == i+1){
				yBus.coeffRef(i,i) = yBus.coeffRef(i,i) + y(j) +newb(j);
			}		
		}
	}
	
	return yBus;

}

void f( const std::vector<double> &x , std::vector<double> &dx, double t)
{
    for(int i = 0;i<5000;i++){
		dx[i] = solution[i];	
	}
}

void ode(int val){
    
    runge_kutta_dopri5<std::vector<double> > stepper;
    
    
    std::vector<double> x(5000/16);
    const double h = 0.001; // time step
    double t = 0.0;
    stepper.do_step(f,x , t, h );

}

int main(int argc, char *argv[])
{

clock_t start = clock();

	solveLinear(solution);

std::cout <<"Solution Time is " <<((double)(clock()-start)/CLOCKS_PER_SEC)*1000 <<" ms" <<std::endl;

	
clock_t s2 = clock();

    for(int i = 0;i<16;i++){
        std::thread t(ode,i);
        t.join();
    }
    
    std::cout <<"ODE Time is " <<((double)(clock()-s2)/CLOCKS_PER_SEC)*1000 <<" ms" <<std::endl;
	return 0;
}

