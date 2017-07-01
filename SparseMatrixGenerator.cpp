//
//  main.cpp
//  SparseMatrixGenerator
//
//  Created by Ivan zhang on 10/4/16.
//  Copyright © 2016 Ivan zhang. All rights reserved.
//

#include <iostream>
#include "Eigen/Core"
#include "unsupported/Eigen/SparseExtra"


Eigen::MatrixXd x;
void generateSparseMatrix(int rows,int nnz);

void generateSparseMatrix(int rows,int nnz){
    Eigen::SparseMatrix<double,Eigen::ColMajor> mat(rows,rows);
    
    std::cout << mat.innerNonZeroPtr() <<std::endl;
    
    
    for(int i=0; i<rows;i++){
        for(int j = 0;j<nnz;j++){
            int c = rand() % rows;
            mat.coeffRef(i, c) = (double)(rand()% 50 + 1);
        }
    }
    
    Eigen::saveMarket(mat,"/Users/Serendipity/Documents/C++/SuperLU/SuperLU/1000Tests.mtx");
    
}

int main(int argc, const char * argv[]) {
    // insert code here...
    
    
    return 0;
}
