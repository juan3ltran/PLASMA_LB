#include<iostream>
#include<cmath>
#include "vector.hpp"

const int a=4;
const int b=5;
const int c=3;
const int d=1;
const int e=1;

int n(int i,int j,int k,int z)
{
    double index;
    index = i + a*j + a*b*k + a*b*c*z;
    return index;
}


int main()
{
    
    
    int matriz[a][b][c][d];
    double vector[a*b*c*d]; 
    double time=0;

    for (int k = 0; k < c; k++)
    {
        for (int j = 0; j < b; j++)
        {
            for (int i = 0; i < a; i++)
            {
                for (int z = 0; z < d; z++)
                {
                    matriz[i][j][k][z] = time;
                    std::cout<<"Elemento: "<<i<<j<<k<<z<<"Index: "<<n(i,j,k,z)<<std::endl;
                    std::cout<<"Matriz: "<<matriz[i][j][k][z]<<"Vector: "<<vector[n(i,j,k,z)]<<std::endl;
                    time++;
                }
                
            }
            
        }
        
    }
    
    

    

    return 0;
}