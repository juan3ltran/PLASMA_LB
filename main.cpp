#include "Lattice.hpp"     


int main(void){
    LatticeBoltzmann plasma;

    vector3D velocity0;
    velocity0.load(0,0,0);
    int t, taux=0, tmax=3000;



    plasma.Start(velocity0); 
    
   
    for ( t = 0; t < tmax; t++)
    {   
        
        plasma.Collision();
        
        
        plasma.ImposeFields();
        
        
        plasma.Advection();
        

        
    }
    
    

  plasma.Print();
    return 0;
}