#include "Lattice.hpp"

void LatticeBoltzmann::Advection(void){
    int ixnew, iynew, iznew;
    for(int ix=0; ix<Lx;ix++){
        for(int iy=0;iy<Ly;iy++){
            for(int iz=0;iz<Lz;iz++){
                for(int p=0;p<3;p++){
                    for(int i=0; i<6;i++){
                        ixnew=(ix+V[0][p][i]+Lx)%Lx; iynew=(iy+V[1][p][i]+Ly)%Ly; iznew=(iz+V[2][p][i]+Lz)%Lz;
                        for(int s_j=0;s_j<2;s_j++){
                           if (i<4){
                                G[nG(ixnew,iynew,iznew,p,i,s_j)] = G_new[nG(ix,iy,iz,p,i,s_j)];
                                G_0[ng0(ix,iy,iz)] = G_0_new[ng0(ix,iy,iz)];
                            }
                            f[nf(ixnew,iynew,iznew,p,i,s_j)] = f_new[nf(ix,iy,iz,p,i,s_j)];
                            f_0[nf0(ix,iy,iz,s_j)] = f_0_new[nf0(ix,iy,iz,s_j)];
                            
                        }                        
                    }
                }
            }
        }
    }
}

void LatticeBoltzmann::Collision(void){
    vector3D E0, B0, Vmed0, Fs, Jmed0;
    double FEQS, FEQ0S, GEQ, GEQ0, rhoS;
    for(int ix=0; ix<Lx ; ix++){
        for(int iy=0; iy<Ly ; iy++){
            for(int iz=0; iz<Lz ; iz++){
                E0 = LatticeBoltzmann::E_med(ix,iy,iz,false); 
                B0 = LatticeBoltzmann::B(ix,iy,iz,false);  
               
                GEQ0 = LatticeBoltzmann::Geq0();
                G_0_new[ng0(ix,iy,iz)] = G_0[ng0(ix,iy,iz)] - (1/tau2)*(G_0[ng0(ix,iy,iz)]-GEQ0);              
                 for(int p=0;p<3;p++){
                    for(int i=0; i<6;i++){
                        for(int s_j=0;s_j<2;s_j++){
                            if (i<4){
                               Jmed0 = LatticeBoltzmann::J_med(ix,iy,iz,false);
                               GEQ = LatticeBoltzmann::Geq(E0,B0,p,i,s_j);
                               G_new[nG(ix,iy,iz,p,i,s_j)] = G[nG(ix,iy,iz,p,i,s_j)]  - (1/tau2)*(G[nG(ix,iy,iz,p,i,s_j)]-GEQ) - ((2*tau2-1)/(16*tau2))*mu0*e[p][i][s_j]*Jmed0;
                               
                            }
                            rhoS = LatticeBoltzmann::rho_s(ix,iy,iz,s_j,false);
                            Vmed0 = LatticeBoltzmann::V_s_med(ix,iy,iz,s_j,rhoS,false);
                            FEQS =  LatticeBoltzmann::feq(rhoS,Vmed0,p,i,s_j);
                            FEQ0S = LatticeBoltzmann::feq0(rhoS, Vmed0,s_j);
                            Fs = LatticeBoltzmann::F_s(ix,iy,iz,s_j,false);
                            f_0_new[nf0(ix,iy,iz,s_j)] = f_0[nf0(ix,iy,iz,s_j)] - (1.0/taus)*(f_0[nf0(ix,iy,iz,s_j)] - FEQ0S);
                            f_new[nf(ix,iy,iz,p,i,s_j)] = f[nf(ix,iy,iz,p,i,s_j)]  - (1.0/taus)*(f[nf(ix,iy,iz,p,i,s_j)]-FEQS)+((2*taus-1)/(20*taus))*v[p][i]*Fs; 
                        }                        
                    }
                }
                            
            }
        }
    }        
}

void LatticeBoltzmann::Start( vector3D Vmed0)
{

    double rho0 = 1.0;
    double rho1 = 1820;
    double B0norm = magnetic_field_y;//0.0-0.5
    vector3D E0, B0;
    E0.load(0,0,0); B0.load(0,B0norm,0);
    for (int ix = 0; ix < Lx; ix++) //for each cell
        for (int iy = 0; iy < Ly; iy++)
            for (int iz = 0; iz < Lz; iz++){
                for (int p = 0; p < 3; p++){ //in each plane
                    for (int i = 0; i < 6; i++){// at each direcction                                                 
                        
                        f[nf(ix,iy,iz,p,i,0)] = feq(rho0, Vmed0, p, i, 0);
                        f_0[nf0(ix,iy,iz,0)] = feq0(rho0, Vmed0, 0);
                        f[nf(ix,iy,iz,p,i,1)] = feq(rho1, Vmed0, p, i, 1);
                        f_0[nf0(ix,iy,iz,1)] = feq0(rho1, Vmed0, 1);
                        if (i<4){
                            G[nG(ix,iy,iz,p,i,0)] = Geq(E0,B0,p,i,0);
                            G[nG(ix,iy,iz,p,i,1)] = Geq(E0,B0,p,i,1);
                            G_0[ng0(ix,iy,iz)] = Geq0();                            
                        }
                        
                    }
                }
        }    
}
void LatticeBoltzmann::ImposeFields(void)
{
    int ix,iy,iz,p,i,s;
    vector3D V0; V0.load(0,0,0);
    iy=0; 
    for(ix=0;ix<Lx;ix++){
        for(iz=0;iz<Lz;iz++){
            for(p=0;p<3;p++){
                for(i=0;i<6;i++){
                    for(s=0;s<2;s++){
                        double rho0=rho_s(ix,iy,iz,s,false);                        
                        f_0_new[nf0(ix,iy,iz,s)] = feq0(rho0,V0,s);
                        f_new[nf(ix,iy,iz,p,i,s)]=feq(rho0,V0,p,i,s);
                    }
                }
            }
        }
    }
  
    iy=Ly-1; 
      for(ix=0;ix<Lx;ix++){
        for(iz=0;iz<Lz;iz++){
            for(p=0;p<3;p++){
                for(i=0;i<6;i++){
                    for(s=0;s<2;s++){
                        double rho0=rho_s(ix,iy,iz,s,false);                        
                        f_0_new[nf0(ix,iy,iz,s)] = feq0(rho0,V0,s);
                        f_new[nf(ix,iy,iz,p,i,s)]=feq(rho0,V0,p,i,s);
                    }
                }
            }
          
        }
      
      }

}

void LatticeBoltzmann::Print(void)
{
    double rho0 ,rho1, rho;
    int ix,iy,iz;
    vector3D velocity0, velocity1, velocity;
    for (ix = 0; ix < Lx; ix++)
    {
        for (iy = 0; iy < Ly; iy++){

        for (iz=0;iz<Lz;iz++){
            rho0 = rho_s(ix, iy, iz, 0, true);
            rho1 = rho_s(ix,iy,iz,1,true);
            rho = rho0 + rho1;

            velocity0 = V_s_med(ix, iy, iz, 0, rho0, true);
            velocity1 = V_s_med(ix, iy, iz, 1, rho0, true);
            velocity = (rho*velocity0 + rho1*velocity1)/rho;

            vector3D E0=E_med(ix,iy,iz,true);
            vector3D B0 = B(ix,iy,iz,true);
            
            std::cout<<iy<<" "<<velocity.x()<<" "<<B0.x()<<" "<<B0.y()<<std::endl;
            
            }}
        std::cout<<std::endl;
    }
    std::cout << "e" << std::endl; 
}