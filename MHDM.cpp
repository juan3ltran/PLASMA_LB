#include<iostream>
#include<cmath>
#include "vector.hpp"


const int Lx = 1;  
const int Ly = 45;   
const int Lz = 1; 
const double m0=0.01,m1=18.20,q0=-0.001,q1=0.001;
const double nu=1;  
const double g=0.00005;
const double mu0=0.01;
const double Gamma=1.0;
const double xi = 0.5;
const double taus = 1.0;
const double tau2 = 0.5;
const double magnetic_field_y = 0.04; //0-0.05

class LatticeBoltzmann{
    private:
        int V[3][3][6]; /*V[xyz][p][i]*/  vector3D v[3][6], v0; //v[p][i]
        vector3D e[3][4][2]; vector3D e0;//e[p][i][j] 
        vector3D b[3][4][2]; vector3D b0;//b[p][i][j]        
        double f[Lx*Ly*Lz*3*6*2]; // f[x][y][z][p][i][s]
        double f_new[Lx*Ly*Lz*3*6*2]; 
        double G[Lx*Ly*Lz*3*4*2];//G[Lx][Ly][Lz][p][i][j]
        double G_new[Lx*Ly*Lz*3*4*2];
        double f_0[Lx*Ly*Lz*2], f_0_new[Lx*Ly*Lz*2]; //f[x][y][z][s]
        double G_0[Lx*Ly*Lz], G_0_new[Lx*Ly*Lz]; //G[x][y][z]
        double q[2], m[2];
        double w[6];
        double w0 = 1.0/3.0;
        
    public:
        LatticeBoltzmann(void);
        int nf(int ix, int iy, int iz, int p, int i, int s) {
        return s + 2 * i + 2 * 6 * p + 2 * 6 * 3 * iz + 2 * 6 * 3 * Lz * iy + 2 * 6 * 3 * Lz * Ly * ix;}
        int nG(int x, int y, int z, int p, int i, int j)
        {return x + Lx*y + Lx*Ly*z + Lx*Ly*Lz*p + Lx*Ly*Lz*3*i + Lx*Ly*Lz*3*4*j;};
        int nf0(int ix, int iy, int iz,int s)
        {return s + 2* iz + 2*Lz * iy + 2*Lz * Ly *  ix;};
        int ng0(int x, int y, int z)
        {return x + Lx*y + Lx*Ly*z;};
        double rho_s(int ix, int iy, int iz, int s, bool use_new);
        double feq(double rhoS, vector3D& Vmeds, int p, int i, int s);
        double feq0(double rhoS, vector3D& Vmeds, int s);
        double Geq(vector3D &E_med, vector3D &B, int p, int i, int j);
        double Geq0(void);
        //Fields and magnituds
        vector3D V_s(int ix, int iy, int iz, int s, double& rho_s, bool use_new);
        vector3D V_s_med(int ix, int iy, int iz, int s, double& rho_s, bool use_new);
        vector3D E(int ix, int iy, int iz, bool use_new);
        vector3D E_med(int ix, int it, int iz, bool use_new);
        vector3D B(int ix, int iy, int iz, bool use_new);
        vector3D J_med(int ix, int iy, int iz, bool use_new);
        vector3D F_s(int ix, int iy, int iz, int s, bool use_new);
        //Evolution dynamics
        void Start(vector3D Vmed0);
        void Advection(void);
        void Collision(void);
        void ImposeFields(void);

        //Extract data
        void Print(const char * NameFile);
        void test(void);
        
};

LatticeBoltzmann::LatticeBoltzmann(void){
    v0.load(0,0,0);

    V[0][0][0]=1;  V[0][1][0]=1;  V[0][2][0]=0;
    V[1][0][0]=1;  V[1][1][0]=0;  V[1][2][0]=1;
    V[2][0][0]=0;  V[2][1][0]=1;  V[2][2][0]=1;

    V[0][0][1]=-1; V[0][1][1]=-1; V[0][2][1]=0;
    V[1][0][1]=1;  V[1][1][1]=0;  V[1][2][1]=-1;
    V[2][0][1]=0;  V[2][1][1]=1;  V[2][2][1]=1;

    V[0][0][2]=-1; V[0][1][2]=-1; V[0][2][2]=0;
    V[1][0][2]=-1; V[1][1][2]=0;  V[1][2][2]=-1;
    V[2][0][2]=0;  V[2][1][2]=-1; V[2][2][2]=-1;

    V[0][0][3]=1;  V[0][1][3]=1;  V[0][2][3]=0;
    V[1][0][3]=-1; V[1][1][3]=0;  V[1][2][3]=1;
    V[2][0][3]=0;  V[2][1][3]=-1; V[2][2][3]=-1;

    V[0][0][4]=-1; V[0][1][4]=0;  V[0][2][4]=0;
    V[1][0][4]=0;  V[1][1][4]=-1; V[1][2][4]=0;
    V[2][0][4]=0;  V[2][1][4]=0;  V[2][2][4]=-1;

    V[0][0][5]=1;  V[0][1][5]=0;  V[0][2][5]=0;
    V[1][0][5]=0;  V[1][1][5]=1;  V[1][2][5]=0;
    V[2][0][5]=0;  V[2][1][5]=0;  V[2][2][5]=1;

    for(int p=0;p<3;p++){
        for(int i=0;i<6;i++){
            v[p][i].load(V[0][p][i],V[1][p][i],V[2][p][i]);
        }
    }

    e0.load(0,0,0);
    for(int p=0;p<3;p++){
        for(int i=0; i<4;i++){
            e[p][i][0] = 0.5*v[p][(i+1)%4];
            e[p][i][1] = 0.5*v[p][(i+3)%4];
        }
    }

    b0.load(0,0,0);
    for(int p=0;p<3;p++){
        for(int i=0;i<4;i++){
            for(int j=0;j<2;j++){
	            b[p][i][j]=(v[p][i]^e[p][i][j]);
            } 
        }
    }
    q[0]=q0; q[1]=q1;
    m[0]=m0; m[1]=m1;

    w[0]=w[1]=w[2]=w[3] = 1.0/36.0;
    w[4]=w[5] = 1.0/18.0;

}


double LatticeBoltzmann::rho_s(int ix, int iy, int iz, int s, bool Usenew){
    double sum = 0.0;
    if (Usenew){
        sum += f_0_new[nf0(ix,iy,iz,s)];
        for (int p=0; p<3;p++){
            for (int i=0; i<6;i++){
                sum += f_new[nf(ix,iy,iz,p,i,s)];
            }
        }
    }
    else{
        sum += f_0[nf0(ix,iy,iz,s)];
        for (int p=0; p<3;p++){
            for (int i=0; i<6;i++){
                sum += f[nf(ix,iy,iz,p,i,s)];
            }
        } 
    }
    return sum;
}

vector3D LatticeBoltzmann::V_s(int ix, int iy, int iz, int s, double& rho_s, bool Usenew){
    vector3D sum; sum.load(0,0,0);
    if (Usenew){
        for (int p=0; p<3;p++){
            for (int i=0; i<6;i++){
                sum += f_new[nf(ix,iy,iz,p,i,s)]*v[p][i];
            }
        }
    }
    else{
        for (int p=0; p<3;p++){
            for (int i=0; i<6;i++){
                sum += f[nf(ix,iy,iz,p,i,s)]*v[p][i];
            }
        } 
    }
    return sum /rho_s;
}

vector3D LatticeBoltzmann::E(int ix, int iy, int iz, bool Usenew){
    vector3D sum;
    sum.load(0,0,0);
    if (Usenew){
        for(int p=0;p<3;p++){
            for(int i=0; i<4;i++){
                for(int j=0;j<2;j++){
                    sum+=G_new[nG(ix,iy,iz,p,i,j)]*e[p][i][j];
                }
            }
        }
    }
    else{
        for(int p=0;p<3;p++){
            for(int i=0; i<4;i++){
                for(int j=0;j<2;j++){
                    sum+=G[nG(ix,iy,iz,p,i,j)]*e[p][i][j];
                }
            }
        }
    }
    return sum;
}

vector3D LatticeBoltzmann::B(int ix, int iy, int iz, bool Usenew){
    vector3D sum;
    sum.load(0,0,0);
    if (Usenew){
        for(int p=0;p<3;p++){
            for(int i=0; i<4;i++){
                for(int j=0;j<2;j++){
                    sum+=G_new[nG(ix,iy,iz,p,i,j)]*b[p][i][j];
                }
            }
        }
    }
    else{
        for(int p=0;p<3;p++){
            for(int i=0; i<4;i++){
                for(int j=0;j<2;j++){
                    sum+=G[nG(ix,iy,iz,p,i,j)]*b[p][i][j];
                }
            }
        }
    }
    return sum;
}

vector3D LatticeBoltzmann::J_med(int ix, int iy, int iz, bool use_new){
    vector3D sum;
    sum.load(0,0,0);
    for (int s=0;s<2;s++){
        double rhoS = LatticeBoltzmann::rho_s(ix,iy,iz,s,use_new);
        sum += (q[s]/m[s])*rhoS*LatticeBoltzmann::V_s_med(ix,iy,iz,s,rhoS,use_new);
    }
    return sum;
}

vector3D LatticeBoltzmann::F_s(int ix, int iy, int iz, int s, bool use_new){
    vector3D E0, B0, VS, VSm1, F0;
    E0 = LatticeBoltzmann::E(ix,iy,iz,false); B0 = LatticeBoltzmann::B(ix,iy,iz,false);
    double rhoS = LatticeBoltzmann::rho_s(ix,iy,iz,s,false); double rhoSm1 = LatticeBoltzmann::rho_s(ix,iy,iz,(s+1)%2,false);
    VS = LatticeBoltzmann::V_s(ix,iy,iz,s,rhoS,false); VSm1 = LatticeBoltzmann::V_s(ix,iy,iz,(s+1)%2,rhoSm1,false);
    //F0.load(0,0,0);
    F0.load(rhoS*g,0,0);
    return  F0 + (q[s]/m[s])*rhoS*(E0+(VS^B0))- nu*rhoS*(VS-VSm1); 
}

vector3D LatticeBoltzmann::V_s_med(int ix, int iy, int iz, int s, double& rho_s, bool use_new){
    vector3D VS, FS;
    VS = LatticeBoltzmann::V_s(ix,iy,iz,s,rho_s,use_new);
    FS = LatticeBoltzmann::F_s(ix,iy,iz,s,use_new);
    return VS + 0.5*FS/rho_s;
}

vector3D LatticeBoltzmann::E_med(int ix, int iy, int iz, bool use_new){
    vector3D E0, J0;
    E0 = LatticeBoltzmann::E(ix,iy,iz,use_new);
    J0 = LatticeBoltzmann::J_med(ix,iy,iz,use_new);
    return E0 - 0.25*mu0*J0;
}

double LatticeBoltzmann::feq(double rhoS, vector3D& Vmeds, int p, int i, int s){
   return w[i]*rhoS*(3*xi*pow(rhoS,Gamma-1) + 3*v[p][i]*Vmeds + 9*pow(v[p][i]*Vmeds,2) - 1.5*Vmeds.norm2());
}

double LatticeBoltzmann::feq0(double rhoS, vector3D &Vmeds, int s){
    return 3*rhoS*(1-0.5*(4*xi*pow(rhoS,Gamma-1) + Vmeds.norm2()));
}

double LatticeBoltzmann::Geq(vector3D &E_med, vector3D &B, int p, int i, int j){
    return 0.25*E_med*e[p][i][j] + 0.125*B*b[p][i][j];
}

double LatticeBoltzmann::Geq0(void){
    return 0;
}

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

void LatticeBoltzmann::Print(const char * NameFile)
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

     

void Animation(void)
{
    std::cout<<"set terminal gif animate "<<std::endl;
    std::cout<<"set output 'campo_vectorial.gif'"<<std::endl;
    std::cout<<"set xrange [0:256]    # Ajusta según tus datos"<<std::endl;
    std::cout<<"set yrange [0:64]    # Ajusta según tus datos"<<std::endl;
    std::cout<<"set size ratio -1    # Para mantener la proporción"<<std::endl;
}

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
    
    

  plasma.Print("null.dat");
    return 0;
}