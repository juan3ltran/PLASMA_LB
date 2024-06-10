#include<iostream>
#include<cmath>
#include "vector.hpp"

const int Lx = 1;   //
const int Ly = 1;   //
const int Lz = 1000; //
const double m0=0,m1=0,q0=0,q1=0;
const double nu_s=0;  
const double g=0;
const double mu0=0;
const double Gamma=0;

class LatticeBoltzmann{
    private:
        int V[3][3][6]; /*V[xyz][p][i]*/  vector3D v[3][6], v0; //v[p][i]
        vector3D e[3][4][2]; vector3D e0;//e[p][i][j] 
        vector3D b[3][4][2]; vector3D b0;//b[p][i][j]        
        double f[Lx][Ly][Lz][3][6][2]; // f[x][y][z][p][i][s]
        double f_new[Lx][Ly][Lz][3][6][2]; 
        double G[Lx][Ly][Lz][3][4][2];//G[Lx][Ly][Lz][p][i][j]
        double G_new[Lx][Ly][Lz][3][4][2];
        double f_0[Lx][Ly][Lz][2], f_0_new[Lx][Ly][Lz][2];
        double G_0[Lx][Ly][Lz], G_0_new[Lx][Ly][Lz];
        double q[2], m[2];
        double w[6];
        double w0 = 1.0/3.0;
        double xi[2];
        double taus[2];
        double tau2;
    public:
        LatticeBoltzmann(void);
        double rho_s(int ix, int iy, int iz, int s, bool use_new);
        double feq(double rhoS, vector3D& Vmeds, int p, int i, int s);
        double feq0(double rhoS, vector3D& Vmeds, int s);
        double Geq(vector3D &E_med, vector3D &B, int p, int i, int j);
        double Geq0(void);
        vector3D V_s(int ix, int iy, int iz, int s, double& rho_s, bool use_new);
        vector3D V_s_med(int ix, int iy, int iz, int s, double& rho_s, bool use_new);
        vector3D E(int ix, int iy, int iz, bool use_new);
        vector3D E_med(int ix, int it, int iz, bool use_new);
        vector3D B(int ix, int iy, int iz, bool use_new);
        vector3D J_med(int ix, int iy, int iz, bool use_new);
        vector3D F_s(int ix, int iy, int iz, int s, bool use_new);
        void Advection(void);
        void Collision(void);
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

    taus[0] = taus[1] = 0;
}


double LatticeBoltzmann::rho_s(int ix, int iy, int iz, int s, bool Usenew){
    double sum = 0;
    if (Usenew){
        sum += f_0_new[ix][iy][iz][s];
        for (int p=0; p<3;p++){
            for (int i=0; i<6;i++){
                sum += f_new[ix][iy][iz][p][i][s];
            }
        }
    }
    else{
        sum += f_0[ix][iy][iz][s];
        for (int p=0; p<3;p++){
            for (int i=0; i<6;i++){
                sum += f[ix][iy][iz][p][i][s];
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
                sum += f_new[ix][iy][iz][p][i][s]*v[p][i];
            }
        }
    }
    else{
        for (int p=0; p<3;p++){
            for (int i=0; i<6;i++){
                sum += f[ix][iy][iz][p][i][s]*v[p][i];
            }
        } 
    }
    return sum/rho_s;
}

vector3D LatticeBoltzmann::E(int ix, int iy, int iz, bool Usenew){
    vector3D sum;
    sum.load(0,0,0);
    if (Usenew){
        for(int p=0;p<3;p++){
            for(int i=0; i<4;i++){
                for(int j=0;j<2;j++){
                    sum+=G_new[ix][iy][iz][p][i][j]*e[p][i][j];
                }
            }
        }
    }
    else{
        for(int p=0;p<3;p++){
            for(int i=0; i<4;i++){
                for(int j=0;j<2;j++){
                    sum+=G[ix][iy][iz][p][i][j]*e[p][i][j];
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
                    sum+=G_new[ix][iy][iz][p][i][j]*b[p][i][j];
                }
            }
        }
    }
    else{
        for(int p=0;p<3;p++){
            for(int i=0; i<4;i++){
                for(int j=0;j<2;j++){
                    sum+=G[ix][iy][iz][p][i][j]*b[p][i][j];
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
    E0 = LatticeBoltzmann::E(ix,iy,iz,use_new); B0 = LatticeBoltzmann::B(ix,iy,iz,use_new);
    double rhoS = LatticeBoltzmann::rho_s(ix,iy,iz,s,use_new); double rhoSm1 = LatticeBoltzmann::rho_s(ix,iy,iz,(s+1)%2,use_new);
    VS = LatticeBoltzmann::V_s(ix,iy,iz,s,rhoS,use_new); VSm1 = LatticeBoltzmann::V_s(ix,iy,iz,(s+1)%2,rhoSm1,use_new);
    F0.load(0,0,0);//F0.load(rhoS*g,0,0);
    return (q[s]/m[s])*rhoS*(E0+(VS^B0)) - nu_s*rhoS*(VS-VSm1) + F0;
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
    return w[i]*rhoS*(3*xi[s]*pow(rhoS,Gamma-1) + 3*v[p][i]*Vmeds + 4.5*pow(v[p][i]*Vmeds,2) - 1.5*Vmeds.norm2());
}

double LatticeBoltzmann::feq0(double rhoS, vector3D &Vmeds, int s){
    return 3*rhoS*(1-0.5*(4*xi[s]*pow(rhoS,Gamma-1) + Vmeds.norm2()));
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
                                G[ixnew][iynew][iznew][p][i][s_j] = G_new[ix][iy][iz][p][i][s_j];
                                G_0[ix][iy][iz] = G_0_new[ix][iy][iz];
                            }
                            f[ixnew][iynew][iznew][p][i][s_j] = f_new[ix][iy][iz][p][i][s_j];
                            f_0[ix][iy][iz][s_j] = f_0_new[ix][iy][iz][s_j];
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
    for(int ix=0; ix<Lx;ix++){
        for(int iy=0;iy<Ly;iy++){
            for(int iz=0;iz<Lz;iz++){
                E0 = LatticeBoltzmann::E_med(ix,iy,iz,false); 
                B0 = LatticeBoltzmann::B(ix,iy,iz,false);  
                GEQ0 = LatticeBoltzmann::Geq0();
                G_0_new[ix][iy][iz] = G_0[ix][iy][iz] - (1/tau2)*(G_0[ix][iy][iz]-GEQ0);              
                 for(int p=0;p<3;p++){
                    for(int i=0; i<6;i++){
                        for(int s_j=0;s_j<2;s_j++){
                            if (i<4){
                               Jmed0 = LatticeBoltzmann::J_med(ix,iy,iz,false);
                               GEQ = LatticeBoltzmann::Geq(E0,B0,p,i,s_j);
                               G_new[ix][iy][iz][p][i][s_j] = G[ix][iy][iz][p][i][s_j] + ((2*tau2-1)/(16*tau2))*mu0*e[p][i][s_j]*Jmed0 - (1/tau2)*(G[ix][iy][iz][p][i][s_j]-GEQ);
                               
                            }
                            rhoS = LatticeBoltzmann::rho_s(ix,iy,iz,s_j,false);
                            Vmed0 = LatticeBoltzmann::V_s_med(ix,iy,iz,s_j,rhoS,false);
                            FEQS = LatticeBoltzmann::feq(rhoS,Vmed0,p,i,s_j);
                            FEQ0S = LatticeBoltzmann::feq0(rhoS, Vmed0,s_j);
                            Fs = LatticeBoltzmann::F_s(ix,iy,iz,s_j,false);
                            f_0_new[ix][iy][iz][s_j] = f_0[ix][iy][iz][s_j] - (1/taus[s_j])*(f_0[ix][iy][iz][s_j] - FEQ0S);
                            f_new[ix][iy][iz][p][i][s_j] = f[ix][iy][iz][p][i][s_j] + ((2*taus[s_j]-1)/(20*taus[s_j]))*v[p][i]*Fs  - (1/taus[s_j])*(f[ix][iy][iz][p][i][s_j]-FEQS);
                        }                        
                    }
                }
            }
        }
    }        
}

int main(void){
    LatticeBoltzmann a;
    a.Advection();
}