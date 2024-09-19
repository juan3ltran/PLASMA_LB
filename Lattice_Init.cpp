#include "Lattice.hpp"
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

