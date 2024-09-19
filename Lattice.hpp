#include "constants.hpp"
#include "vector.hpp"

#ifndef LATTICE_H
#define LATTICE_H

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
        void Print(void);
        
};

#endif //LATTICE_H

