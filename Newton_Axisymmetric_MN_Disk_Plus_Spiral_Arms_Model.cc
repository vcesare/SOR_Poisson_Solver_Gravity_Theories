/*
  © 2025 Valentina Cesare
  Licensed under the MIT License.
  See LICENSE file for full license text.
*/

// This code finds the Newtonian gravitational potential, φN, given a mass density distribution rho, by numerically solving Newtonian Poisson Equation (∇2φN = 4πGρ; e.g., Eq. (3) in Cesare, et al. (2020) (https://ui.adsabs.harvard.edu/abs/2020A%26A...637A..70C/abstract) over a 3D grid in cylindrical coordinates, (R,z,phi), with Successive-Over-Relaxation (SOR) Poisson Solver. This Poisson Solver is a 3D extension of 2D Poisson Solver described in Appendix B of Cesare, et al. (2020).
// Author: Valentina Cesare.
// E-Mail address: valentina.cesare@inaf.it
// Creation date: July 31th 2025.

// The Poisson Solver for the density distribution considered in this code converges in 2599 iterations.

#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <cmath>
#include <random>
#include <iomanip>
#include <fstream>

using namespace std;


double sech(double x) {
    return 1.0 / cosh(x);
}

double smooth_cutoff(double R, double Rmin_cutoff) {
    return 1.0 - exp(-R*R/(Rmin_cutoff*Rmin_cutoff));
}

double Error(double ***V, double ***V_New);


//    Number of grid points
#define NR 301
#define NZ 300
#define NPHI 100

int main()
{
//    pi
    const double CONST_PI = 3.14159265358979;
// Universal gravitational constant is in (kpc/Msun)*(km/s)^2
    const double G = 4.3011e-6;
    
//    The grid along R and z is in kpc.
    double xbeg = -30.0, xend = 30.0;
    double zbeg = -30.0, zend = 30.0;
//    The grid along phi is in radians.
    double phibeg = 0.0, phiend = 2.0*CONST_PI;
    
    double dx = (xend-xbeg)/double(NR);
    double dz = (zend-zbeg)/double(NZ);
    double dphi = (phiend-phibeg)/double(NPHI);
    
    double dx2 = dx*dx;
    double dz2 = dz*dz;
    double dphi2 = dphi*dphi;
    
    double x[NR + 1], z[NZ + 1], phi[NPHI];
    
    int Nmin;
    Nmin = NR;
    if (Nmin > NZ)
        Nmin = NZ;
        
    if (Nmin > NPHI)
        Nmin = NPHI;
    
    cout<<"Nmin = "<<Nmin<<endl;
    
//    Parameter that regulates the convergence and the speed of SOR Poisson Solver.
    double wsor = 2.0/(1.0 + CONST_PI/Nmin);
    cout<<"wsor = "<<wsor<<endl;
    
    // Initialize grid (x and z grids are in kpc, phi grid is in radians).
    for (int i = 0; i <= NR; i++) x[i] = xbeg + i*dx;
    for (int k = 0; k <= NZ; k++) z[k] = zbeg + k*dz;
    for (int j = 0; j < NPHI; j++) phi[j] = phibeg + j*dphi;
    
    for (int i = 0; i <= NR; i++)
        cout<<"x["<<i<<"] = "<<x[i]<<endl;
    
    for (int k = 0; k <= NZ; k++)
        cout<<"z["<<k<<"] = "<<z[k]<<endl;
    
    for (int j = 0; j < NPHI; j++)
        cout<<"phi["<<j<<"] = "<<phi[j]<<endl;
    
    double R[NR + 1];
    for (int i = 0; i <= NR; i++) R[i] = abs(x[i]);
    
//    Files writing are for plots (if you do not need them, you can comment them out).
    ofstream outfileRHalf ("R_Half.txt");
    {
        for(int i = (NR - 1)/2 + 1; i <= NR; i++) outfileRHalf<<R[i]<<endl;
    }
    outfileRHalf.close();
    
    ofstream outfileR ("R.txt");
    {
        for(int i = 0; i <= NR; i++) outfileR<<R[i]<<endl;
    }
    outfileR.close();
    
    ofstream outfilez ("z.txt");
    {
        for(int k = 0; k <= NZ; k++) outfilez<<z[k]<<endl;
    }
    outfilez.close();
    
    ofstream outfilephi ("phi.txt");
    {
        for(int j = 0; j < NPHI; j++) outfilephi<<phi[j]<<endl;
    }
    outfilephi.close();
    
    double dR = dx;
    double dR2 = dx2;
    
    double **rhoAxisymm, **rhoPertAmplitude, **phiNAxisymmTheo;
    rhoAxisymm=new double*[NR+1];
    rhoPertAmplitude=new double*[NR+1];
    phiNAxisymmTheo=new double*[NR+1];
    for(int i = 0; i <= NR; i++)
    {
        rhoAxisymm[i]=new double[NZ+1];
        rhoPertAmplitude[i]=new double[NZ+1];
        phiNAxisymmTheo[i]=new double[NZ+1];
    }
    
    double **gamma;
    gamma=new double*[NR+1];
    for(int i = 0; i <= NR; i++)
    {
        gamma[i]=new double[NPHI+1];
    }
    
    double ***rhoPert, ***rho, ***SN, ***phiNPertTheo, ***phiNTheo, ***phi0N, ***phi1N;
    rhoPert=new double**[NR+1];
    rho=new double**[NR+1];
    SN=new double**[NR+1];
    phiNPertTheo=new double**[NR+1];
    phiNTheo=new double**[NR+1];
    phi0N=new double**[NR+1];
    phi1N=new double**[NR+1];
    for(int i = 0; i <= NR; i++)
    {
        rhoPert[i]=new double*[NZ+1];
        rho[i]=new double*[NZ+1];
        SN[i]=new double*[NZ+1];
        phiNPertTheo[i]=new double*[NZ+1];
        phiNTheo[i]=new double*[NZ+1];
        phi0N[i]=new double*[NZ+1];
        phi1N[i]=new double*[NZ+1];
        for(int k = 0; k <= NZ; k++)
        {
            rhoPert[i][k]=new double[NPHI+1];
            rho[i][k]=new double[NPHI+1];
            SN[i][k]=new double[NPHI+1];
            phiNPertTheo[i][k]=new double[NPHI+1];
            phiNTheo[i][k]=new double[NPHI+1];
            phi0N[i][k]=new double[NPHI+1];
            phi1N[i][k]=new double[NPHI+1];
        }
    }
    
//    As axisymmetric density distribution, a Miyamoto-Nagai disk is considered.
    
//    Masses are in solar masses.
    double MMN = 1.0e10;
    
//    Distances are in kpc.
    double aMN = 3.0;
    double bMN = 1.0;
    
//    Mass density of a Miyamoto-Nagai disk (Eq. (B.12) in Cesare, et al. (2020)).
    for (int i = 0; i <= NR; i++)
    {
        for (int k = 0; k <= NZ; k++)
        {
            rhoAxisymm[i][k] = (MMN*bMN*bMN*(aMN*R[i]*R[i] + (aMN + 3*sqrt(z[k]*z[k] + bMN*bMN))*(aMN + sqrt(z[k]*z[k] + bMN*bMN))*(aMN + sqrt(z[k]*z[k] + bMN*bMN))))/(4.0*CONST_PI*pow(R[i]*R[i] + (aMN + sqrt(z[k]*z[k] + bMN*bMN))*(aMN + sqrt(z[k]*z[k] + bMN*bMN)),2.5)*pow(z[k]*z[k] + bMN*bMN,1.5));
        }
    }
    
//    These numbers refer to the first example in Section 4 of "Cox and Gomez (2002)" (https://ui.adsabs.harvard.edu/abs/2002ApJS..142..261C/abstract) (*).
//    Number of spiral arms.
    int N = 3;  // Only exception from the considered example, where N = 2.
//    Pitch angle (in radians)
    double alpha = CONST_PI/12;  // It corresponds to alpha = 15°.
//    Radial scale length of the dropoff in density amplitude of the arms (in kpc).
    double Rs = 7.0;
//    Fiducial radius (in kpc).
    double r0 = 8.0;
//    Midplane arm number density at fiducial radius r0 (in atom/kpc^3).
    double n0 = 3.086*3.086*3.086e54;  // It corresponds to n0 = 1 atom/cm^3.
//    Hydrogen mass (in solar masses).
    double mH = 1.0e-57;
    mH = (1.6735575/1.989)*mH;  // It corresponds to 1,6735575 × 10^-27 kg.
//    Average mass per atom (in solar masses).
    double m = (14/11)*mH;
//    Midplane arm mass density at fiducial radius r0 (in solar masses/kpc^3).
    double rho0 = m*n0;
//    Scale height of the stellar arm perturbation (in kpc).
    double H = 0.18;  // Chosen to match the scaleheight of the thin stellar disk of Dehnen and Binney, 1998, Model 2 (*).
//    phip(r0) in Eq. (3) in Cox & Gomez, 2002 is assumed to be equal to 0 radians.
    double phipr0 = 0.0;
    
//    Amplitude of the spiral density distribution (Eq. (1) in Cox and Gomez (2002)):
    for (int i = 0; i <= NR; i++)
    {
        for (int k = 0; k <= NZ; k++)
        {
            rhoPertAmplitude[i][k] = rho0*exp(-(R[i] - r0)/Rs)*sech(z[k]/H)*sech(z[k]/H);
        }
    }
    
    
//    gamma (Eq. (3) in Cox & Gomex, 2002).
    for (int i = 0; i <= NR; i++)
    {
        for (int j = 0; j < NPHI; j++)
        {
            gamma[i][j] = N*(phi[j] - phipr0 - log(R[i]/r0)/tan(alpha));
        }
    }
    
//    Given that the number of spiral arms is N = 3, the equation for the mass density is given by Eq. (4) in Cox & Gomez (2002), summing over n from 1 to N = 3, with C1 = 8/(3π), C2 = 1/2, and C3 = 8/(15π), as defined below Eq. (4).
    double C1 = 8.0/(3.0*CONST_PI);
    double C2 = 0.5;
    double C3 = 8.0/(15.0*CONST_PI);
    
    for (int i = 0; i <= NR; i++)
    {
        for (int k = 0; k <= NZ; k++)
        {
            for (int j = 0; j < NPHI; j++)
            {
                rhoPert[i][k][j] = rhoPertAmplitude[i][k]*(C1*cos(1.0*gamma[i][j]) + C2*cos(2.0*gamma[i][j]) + C3*cos(3.0*gamma[i][j]));
            }
        }
    }
    
// OPTIONAL:
//    Introduction of a smooth cutoff in rhoPertModulation to remove the divergent behaviour of log(R) and make the system physical close to the galaxy center.
//    Phisical interpretation of smooth cutoff:
//    Function smooth_cutoff is:
//    smooth_cutoff(R) -> 0 for R -> 0,
//    which means that the perturbative density slowly gets to zero close to the galaxy center;
//    smooth_cutoff(R) -> 1 for R >> Rmin_cutoff,
//    which means that the perturbation gets full amplitude.
//    Physical sense of this smooth cutoff in rhoPertModulation:
//    (1) The spiral arms do not appear from a ray defined "magically": their intensity grows gradually with the ray.
//    (2) There is a central bulge often dominated by different (symmetrical) dynamics where spirals do not develop.
//    (3) This function naturally modulates the amplitude of the spirals, respecting this transition.
    double Rmin_cutoff = 2.0; // kpc
    
    for (int i = 0; i <= NR; i++)
        cout<<"R["<<i<<"] = "<<R[i]<<" smooth_cutoff = "<<smooth_cutoff(R[i],Rmin_cutoff)<<endl;
    
    for (int i = 0; i <= NR; i++)
    {
        for (int k = 0; k <= NZ; k++)
        {
            for (int j = 0; j < NPHI; j++)
            {
                rhoPert[i][k][j] *= smooth_cutoff(R[i],Rmin_cutoff);
            }
        }
    }
    
//    Total mass density, given by the sum of the Miyamoto-Nagai disk mass density (Eq. (B.12) in Cesare, et al., 2020) and the perturbed density (Eq. (4) in Cox & Gomex, 2002).
    for (int i = 0; i <= NR; i++)
    {
        for (int k = 0; k <= NZ; k++)
        {
            for (int j = 0; j < NPHI; j++)
            {
                rho[i][k][j] = rhoAxisymm[i][k] + rhoPert[i][k][j];
            }
        }
    }
    
//    Newtonian source term of SOR Poisson Solver (Eq. (B.6) in Cesare, et al. (2020); https://ui.adsabs.harvard.edu/abs/2020A%26A...637A..70C/abstract).
    for (int i = 0; i <= NR; i++)
    {
        for (int k = 0; k <= NZ; k++)
        {
            for (int j = 0; j < NPHI; j++)
            {
                SN[i][k][j] = -4.0*CONST_PI*G*rho[i][k][j];
            }
        }
    }
    
    
    
//    ANALYTICAL NEWTONIAN POTENTIAL:
    
//    Theoretical Newtonian Miyamoto-Nagai axisymmetric potential, given by Eq. (B.13) in Cesare, et al. (2020):
    for (int i = 0; i <= NR; i++)
    {
        for (int k = 0; k <= NZ; k++)
        {
            phiNAxisymmTheo[i][k] = -G*MMN/sqrt(R[i]*R[i] + (aMN + sqrt(z[k]*z[k] + bMN*bMN))*(aMN + sqrt(z[k]*z[k] + bMN*bMN)));
        }
    }
    
//    Calculation of theoretical Newtonian potential generated by the spiral pattern for Poisson Solver boundary conditions, for N = 3 spiral arms.
    
    
    double K1[NR + 1], K2[NR + 1], K3[NR + 1];
    double beta1[NR + 1], beta2[NR + 1], beta3[NR + 1];
    double D1[NR + 1], D2[NR + 1], D3[NR + 1];
    
//    Eq. (5) in Cox & Gomez (2002).
    for (int i = 0; i <= NR; i++)
    {
        K1[i] = 1.0*N/(R[i]*sin(alpha));
        K2[i] = 2.0*N/(R[i]*sin(alpha));
        K3[i] = 3.0*N/(R[i]*sin(alpha));
    }
//    Eq. (6) in Cox & Gomez (2002).
    for (int i = 0; i <= NR; i++)
    {
        beta1[i] = K1[i]*H*(1.0 + 0.4*K1[i]*H);
        beta2[i] = K2[i]*H*(1.0 + 0.4*K2[i]*H);
        beta3[i] = K3[i]*H*(1.0 + 0.4*K3[i]*H);
    }
//    Eq. (7) in Cox & Gomez (2002).
    for (int i = 0; i <= NR; i++)
    {
        D1[i] = (1.0 + K1[i]*H + 0.3*(K1[i]*H)*(K1[i]*H))/(1.0 + 0.3*K1[i]*H);
        D2[i] = (1.0 + K2[i]*H + 0.3*(K2[i]*H)*(K2[i]*H))/(1.0 + 0.3*K2[i]*H);
        D3[i] = (1.0 + K3[i]*H + 0.3*(K3[i]*H)*(K3[i]*H))/(1.0 + 0.3*K3[i]*H);
    }
    
//    Equation (8) in Cox & Gomez (2002), summing over n from 1 to N = 3.
    for (int i = 0; i <= NR; i++)
    {
        for (int k = 0; k <= NZ; k++)
        {
            for (int j = 0; j < NPHI; j++)
            {
                phiNPertTheo[i][k][j] = -4.0*CONST_PI*G*H*rho0*exp(-(R[i] - r0)/Rs)*(
                    (C1/(K1[i]*D1[i]))*cos(1.0*gamma[i][j]*pow(sech(K1[i]*z[k]/beta1[i]),beta1[i])) +
                    (C2/(K2[i]*D2[i]))*cos(2.0*gamma[i][j]*pow(sech(K2[i]*z[k]/beta2[i]),beta2[i])) +
                    (C3/(K3[i]*D3[i]))*cos(3.0*gamma[i][j]*pow(sech(K3[i]*z[k]/beta3[i]),beta3[i]))
                                                                );
            }
        }
    }
    
    
//    OPTIONAL:
//    Application of the same smooth_cutoff() function as for rhoPertModulation:
    for (int i = 0; i <= NR; i++)
    {
        for (int k = 0; k <= NZ; k++)
        {
            for (int j = 0; j < NPHI; j++)
            {
                phiNPertTheo[i][k][j] *= smooth_cutoff(R[i], Rmin_cutoff);
            }
        }
    }
    
//    Full (axisymmetric Miyamoto-Nagai -> Eq. (B.13 in Cesare, et al., 2020) + spiral arms -> Eq. (8) in Cox & Gomez, 2002) theoretical Newtonian potential.
    for (int i = 0; i <= NR; i++)
    {
        for (int k = 0; k <= NZ; k++)
        {
            for (int j = 0; j < NPHI; j++)
            {
                phiNTheo[i][k][j] = phiNAxisymmTheo[i][k] + phiNPertTheo[i][k][j];
            }
        }
    }
    
    
//  FILE WRITING FOR ANALYTICAL NEWTONIAN POTENTIAL.
    
//    The three files writing below are for plots (you can comment them out if you do not need them).
    ofstream outfilephiNTheoR ("phiNTheo_MNPlusSpiral_R_zeq0_phieq0_Smooth_Cutoff.txt");
    {
        for(int i = (NR - 1)/2 + 1; i <= NR; i++) outfilephiNTheoR<<phiNTheo[i][NZ/2][0]<<endl;
    }
    outfilephiNTheoR.close();
    
    ofstream outfilephiNTheoz ("phiNTheo_MNPlusSpiral_z_Rsimeq0_phieq0_Smooth_Cutoff.txt");
    {
        for(int k = 0; k <= NZ; k++) outfilephiNTheoz<<phiNTheo[(NR - 1)/2 + 1][k][0]<<endl;
    }
    outfilephiNTheoz.close();
    
    ofstream outfilephiNTheophi ("phiNTheo_MNPlusSpiral_phi_Rsimeq0_zeq0_Smooth_Cutoff.txt");
    {
        for(int j = 0; j < NPHI; j++) outfilephiNTheophi<<phiNTheo[(NR - 1)/2 + 1][NZ/2][j]<<endl;
    }
    outfilephiNTheophi.close();
    
//    File writing of full (axisymmetric Miyamoto-Nagai + spiral arms) theoretical Newtonian potential. -> NECESSARY AS INPUT FOR QUMOND POISSON SOLVER (SEE CODE "QUMOND_Axisymmetric_MN_Disk_Plus_Spiral_Arms_Model.cc").
    ofstream outfilephiNTheo ("phiNTheo_MNPlusSpiral_Smooth_Cutoff.txt");
    {
        for(int i = 0; i <= NR; i++)
        {
            for(int k = 0; k <= NZ; k++)
            {
                for(int j = 0; j < NPHI; j++)
                {
                    outfilephiNTheo<<phiNTheo[i][k][j]<<endl;
                }
            }
        }
    }
    outfilephiNTheo.close();
    
//    exit(0);  // PUT AN exit (0); HERE IF YOU ARE NOT INTERESTED IN EXECUTING THE NEWTONIAN POISSON SOLVER BUT YOU ONLY WANT TO WRITE TO FILE THE FULL THEORETICAL NEWTONIAN POTENTIAL TO GIVE AS INPUT TO QUMOND POISSON SOLVER!!!
    
//    NUMERICAL NEWTONIAN POTENTIAL:
    
//    POISSON SOLVER BASED ON SOR METHOD:
    
    int counter = 0;
    double errN = 1.0e30;
    double tol = 1.0e-9;
    
//    Newtonian potential inizialization.
    for (int i = 0; i <= NR; i++)
    {
        for (int k = 0; k <= NZ; k++)
        {
            for (int j = 0; j < NPHI; j++)
            {
                phi0N[i][k][j] = 0.0;
            }
        }
    }
    
//    Dirichlet boundary conditions for Newtonian potential
    // Bottom
    for (int i = 0; i <= NR; i++)
    {
        for (int k = 0; k <= NZ; k++)
        {
            phi0N[i][k][0] = phiNTheo[i][k][0];
        }
    }
    
//    !!! THIS BOUNDARY CONDITION IS NO MORE NECESSARY SINCE THE POTENTIAL IS PERIODIC ALONG THE PHI-DIRECTION OF THE DOMAIN.
//    // Top
//    for (int i = 0; i <= NR; i++)
//    {
//        for (int k = 0; k <= NZ; k++)
//        {
//            phi0N[i][k][NPHI - 1] = phiNTheo[i][k][NPHI - 1];
//        }
//    }
    
    // Periodicity boundary condition
    for (int i = 0; i <= NR; i++) {
        for (int k = 0; k <= NZ; k++) {
            phi0N[i][k][NPHI - 1] = phi0N[i][k][0];
        }
    }
    
    // Left
    for (int k = 0; k <= NZ; k++)
    {
        for (int j = 0; j < NPHI; j++)
        {
            phi0N[0][k][j] = phiNTheo[0][k][j];
        }
    }
    
    // Right
    for (int k = 0; k <= NZ; k++)
    {
        for (int j = 0; j < NPHI; j++)
        {
            phi0N[NR][k][j] = phiNTheo[NR][k][j];
        }
    }
    
    // Back
    for (int i = 0; i <= NR; i++)
    {
        for (int j = 0; j < NPHI; j++)
        {
            phi0N[i][0][j] = phiNTheo[i][0][j];
        }
    }
    
    // Front
    for (int i = 0; i <= NR; i++)
    {
        for (int j = 0; j < NPHI; j++)
        {
            phi0N[i][NZ][j] = phiNTheo[i][NZ][j];
        }
    }
    
//    SOR POISSON SOLVER WHILE LOOP UP TO CONVERGENCE:
    
    while (errN > tol) {
        for(int i = 0; i <= NR; i++)
        {
            for(int k = 0; k <= NZ; k++)
            {
                for (int j = 0; j < NPHI; j++)
                {
                    phi1N[i][k][j] = phi0N[i][k][j];
                }
            }
        }
        
        for(int i = 1; i < NR; i++)
        {
            for(int k = 1; k < NZ; k++)
            {
                for (int j = 1; j < NPHI - 1; j++)
                {
                    int jplus  = (j + 1) % NPHI;
                    int jminus = (j - 1 + NPHI) % NPHI;
                    phi1N[i][k][j] = phi0N[i][k][j]*(1.0 - wsor) + wsor/(2.0 + 2.0*dR2/dz2 + 2.0*dR2/(R[i]*R[i]*dphi2))*((1.0 + dR/(2.0*R[i]))*phi1N[i+1][k][j] + (1.0 - dR/(2.0*R[i]))*phi1N[i-1][k][j] + dR2/(dphi2*R[i]*R[i])*(phi1N[i][k][jplus] + phi1N[i][k][jminus]) + (dR2/dz2)*(phi1N[i][k+1][j] + phi1N[i][k-1][j]) + dR2*SN[i][k][j]);
                }
            }
        }
        
        // Compute error.
        errN = Error(phi0N, phi1N);
        
        for(int i = 1; i < NR; i++)
        {
            for(int k = 1; k < NZ; k++)
            {
                for (int j = 1; j < NPHI - 1; j++)
                {
                    phi0N[i][k][j] = phi1N[i][k][j];
                }
            }
        }
        counter++;
        cout <<"counter = "<< counter << "; errN = " << setprecision(8) << errN << endl;
    }
    
//    FILE WRITING FOR NUMERICAL NEWTONIAN POTENTIAL.
//    Files writing are for plots (if you do not need them, you can comment them out).
    ofstream outfilephiNCompR ("phiNComp_MNPlusSpiral_R_zeq0_phieq0_Smooth_Cutoff.txt");
    {
        for(int i = (NR - 1)/2 + 1; i <= NR; i++) outfilephiNCompR<<phi1N[i][NZ/2][0]<<endl;
    }
    outfilephiNCompR.close();
    
    ofstream outfilephiNCompz ("phiNComp_MNPlusSpiral_z_Rsimeq0_phieq0_Smooth_Cutoff.txt");
    {
        for(int k = 0; k <= NZ; k++) outfilephiNCompz<<phi1N[(NR - 1)/2 + 1][k][0]<<endl;
    }
    outfilephiNCompz.close();
    
    ofstream outfilephiNCompphi ("phiNComp_MNPlusSpiral_phi_Rsimeq0_zeq0_Smooth_Cutoff.txt");
    {
        for(int j = 0; j < NPHI; j++) outfilephiNCompphi<<phi1N[(NR - 1)/2 + 1][NZ/2][j]<<endl;
    }
    outfilephiNCompphi.close();
    
    ofstream outfilephiNComp ("phiNComp_MNPlusSpiral_Smooth_Cutoff.txt");
    {
        for(int i = 0; i <= NR; i++)
        {
            for(int k = 0; k <= NZ; k++)
            {
                for(int j = 0; j < NPHI; j++)
                {
                    outfilephiNComp<<phi1N[i][k][j]<<endl;
                }
            }
        }
    }
    outfilephiNComp.close();
    
    
//    MEMORY DEALLOCATION:
    for(int i = 0; i <= NR; i++)
    {
        delete[] rhoAxisymm[i];
        delete[] rhoPertAmplitude[i];
        delete[] phiNAxisymmTheo[i];
    }
    delete[] rhoAxisymm;
    delete[] rhoPertAmplitude;
    delete[] phiNAxisymmTheo;
    
    for(int i = 0; i <= NR; i++)
    {
        delete[] gamma[i];
    }
    delete[] gamma;
    
    for(int i = 0; i <= NR; i++)
    {
        for(int k = 0; k <= NZ; k++)
        {
            delete[] rhoPert[i][k];
            delete[] rho[i][k];
            delete[] SN[i][k];
            delete[] phiNPertTheo[i][k];
            delete[] phiNTheo[i][k];
            delete[] phi0N[i][k];
            delete[] phi1N[i][k];
        }
    }
    
    delete[] rhoPert;
    delete[] rho;
    delete[] SN;
    delete[] phiNPertTheo;
    delete[] phiNTheo;
    delete[] phi0N;
    delete[] phi1N;
    
    
    
    return 0;
}

////////////////////////////////////////////////////////////////////////
double Error(double ***V0, double ***V1)
//
// Compute error
////////////////////////////////////////////////////////////////////////
{
    double err;
    double err_L1  = 0.0;   // average relative error per lattice point
    double err_max = 0.0;
    
    for (int i = 0; i <= NR; i++)
    {
        for (int k = 0; k <= NZ; k++)
        {
            for (int j = 0; j < NPHI; j++)
            {
                err      = fabs(V1[i][k][j] - V0[i][k][j]); //(fabs(V0[i][j]) + 1.e-40);
                err_L1  += err;
                err_max  = (err > err_max ? err:err_max);
            }
        }
    }
    
    err_L1 /= (double)(NR - 1)*(NZ - 1)*(NPHI - 1);
    return err_L1;
}
////////////////////////////////////////////////////////////////////////
