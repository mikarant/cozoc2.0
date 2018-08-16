static inline void omega_stencil(int i,int j,int k,
    PetscScalar hx,PetscScalar hy,PetscScalar hz,PetscScalar r,
    PetscScalar *f,PetscScalar *phi,const PetscScalar ***sigma,
    const PetscScalar ***zeta,const PetscScalar ****V,int *n,
    MatStencil *col,PetscScalar *w)
{
*n = 15;
col[ 0].k = k - 1 ; col[ 0].j = j     ; col[ 0].i = i - 1 ;
col[ 1].k = k     ; col[ 1].j = j     ; col[ 1].i = i - 1 ;
col[ 2].k = k + 1 ; col[ 2].j = j     ; col[ 2].i = i - 1 ;
col[ 3].k = k - 1 ; col[ 3].j = j - 1 ; col[ 3].i = i     ;
col[ 4].k = k     ; col[ 4].j = j - 1 ; col[ 4].i = i     ;
col[ 5].k = k + 1 ; col[ 5].j = j - 1 ; col[ 5].i = i     ;
col[ 6].k = k - 1 ; col[ 6].j = j     ; col[ 6].i = i     ;
col[ 7].k = k     ; col[ 7].j = j     ; col[ 7].i = i     ;
col[ 8].k = k + 1 ; col[ 8].j = j     ; col[ 8].i = i     ;
col[ 9].k = k - 1 ; col[ 9].j = j + 1 ; col[ 9].i = i     ;
col[10].k = k     ; col[10].j = j + 1 ; col[10].i = i     ;
col[11].k = k + 1 ; col[11].j = j + 1 ; col[11].i = i     ;
col[12].k = k - 1 ; col[12].j = j     ; col[12].i = i + 1 ;
col[13].k = k     ; col[13].j = j     ; col[13].i = i + 1 ;
col[14].k = k + 1 ; col[14].j = j     ; col[14].i = i + 1 ;
w[ 0] = 0.125*hy*r*f[j]*V[k-1][j][i][1]/hz - 0.125*hy*r*f[j]*V[k+1][j][i][1]/hz;
w[ 1] = 0.5*hy*r*f[j]*V[k-1][j][i][1]/hz + 0.5*hy*r*f[j]*V[k+1][j][i][1]/hz - 1.0*hy*r*f[j]*V[k][j][i][1]/hz + 1.0*hy*hz*sigma[k][j][i-1]/(hx*cos(phi[j]));
w[ 2] = -0.125*hy*r*f[j]*V[k-1][j][i][1]/hz + 0.125*hy*r*f[j]*V[k+1][j][i][1]/hz;
w[ 3] = -1.0L/8.0L*hx*r*cos(phi[j])*f[j]*V[k-1][j][i][0]/hz + (1.0L/8.0L)*hx*r*cos(phi[j])*f[j]*V[k+1][j][i][0]/hz;
w[ 4] = 0.5*hx*hz*cos(phi[j])*tan(phi[j])*sigma[k][j-1][i] - 1.0L/2.0L*hx*r*cos(phi[j])*f[j]*V[k-1][j][i][0]/hz - 1.0L/2.0L*hx*r*cos(phi[j])*f[j]*V[k+1][j][i][0]/hz + hx*r*cos(phi[j])*f[j]*V[k][j][i][0]/hz + 1.0*hx*hz*cos(phi[j])*sigma[k][j-1][i]/hy;
w[ 5] = (1.0L/8.0L)*hx*r*cos(phi[j])*f[j]*V[k-1][j][i][0]/hz - 1.0L/8.0L*hx*r*cos(phi[j])*f[j]*V[k+1][j][i][0]/hz;
w[ 6] = hx*hy*r*r*cos(phi[j])*f[j]*f[j]/hz + hx*hy*r*r*cos(phi[j])*f[j]*zeta[k][j][i]/hz;
w[ 7] = -2*hx*hy*r*r*cos(phi[j])*f[j]*f[j]/hz - hx*hy*r*r*cos(phi[j])*f[j]*zeta[k-1][j][i]/hz - hx*hy*r*r*cos(phi[j])*f[j]*zeta[k+1][j][i]/hz - 2.0*hx*hz*cos(phi[j])*sigma[k][j][i]/hy - 2.0*hy*hz*sigma[k][j][i]/(hx*cos(phi[j]));
w[ 8] = hx*hy*r*r*cos(phi[j])*f[j]*f[j]/hz + hx*hy*r*r*cos(phi[j])*f[j]*zeta[k][j][i]/hz;
w[ 9] = (1.0L/8.0L)*hx*r*cos(phi[j])*f[j]*V[k-1][j][i][0]/hz - 1.0L/8.0L*hx*r*cos(phi[j])*f[j]*V[k+1][j][i][0]/hz;
w[10] = -0.5*hx*hz*cos(phi[j])*tan(phi[j])*sigma[k][j+1][i] + (1.0L/2.0L)*hx*r*cos(phi[j])*f[j]*V[k-1][j][i][0]/hz + (1.0L/2.0L)*hx*r*cos(phi[j])*f[j]*V[k+1][j][i][0]/hz - hx*r*cos(phi[j])*f[j]*V[k][j][i][0]/hz + 1.0*hx*hz*cos(phi[j])*sigma[k][j+1][i]/hy;
w[11] = -1.0L/8.0L*hx*r*cos(phi[j])*f[j]*V[k-1][j][i][0]/hz + (1.0L/8.0L)*hx*r*cos(phi[j])*f[j]*V[k+1][j][i][0]/hz;
w[12] = -0.125*hy*r*f[j]*V[k-1][j][i][1]/hz + 0.125*hy*r*f[j]*V[k+1][j][i][1]/hz;
w[13] = -0.5*hy*r*f[j]*V[k-1][j][i][1]/hz - 0.5*hy*r*f[j]*V[k+1][j][i][1]/hz + 1.0*hy*r*f[j]*V[k][j][i][1]/hz + 1.0*hy*hz*sigma[k][j][i+1]/(hx*cos(phi[j]));
w[14] = 0.125*hy*r*f[j]*V[k-1][j][i][1]/hz - 0.125*hy*r*f[j]*V[k+1][j][i][1]/hz;
}
