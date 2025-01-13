# include "utils.hpp"
# include "cubic.hpp"
# include "linear.hpp"

void build_original_function(Function * function, int Nx, int Ny, int Nz)
{
    float xi = -5.0f;
    float yi = -5.0f;
    float zi = -5.0f;

    float xf = 5.0f;
    float yf = 5.0f;
    float zf = 5.0f;

    float * x = new float[Nx]();
    float * y = new float[Ny]();
    float * z = new float[Nz]();

    linspace(x, xi, xf, Nx);
    linspace(y, yi, yf, Ny);
    linspace(z, zi, zf, Nz);

    for (int i = 0; i < Nz; i++)
    {
        for (int j = 0; j < Nx; j++)
        {
            for (int k = 0; k < Ny; k++)
            {
                function[i + j*Nz + k*Nx*Nz].z = z[i];
                function[i + j*Nz + k*Nx*Nz].x = x[j];
                function[i + j*Nz + k*Nx*Nz].y = y[k];

                function[i + j*Nz + k*Nx*Nz].v = sinf(z[i]) + sinf(x[j]) + sinf(y[k]);
            }
        }
    }

    delete[] x;
    delete[] y;
    delete[] z;
}

void perform_linear_interpolation(Function * F_in, Function * F_out, int skip_x, int skip_y, int skip_z, int Nx, int Ny, int Nz)
{
    float P[2][2][2];

    int rNx = (int)((Nx-1)/skip_x) + 1;
    int rNz = (int)((Nz-1)/skip_z) + 1;

    for (int i = skip_z+1; i < Nz-skip_z-1; i++)
    {   
        for (int j = skip_x+1; j < Nx-skip_x-1; j++)
        {   
            for (int k = skip_y+1; k < Ny-skip_y-1; k++)
            {   
                int idz = (int)(i/skip_z);
                int idx = (int)(j/skip_x);
                int idy = (int)(k/skip_y);

                float dz = (F_out[i + j*Nz + k*Nx*Nz].z - F_in[idz + idx*rNz + idy*rNx*rNz].z) / (F_in[(idz+1) + idx*rNz + idy*rNx*rNz].z - F_in[idz + idx*rNz + idy*rNx*rNz].z);  
                float dx = (F_out[i + j*Nz + k*Nx*Nz].x - F_in[idz + idx*rNz + idy*rNx*rNz].x) / (F_in[idz + (idx+1)*rNz + idy*rNx*rNz].x - F_in[idz + idx*rNz + idy*rNx*rNz].x);  
                float dy = (F_out[i + j*Nz + k*Nx*Nz].y - F_in[idz + idx*rNz + idy*rNx*rNz].y) / (F_in[idz + idx*rNz + (idy+1)*rNx*rNz].y - F_in[idz + idx*rNz + idy*rNx*rNz].y);  

                P[0][0][0] = F_in[idz + idx*rNz + idy*rNx*rNz].v;                
                P[0][0][1] = F_in[(idz+1) + idx*rNz + idy*rNx*rNz].v;                
                P[0][1][0] = F_in[idz + idx*rNz + (idy+1)*rNx*rNz].v;               
                P[0][1][1] = F_in[(idz+1) + idx*rNz + (idy+1)*rNx*rNz].v;

                P[1][0][0] = F_in[idz + (idx+1)*rNz + idy*rNx*rNz].v;                
                P[1][0][1] = F_in[(idz+1) + (idx+1)*rNz + idy*rNx*rNz].v;                
                P[1][1][0] = F_in[idz + (idx+1)*rNz + (idy+1)*rNx*rNz].v;               
                P[1][1][1] = F_in[(idz+1) + (idx+1)*rNz + (idy+1)*rNx*rNz].v;

                F_out[i + j*Nz + k*Nx*Nz].v = linear3D(P, dx, dy, dz);   
            }
        }
    }
}

void perform_cubic_interpolation(Function * F_in, Function * F_out, int skip_x, int skip_y, int skip_z, int Nx, int Ny, int Nz)
{
    float P[4][4][4];

    int rNx = (int)((Nx-1)/skip_x) + 1;
    int rNz = (int)((Nz-1)/skip_z) + 1;

    for (int i = skip_z+1; i < Nz-skip_z-1; i++)
    {   
        for (int j = skip_x+1; j < Nx-skip_x-1; j++)
        {   
            for (int k = skip_y+1; k < Ny-skip_y-1; k++)
            {   
                int idz = (int)(i/skip_z);
                int idx = (int)(j/skip_x);
                int idy = (int)(k/skip_y);

                float dz = (F_out[i + j*Nz + k*Nx*Nz].z - F_in[idz + idx*rNz + idy*rNx*rNz].z) / (F_in[(idz+1) + idx*rNz + idy*rNx*rNz].z - F_in[idz + idx*rNz + idy*rNx*rNz].z);  
                float dx = (F_out[i + j*Nz + k*Nx*Nz].x - F_in[idz + idx*rNz + idy*rNx*rNz].x) / (F_in[idz + (idx+1)*rNz + idy*rNx*rNz].x - F_in[idz + idx*rNz + idy*rNx*rNz].x);  
                float dy = (F_out[i + j*Nz + k*Nx*Nz].y - F_in[idz + idx*rNz + idy*rNx*rNz].y) / (F_in[idz + idx*rNz + (idy+1)*rNx*rNz].y - F_in[idz + idx*rNz + idy*rNx*rNz].y);  

                for (int pIdz = 0; pIdz < 4; pIdz++)
                {
                    for (int pIdx = 0; pIdx < 4; pIdx++)
                    {
                        for (int pIdy = 0; pIdy < 4; pIdy++)
                        {
                            P[pIdx][pIdy][pIdz] = F_in[(idz+pIdz-1) + (idx+pIdx-1)*rNz + (idy+pIdy-1)*rNx*rNz].v;                
                        }
                    }
                }
                
                F_out[i + j*Nz + k*Nx*Nz].v = cubic3D(P, dx, dy, dz);   
            }
        }
    }
}

int main()
{
    int original_Nx = 501;
    int original_Ny = 501;
    int original_Nz = 501;
    
    int rescaled_Nx = 26;
    int rescaled_Ny = 26;
    int rescaled_Nz = 26;

    int original_N = original_Nx*original_Ny*original_Nz;
    int rescaled_N = rescaled_Nx*rescaled_Ny*rescaled_Nz;

    auto * original_F = new Function[original_N]();
    auto * rescaled_F = new Function[rescaled_N]();

    build_original_function(original_F, original_Nx, original_Ny, original_Nz);

    int skip_x = (int)((original_Nx-1)/(rescaled_Nx-1));
    int skip_y = (int)((original_Ny-1)/(rescaled_Ny-1));
    int skip_z = (int)((original_Nz-1)/(rescaled_Nz-1));

    for (int i = 0; i < rescaled_Nz; i++)
    {
        for (int j = 0; j < rescaled_Nx; j++)
        {
            for (int k = 0; k < rescaled_Ny; k++)
            {
                int oId = i*skip_y + j*skip_x*original_Nz + k*skip_y*original_Nx*original_Nz;    
                int rId = i + j*rescaled_Nz + k*rescaled_Nx*rescaled_Nz;

                rescaled_F[rId].x = original_F[oId].x;
                rescaled_F[rId].y = original_F[oId].y;
                rescaled_F[rId].z = original_F[oId].z;
                rescaled_F[rId].v = original_F[oId].v;
            }
        }
    }

    auto * cubic_F = new Function[original_N]();
    auto * linear_F = new Function[original_N]();

    for (int n = 0; n < original_N; n++)
    {
        cubic_F[n].x = original_F[n].x;
        cubic_F[n].y = original_F[n].y;
        cubic_F[n].z = original_F[n].z;

        linear_F[n].x = original_F[n].x;
        linear_F[n].y = original_F[n].y;
        linear_F[n].z = original_F[n].z;
    }

    perform_cubic_interpolation(rescaled_F, cubic_F, skip_x, skip_y, skip_z, original_Nx, original_Ny, original_Nz);
    perform_linear_interpolation(rescaled_F, linear_F, skip_x, skip_y, skip_z, original_Nx, original_Ny, original_Nz);

    export_3d("outputs/cubic3D.bin", cubic_F, original_N);
    export_3d("outputs/linear3D.bin", linear_F, original_N);
    export_3d("outputs/original3D.bin", original_F, original_N);
    export_3d("outputs/rescaled3D.bin", rescaled_F, rescaled_N);

    return 0;
}