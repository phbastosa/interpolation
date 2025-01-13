# include "utils.hpp"
# include "cubic.hpp"
# include "linear.hpp"

void build_original_function(Function * function, int Nx, int Ny)
{
    float xi = -10.0f;
    float yi = -10.0f;

    float xf = 10.0f;
    float yf = 10.0f;

    float * x = new float[Nx]();
    float * y = new float[Ny]();

    linspace(x, xi, xf, Nx);
    linspace(y, yi, yf, Ny);

    for (int i = 0; i < Ny; i++)
    {
        for (int j = 0; j < Nx; j++)
        {
            function[i + j*Ny].x = x[j];
            function[i + j*Ny].y = y[i];

            function[i + j*Ny].v = sinf(y[i]) + sinf(x[j]);
        }
    }
}

void perform_linear_interpolation(Function * F_in, Function * F_out, int skip_x, int skip_y, int Nx, int Ny)
{
    float P[2][2];

    int rNy = (int)((Ny-1)/skip_y) + 1;

    for (int i = skip_y+1; i < Ny-skip_y-1; i++)
    {   
        for (int j = skip_x+1; j < Nx-skip_x-1; j++)
        {   
            int idx = (int)(j/skip_x);
            int idy = (int)(i/skip_y);

            float dx = (F_out[i + j*Ny].x - F_in[idy + idx*rNy].x) / (F_in[idy + (idx+1)*rNy].x - F_in[idy + idx*rNy].x);  
            float dy = (F_out[i + j*Ny].y - F_in[idy + idx*rNy].y) / (F_in[(idy+1) + idx*rNy].y - F_in[idy + idx*rNy].y);  

            P[0][0] = F_in[idy + idx*rNy].v;                
            P[0][1] = F_in[(idy+1) + idx*rNy].v;                
            P[1][0] = F_in[idy + (idx+1)*rNy].v;               
            P[1][1] = F_in[(idy+1) + (idx+1)*rNy].v;
            
            F_out[i + j*Ny].v = linear2D(P, dx, dy);   
        }
    }
}

void perform_cubic_interpolation(Function * F_in, Function * F_out, int skip_x, int skip_y, int Nx, int Ny)
{
    float P[4][4];

    int rNy = (int)((Ny-1)/skip_y) + 1;

    for (int i = skip_y+1; i < Ny-skip_y-1; i++)
    {   
        for (int j = skip_x+1; j < Nx-skip_x-1; j++)
        {   
            int idx = (int)(j/skip_x);
            int idy = (int)(i/skip_y);

            float dx = (F_out[i + j*Ny].x - F_in[idy + idx*rNy].x) / (F_in[idy + (idx+1)*rNy].x - F_in[idy + idx*rNy].x);  
            float dy = (F_out[i + j*Ny].y - F_in[idy + idx*rNy].y) / (F_in[(idy+1) + idx*rNy].y - F_in[idy + idx*rNy].y);  

            for (int pIdx = 0; pIdx < 4; pIdx++)
            {
                for (int pIdy = 0; pIdy < 4; pIdy++)
                {
                    P[pIdx][pIdy] = F_in[(idy+pIdy-1) + (idx+pIdx-1)*rNy].v;                
                }
            }

            F_out[i + j*Ny].v = cubic2D(P, dx, dy);   
        }
    }
}

int main()
{
    int original_Nx = 1001;
    int original_Ny = 1001;
    
    int rescaled_Nx = 51;
    int rescaled_Ny = 51;

    int original_N = original_Nx*original_Ny;
    int rescaled_N = rescaled_Nx*rescaled_Ny;

    auto * original_F = new Function[original_N]();
    auto * rescaled_F = new Function[rescaled_N]();

    build_original_function(original_F, original_Nx, original_Ny);

    int skip_x = (int)((original_Nx-1)/(rescaled_Nx-1));
    int skip_y = (int)((original_Ny-1)/(rescaled_Ny-1));

    for (int i = 0; i < rescaled_Ny; i++)
    {
        for (int j = 0; j < rescaled_Nx; j++)
        {
            int index = i*skip_y + j*skip_x*original_Ny;    

            rescaled_F[i + j*rescaled_Ny].x = original_F[index].x;
            rescaled_F[i + j*rescaled_Ny].y = original_F[index].y;
            rescaled_F[i + j*rescaled_Ny].v = original_F[index].v;
        }
    }

    auto * cubic_F = new Function[original_N]();
    auto * linear_F = new Function[original_N]();

    for (int n = 0; n < original_N; n++)
    {
        cubic_F[n].x = original_F[n].x;
        cubic_F[n].y = original_F[n].y;

        linear_F[n].x = original_F[n].x;
        linear_F[n].y = original_F[n].y;
    }

    perform_cubic_interpolation(rescaled_F, cubic_F, skip_x, skip_y, original_Nx, original_Ny);
    perform_linear_interpolation(rescaled_F, linear_F, skip_x, skip_y, original_Nx, original_Ny);

    export_2d("outputs/cubic2D.bin", cubic_F, original_N);
    export_2d("outputs/linear2D.bin", linear_F, original_N);
    export_2d("outputs/original2D.bin", original_F, original_N);
    export_2d("outputs/rescaled2D.bin", rescaled_F, rescaled_N);

    return 0;
}