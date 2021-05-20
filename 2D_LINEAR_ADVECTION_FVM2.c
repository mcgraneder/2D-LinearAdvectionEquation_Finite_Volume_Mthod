#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <stdbool.h>
#define PI 3.14
#define max(a,b)             \
({                           \
    __typeof__ (a) _a = (a); \
    __typeof__ (b) _b = (b); \
    _a > _b ? _a : _b;       \
})


void InitialConditions(double*** u_array, double** x_centroid_array, double** y_centroid_array, double** dx_array, double** dy_array, double ARRAY_SIZE_X, int rkStep, int INITIAL_CONDITIONS, FILE* IC_fp, FILE* IC2_fp);
double getTimeStep(double** dx_array, double** dy_array, double INITIAL_VELOCITY_X, double INITIAL_VELOCITY_Y, double COURANT_NUM);
void updateGhostCells(double*** u_array, int ARRAY_SIZE_X, int NX, int NY, int N_GHOST, int rkStep);
void reconstructVariables(double*** u_array, double** uEast, double** uWest, double** uNorth, double** uSouth, double** dx_array, double** dy_array, int NUMERICAL_METHOD, int N_GHOST, int NX, int NY, int rkStep, int SLOPE_LIMITER);
void fluxCalculator(double** uEast, double** uWest, double** uNorth, double** uSouth, double** total_flux, double** dx_array, double** dy_array, int NX, int NY, int N_GHOST, double INITIAL_VELOCITY, int ARRAY_SIZE_X);
void TimeIntegration(int N_RUNGE_KUTTA_STEPS, double*** u_array, double** dx_array, double** dy_array, double** total_flux, double dt, int rkStep, int NX, int NY, int N_GHOST);
void copySolution(double*** u_array, int NX, int NY, int N_GHOST, int N_RUNGE_KUTTA_STEPS, int rkStep);
// void errorAnalysis(double** u_array, double** u_Exact, int NX, int NX2, int N_GHOST, int error, int error_mesh1, int error_mesh2, int rkStep);
void writeSolutionToFile(double*** u_array, int N_GHOST, int NX, int NY);

int main()
{
    //open files to hold initial condition data
    FILE* IC_fp = NULL;
    IC_fp = fopen("./out/Square_Initial_conditions.txt", "w");

    FILE* IC2_fp = NULL;
    IC2_fp = fopen("./out/Sine_Initial_Conditions.txt", "w");

    //error handle
    if (IC_fp == NULL || IC2_fp == NULL)
    {
        printf("\n");
        printf("FILE OPEN FAILED");
        printf("\n");
    }

    int NX = 100;
    int NY = 100;

    int N_GHOST = 4;
    int N_RUNGE_KUTTA_STEPS = 3;
    int INITIAL_CONDITIONS = 1;
    int NUMERICAL_METHOD = 5;
    int SLOPE_LIMITER = 1;
    double STOPPING_TIME = 0.2;
    int MAX_TIME_STEPS = 1000000;
    double MAX_X = 1.0;
    double MIN_X = -1.0;
    double MAX_Y = 1.0;
    double MIN_Y = -1.0;
    double INITIAL_VELOCITY_X = 1;
    double INITIAL_VELOCITY_Y = 1;
    double COURANT_NUM = 0.5;

    double DX = (MAX_X - MIN_X) / NX;
    double DY = (MAX_Y - MIN_Y) / NY;

    int ARRAY_SIZE_X = NX + 2 * N_GHOST;
    int ARRAY_SIZE_Y = NY + 2 * N_GHOST;

    //initialise arrrays we need two arrays since we are using 2nd order
    //runge kutta time integration
    double*** u_array = (double***)calloc(ARRAY_SIZE_X, sizeof(double**));
    double*** u_Exact = (double***)calloc(ARRAY_SIZE_X, sizeof(double**));

    double** uEast = (double**)calloc(NX + 2 * N_GHOST, sizeof(double*));
    double** uWest = (double**)calloc(NX + 2 * N_GHOST, sizeof(double*));
    double** uNorth = (double**)calloc(NX + 2 * N_GHOST, sizeof(double*));
    double** uSouth = (double**)calloc(NX + 2 * N_GHOST, sizeof(double*));

    double** dx_array = (double**)calloc(NX + 2 * N_GHOST, sizeof(double*));
    double** dy_array = (double**)calloc(NY + 2 * N_GHOST, sizeof(double*));
    double** x_centroid_array = (double**)calloc(NX + 2 * N_GHOST, sizeof(double*));
    double** y_centroid_array = (double**)calloc(NY + 2 * N_GHOST, sizeof(double*));

    double** total_flux = (double**)calloc(NX + 2 * N_GHOST, sizeof(double*));

    for (int i = 0; i < ARRAY_SIZE_X; i++)
    {
        u_array[i] = (double**)calloc(ARRAY_SIZE_X, sizeof(double*));
        u_Exact[i] = (double**)calloc(ARRAY_SIZE_X, sizeof(double*));

        dx_array[i] = (double*)calloc(NX + 2 * N_GHOST, sizeof(double));
        dy_array[i] = (double*)calloc(NY + 2 * N_GHOST, sizeof(double));

        x_centroid_array[i] = (double*)calloc(NX + 2 * N_GHOST, sizeof(double));
        y_centroid_array[i] = (double*)calloc(NY + 2 * N_GHOST, sizeof(double));

        uEast[i] = (double*)calloc(NX + 2 * N_GHOST, sizeof(double));
        uWest[i] = (double*)calloc(NX + 2 * N_GHOST, sizeof(double));
        uNorth[i] = (double*)calloc(NX + 2 * N_GHOST, sizeof(double));
        uSouth[i] = (double*)calloc(NX + 2 * N_GHOST, sizeof(double));

        total_flux[i] = (double*)calloc(NX + 2 * N_GHOST, sizeof(double));

        for (int j = 0; j < ARRAY_SIZE_Y; j++)
        {
            u_array[i][j] = (double*)calloc(N_RUNGE_KUTTA_STEPS + 1, sizeof(double));
            u_Exact[i][j] = (double*)calloc(N_RUNGE_KUTTA_STEPS + 1, sizeof(double));
        }
    }

    //populate DX ann cx
    for (int i = 0; i < NX + 2 * N_GHOST; i++ )
    {
        for (int j = 0; j < ARRAY_SIZE_Y; j++)
        {
            dx_array[i][j] = DX;
            dy_array[i][j] = DY;

            x_centroid_array[i][j] = MIN_X + DX * (i + 0.5 - N_GHOST);
            y_centroid_array[j][i] = MIN_Y + DY * (j + 0.5 - N_GHOST);
            //printf(" %14f", x_centroid_array[i][j]);
        }

    }

    int rkStep = 0;
    //initialise ICs wec can choose between square wave ICs or sine wave ICs
    //InitialConditions(u_array, u_Exact, centroid_array, dx_array, ARRAY_SIZE, rkStep, INITIAL_CONDITIONS, IC_fp, IC2_fp);
    InitialConditions(u_array, x_centroid_array, y_centroid_array, dx_array, dy_array, ARRAY_SIZE_X, rkStep, INITIAL_CONDITIONS, IC_fp, IC2_fp);


     double time = 0;
     bool lastTimeStep = false;
     for (int t = 0; t < MAX_TIME_STEPS; t++)
     {
         double dt = getTimeStep(dx_array, dy_array, INITIAL_VELOCITY_X, INITIAL_VELOCITY_Y, COURANT_NUM);
         //printf("%f\n", dt);
    //     //printf("Timestep %d:", t + 1);
    //     //printf("\n");

         if (time + dt > STOPPING_TIME)
         {
             dt = STOPPING_TIME - time;
             lastTimeStep = true;
         }

         for (rkStep = 0; rkStep < N_RUNGE_KUTTA_STEPS; rkStep++)
         {

             updateGhostCells(u_array, ARRAY_SIZE_X, NX, NY, N_GHOST, rkStep);

             reconstructVariables(u_array, uEast, uWest, uNorth, uSouth, dx_array, dy_array, NUMERICAL_METHOD, N_GHOST, NX, NY, rkStep, SLOPE_LIMITER);

             fluxCalculator(uEast, uWest, uNorth, uSouth, total_flux, dx_array, dy_array, NX, NY, N_GHOST, INITIAL_VELOCITY_X, ARRAY_SIZE_X);

             TimeIntegration(N_RUNGE_KUTTA_STEPS, u_array, dx_array, dy_array, total_flux, dt, rkStep, NX, NY, N_GHOST);

             copySolution(u_array, NX, NY, N_GHOST, N_RUNGE_KUTTA_STEPS, rkStep);

         }

         //print each timestep to its own file
         if (t % 5 == 0)
         {
             //writeSolutionToFile(u_array, N_GHOST, NX, NY);
         }


         time += dt;
         if (lastTimeStep)
         {
             break;
         }

     }

    //print final tmestep to screen
    printf("\n");
    for (int i = N_GHOST; i < NX + N_GHOST; i++)
    {
        for (int j = N_GHOST; j < NY + N_GHOST; j++)
        {

            //double dt;
            //u_array[i][rkStep + 1] = u_array[i][rkStep] + dt / dx_array[i] * total_flux[i];
            printf(" %14f", u_array[i][j][N_RUNGE_KUTTA_STEPS]);
            fprintf(IC2_fp, "%f\n", u_array[i][j][rkStep] );
        }

    }


    return EXIT_SUCCESS;
}




void InitialConditions(double*** u_array, double** x_centroid_array, double** y_centroid_array, double** dx_array, double** dy_array, double ARRAY_SIZE_X, int rkStep, int INITIAL_CONDITIONS, FILE* IC_fp, FILE* IC2_fp)
{
    switch (INITIAL_CONDITIONS)
    {
        case 1:
            for (int i = 0; i < ARRAY_SIZE_X; i++)
            {
                for (int j = 0; j < ARRAY_SIZE_X; j++)
                {

                    //u_Exact[i][j][rkStep] = 0.0;


                    //comment our sine wave ICs for the moment
                    //u_array[i] = sin(PI * cx[i]);
                    //
                    //u_array[rkStep] = sin(PI * centroid_array[i]);
                    if (x_centroid_array[i][j] > -0.25 && x_centroid_array[i][j] < 0.25 && y_centroid_array[j][i] > -0.25 && y_centroid_array[j][i] < 0.25)
                    {
                        u_array[i][j][rkStep] = 1.0;
                        //u_Exact[i][j][rkStep] = 1.0;
                    }
                    else
                    {
                        u_array[i][j][rkStep] = 0.0;
                    }
                    //printf(" %14f", u_array[i][j][rkStep]);
                    fprintf(IC_fp, "%f\n", u_array[i][j][rkStep]);
                }

            }
            break;

         case 2:
             for (int i = 0; i < ARRAY_SIZE_X; i++)
             {
                 for (int j = 0; j < ARRAY_SIZE_X; j++)
                 {

                     u_array[i][j][rkStep] = sin(PI * x_centroid_array[i][j]) + 1;
                     //u_Exact[i][j][rkStep] = sin(PI * x_centroid_array[i][j]);
                     //printf("%f\n", u_array[i][rkStep]);
                     //fprintf(IC2_fp, "%f\n", u_array[i][rkStep]);
                 }
             }
             break;

        // case 3:
        //     for (int i = 0; i < ARRAY_SIZE; i++)
        //     {
        //         if (centroid_array[i] <= 0.15)
        //         {
        //             u_array[i][rkStep] = 1.0;
        //             u_Exact[i][rkStep] = 1.0;
        //         }
        //         else
        //         {
        //             u_array[i][rkStep] = 0;
        //             u_Exact[i][rkStep] = 0;
        //         }
        //         //printf("%f\n", u_array[i][rkStep]);
        //     }
        //     break;

    }

}

double getTimeStep(double** dx_array, double** dy_array, double INITIAL_VELOCITY_X, double INITIAL_VELOCITY_Y, double COURANT_NUM)
{
    double dt1 = dx_array[0][0] / INITIAL_VELOCITY_X;
    double dt2 = dy_array[0][0] / INITIAL_VELOCITY_Y;

    double dt = 1.0 / (1.0 / dt1 + 1.0 / dt2);
    return dt * COURANT_NUM;
}



void updateGhostCells(double*** u_array, int ARRAY_SIZE_X, int NX, int NY, int N_GHOST, int rkStep)
{
    for (int j = 0; j < ARRAY_SIZE_X; j++)
    {

        for (int ghostcell = 0; ghostcell < N_GHOST; ghostcell++)
        {

            //we want to apply periodic boundary conditions. Whereby at the end of wach
            //timestep we copy the last two cells to the first ghost cells, and the first two ghost cells
            //with the last two ghost cells
            u_array[ghostcell][j][rkStep] = u_array[NX + ghostcell][j][rkStep];
            u_array[NX + N_GHOST + ghostcell][j][rkStep] = u_array[N_GHOST + ghostcell][j][rkStep];
            //printf("%f\n", u_array[ghostcell][rkStep]);
        }
    }

    for (int i = 0; i < ARRAY_SIZE_X; i++)
    {

        for (int ghostcell = 0; ghostcell < N_GHOST; ghostcell++)
        {

            //we want to apply periodic boundary conditions. Whereby at the end of wach
            //timestep we copy the last two cells to the first ghost cells, and the first two ghost cells
            //with the last two ghost cells
            u_array[i][ghostcell][rkStep] = u_array[i][NY + ghostcell][rkStep];
            u_array[i][NY + N_GHOST + ghostcell][rkStep] = u_array[i][N_GHOST + ghostcell][rkStep];
                //printf("%f\n", u_array[ghostcell][rkStep]);
        }
    }

    return;
}

void reconstructVariables(double*** u_array, double** uEast, double** uWest, double** uNorth, double** uSouth, double** dx_array, double** dy_array, int NUMERICAL_METHOD, int N_GHOST, int NX, int NY, int rkStep, int SLOPE_LIMITER)
{
    //int SLOPE_LIMITERS = 2; Once we have applied the BCs we want to reconstruct our variables
        switch(NUMERICAL_METHOD)
        {
            case 1:
                //we want to index from one cell before the internal cells to one cell after the size of the mesh
                //so that we can account for the fluxes at i - 1/2 and i + 1/2 for the first and last cells
                for (int i = N_GHOST - 1; i < NX + N_GHOST + 1; i++)
                {
                    for (int j = N_GHOST -1; j < NY +N_GHOST + 1; j++)
                    {
                        //first order reconstruction we copy the centroid values from u_array into our i+1/2 and i-1/2 arrays
                        //which are represented by uEast and uWest. uEast and uWest reprresents the interface values of u_array
                        uEast[i][j] = u_array[i][j][rkStep];
                        uWest[i][j] = u_array[i][j][rkStep];
                        uNorth[i][j] = u_array[i][j][rkStep];
                        uSouth[i][j] = u_array[i][j][rkStep];
                    }
                }
                break;

            case 2:


                for (int i = N_GHOST - 1; i < NX + N_GHOST + 1; i++)
                {
                    for (int j = N_GHOST - 1; j < NY + N_GHOST + 1; j++)
                    {

                        double du_dx = (u_array[i][j][rkStep] - u_array[i - 1][j][rkStep]) / dx_array[i][j];
                        double du_dy = (u_array[i][j][rkStep] - u_array[i][j - 1][rkStep]) / dy_array[i][j];

                        uWest[i][j] = u_array[i][j][rkStep] - du_dx * dx_array[i][j] / 2.0;
                        uEast[i][j] = u_array[i][j][rkStep] + du_dx * dx_array[i][j] / 2.0;

                        uNorth[i][j] = u_array[i][j][rkStep] - du_dy * dy_array[i][j] / 2.0;
                        uSouth[i][j] = u_array[i][j][rkStep] + du_dy * dy_array[i][j] / 2.0;
                    }
                }
                break;




            case 3:


                for (int i = N_GHOST - 1; i < NX + N_GHOST + 1; i++)
                {
                    for (int j = N_GHOST - 1; j < NY + N_GHOST + 1; j++)
                    {

                        double du_dx = (u_array[i + 1][j][rkStep] - u_array[i][j][rkStep]) / dx_array[i][j];
                        double du_dy = (u_array[i][j + 1][rkStep] - u_array[i][j][rkStep]) / dy_array[i][j];

                        uWest[i][j] = u_array[i][j][rkStep] - du_dx * dx_array[i][j] / 2.0;
                        uEast[i][j] = u_array[i][j][rkStep] + du_dx * dx_array[i][j] / 2.0;

                        uNorth[i][j] = u_array[i][j][rkStep] - du_dy * dy_array[i][j] / 2.0;
                        uSouth[i][j] = u_array[i][j][rkStep] + du_dy * dy_array[i][j] / 2.0;
                    }
                }
                break;



            case 4:


                for (int i = N_GHOST - 1; i < NX + N_GHOST + 1; i++)
                {
                    for (int j = N_GHOST - 1; j < NY + N_GHOST + 1; j++)
                    {

                        double du_dx = (u_array[i + 1][j][rkStep] - u_array[i - 1][j][rkStep]) / 2.0 / dx_array[i][j];
                        double du_dy = (u_array[i][j + 1][rkStep] - u_array[i][j - 1][rkStep]) / 2.0 / dy_array[i][j];

                        uWest[i][j] = u_array[i][j][rkStep] - du_dx * dx_array[i][j] / 2.0;
                        uEast[i][j] = u_array[i][j][rkStep] + du_dx * dx_array[i][j] / 2.0;

                        uNorth[i][j] = u_array[i][j][rkStep] - du_dy * dy_array[i][j] / 2.0;
                        uSouth[i][j] = u_array[i][j][rkStep] + du_dy * dy_array[i][j] / 2.0;
                    }
                }
                break;

            case 5:

                for (int i = N_GHOST - 1; i < NX + N_GHOST + 1; i++)
                {
                    for (int j = N_GHOST - 1; j < NY + N_GHOST + 1; j++)
                    {
                        double numerator = u_array[i][j] - u_array[i - 1][j];
                        double denominator = u_array[i + 1][j] - u_array[i][j];

                        double r_x = numerator / denominator;
                        if (fabs(denominator) < 1e-20)
                        {
                            r_x = 0.0;
                        }
                        r_x = fmax(0, r_x);
                        double phi_x = (r_x * r_x + r_x) / (r_x * r_x + 1.0);

                        double numerator1 = u_array[i][j] - u_array[i][j - 1];
                        double denominator1 = u_array[i][j + 1] - u_array[i][j];

                        double r_y = numerator1 / denominator1;
                        if (fabs(denominator1) < 1e-20)
                        {
                            r_y = 0.0;
                        }
                        r_y = fmax(0, r_y);
                        double phi_y = (r_y * r_y + r_y) / (r_y * r_y + 1.0);

                        double du_dx = (u_array[i + 1][j][rkStep] - u_array[i][j][rkStep]) / dx_array[i][j];
                        double du_dy = (u_array[i][j + 1][rkStep] - u_array[i][j][rkStep]) / dy_array[i][j];

                        uWest[i][j] = u_array[i][j][rkStep] - phi_x * du_dx * dx_array[i][j] / 2.0;
                        uEast[i][j] = u_array[i][j][rkStep] + phi_x * du_dx * dx_array[i][j] / 2.0;

                        uSouth[i][j] = u_array[i][j][rkStep] - phi_y * du_dy * dy_array[i][j] / 2.0;
                        uNorth[i][j] = u_array[i][j][rkStep] + phi_y * du_dy * dy_array[i][j] / 2.0;
                    }
                }
                break;

            case 6:

                for (int i = N_GHOST - 1; i < NX + N_GHOST + 1; i++)
                {
                    for (int j = N_GHOST - 1; j < NY + N_GHOST + 1; j++)
                    {
                        double numerator = u_array[i + 1][j] - u_array[i][j];
                        double denominator = u_array[i][j] - u_array[i - 1][j];

                        double r_x = numerator / denominator;
                        if (fabs(denominator) < 1e-100)
                        {
                            r_x = 0.0;
                        }
                        r_x = fmax(0, r_x);
                        double phi_x = (r_x * r_x + r_x) / (r_x * r_x + 1.0);

                        double numerator1 = u_array[i][j + 1] - u_array[i][j];
                        double denominator1 = u_array[i][j] - u_array[i][j - 1];

                        double r_y = numerator1 / denominator1;
                        if (fabs(denominator1) < 1e-100)
                        {
                            r_y = 0.0;
                        }
                        r_y = fmax(0, r_y);
                        double phi_y = (r_y * r_y + r_y) / (r_y * r_y + 1.0);

                        double du_dx = (u_array[i][j][rkStep] - u_array[i - 1][j][rkStep]) / dx_array[i][j];
                        double du_dy = (u_array[i][j][rkStep] - u_array[i][j - 1][rkStep]) / dy_array[i][j];

                        uWest[i][j] = u_array[i][j][rkStep] - phi_x * du_dx * dx_array[i][j] / 2.0;
                        uEast[i][j] = u_array[i][j][rkStep] + phi_x * du_dx * dx_array[i][j] / 2.0;

                        uSouth[i][j] = u_array[i][j][rkStep] - phi_y * du_dy * dy_array[i][j] / 2.0;
                        uNorth[i][j] = u_array[i][j][rkStep] + phi_y * du_dy * dy_array[i][j] / 2.0;
                    }
                }
                break;

            case 7:

                for (int i = N_GHOST - 1; i < NX + N_GHOST + 1; i++)
                {
                    for (int j = N_GHOST - 1; j < NY + N_GHOST + 1; j++)
                    {
                        double numerator = u_array[i + 1][j] - u_array[i - 1][j];
                        double denominator = u_array[i - 1][j] - u_array[i + 1][j];

                        double r_x = numerator / denominator;
                        if (fabs(denominator) < 1e-100)
                        {
                            r_x = 0.0;
                        }
                        r_x = fmax(0, r_x);
                        double phi_x = (r_x * r_x + r_x) / (r_x * r_x + 1.0);

                        double numerator1 = u_array[i][j + 1] - u_array[i][j - 1];
                        double denominator1 = u_array[i][j - 1] - u_array[i][j + 1];

                        double r_y = numerator1 / denominator1;
                        if (fabs(denominator1) < 1e-100)
                        {
                            r_y = 0.0;
                        }
                        r_y = fmax(0, r_y);
                        double phi_y = (r_y * r_y + r_y) / (r_y * r_y + 1.0);

                        double du_dx = (u_array[i + 1][j][rkStep] - u_array[i - 1][j][rkStep]) / 2.0 / dx_array[i][j];
                        double du_dy = (u_array[i][j + 1][rkStep] - u_array[i][j - 1][rkStep]) / 2.0 / dy_array[i][j];

                        uWest[i][j] = u_array[i][j][rkStep] - phi_x * du_dx * dx_array[i][j] / 2.0;
                        uEast[i][j] = u_array[i][j][rkStep] + phi_x * du_dx * dx_array[i][j] / 2.0;

                        uSouth[i][j] = u_array[i][j][rkStep] - phi_y * du_dy * dy_array[i][j] / 2.0;
                        uNorth[i][j] = u_array[i][j][rkStep] + phi_y * du_dy * dy_array[i][j] / 2.0;
                    }
                }
                break;

            case 8:

                for (int i = N_GHOST - 1; i < NX + N_GHOST + 1; i++)
                {
                    for (int j = N_GHOST - 1; j < NY + N_GHOST + 1; j++)
                    {
                        double EPS = 1e-12;
                        double uij = u_array[i][j][rkStep];
                        double uip1j = u_array[i + 1][j][rkStep];
                        double uim1j = u_array[i - 1][j][rkStep];
                        double uijp1 = u_array[i][j + 1][rkStep];
                        double uijm1 = u_array[i][j -1][rkStep];

                        double dx = dx_array[i][j];
                        double dy = dy_array[i][j];

                        double du1 = uip1j - uij;
                        double du2 = uim1j - uij;
                        double w1 = 1.0 / (du1 * du1 + EPS);
                        double w2 = 1.0 / (du2 * du2 + EPS);

                        double du_dx = (w1 * du1 - w2 * du2) / (w1 + w2) / dx;
                        uWest[i][j] = u_array[i][j][rkStep] - du_dx * dx / 2.0;
                        uEast[i][j] = u_array[i][j][rkStep] - du_dx * dx / 2.0;

                        du1 = uijp1 - uij;
                        du2 = uijm1 - uij;
                        w1 = 1.0 / (du1 * du1 + EPS);
                        w2 = 1.0 / (du2 * du2 + EPS);

                        double du_dy = (w1 * du1 - w2 * du2) / (w1 + w2) / dy;
                        uSouth[i][j] = u_array[i][j][rkStep] - du_dy * dy / 2.0;
                        uNorth[i][j] = u_array[i][j][rkStep] - du_dy * dy / 2.0; 

                    }
                }
                break;



         }


    return;

}

// //

//after we reconstruct our variables and store the interface fluxes in uEast and uWest we need to
//calculate the fluxes at the interfaces
void fluxCalculator(double** uEast, double** uWest, double** uNorth, double** uSouth, double** total_flux, double** dx_array, double** dy_array, int NX, int NY, int N_GHOST, double INITIAL_VELOCITY, int ARRAY_SIZE_X)
{
    for (int i = 0; i < ARRAY_SIZE_X; i++)
    {
        for (int j = 0; j < ARRAY_SIZE_X; j++)
        {
            total_flux[i][j] = 0.0;
        }
    }

    double dx = dx_array[0][0];
    double dy = dy_array[0][0];

    //we want to look through all of the internal cells
    for (int j = N_GHOST; j < NY + N_GHOST; j++)
    {

        for (int verticalEdgeIndex = N_GHOST; verticalEdgeIndex < NY + N_GHOST + 1; verticalEdgeIndex++ )
        {
            //temporary variables to store the indices for the lax fredrichs flux formula
            double leftValue = uEast[verticalEdgeIndex - 1][j];
            double rightValue = uWest[verticalEdgeIndex][j];

            //lax fredrichs schem for calculating the flux. We are using postive advection velocity
            //so we use upwind scheme
            double flux = 0.5 * INITIAL_VELOCITY * (leftValue + rightValue) - 0.5 * fabs(INITIAL_VELOCITY) * (rightValue - leftValue);

            flux *= dy;
            //whatever flux we ge we need to add it to the cell. The total flux array will store our fluxes
            //over the domain

            //left side of the interface subtracted. right side of the flux added
            total_flux[verticalEdgeIndex - 1][j] -= flux;
            total_flux[verticalEdgeIndex][j] += flux;
            //printf("%f\n", total_flux[verticalEdgeIndex]);
        }
    }

    for (int i = N_GHOST; i < NX + N_GHOST; i++)
    {

        for (int horizontalEdgeIndex = N_GHOST; horizontalEdgeIndex < NY + N_GHOST + 1; horizontalEdgeIndex++ )
        {
            //temporary variables to store the indices for the lax fredrichs flux formula
            double leftValue = uNorth[i][horizontalEdgeIndex - 1];
            double rightValue = uSouth[i][horizontalEdgeIndex];

            //lax fredrichs schem for calculating the flux. We are using postive advection velocity
            //so we use upwind scheme
            double flux = 0.5 * INITIAL_VELOCITY * (leftValue + rightValue) - 0.5 * fabs(INITIAL_VELOCITY) * (rightValue - leftValue);

            flux *= dx;
            //whatever flux we ge we need to add it to the cell. The total flux array will store our fluxes
            //over the domain

            //left side of the interface subtracted. right side of the flux added
            total_flux[i][horizontalEdgeIndex - 1] -= flux;
            total_flux[i][horizontalEdgeIndex] += flux;
            //printf("%f\n", total_flux[horizontalEdgeIndex]);
        }
    }


    return;
}

void TimeIntegration(int N_RUNGE_KUTTA_STEPS, double*** u_array, double** dx_array, double** dy_array, double** total_flux, double dt, int rkStep, int NX, int NY, int N_GHOST)
{
    double area = dx_array[0][0] * dy_array[0][0];
    //printf("%f\n", area);
    switch(N_RUNGE_KUTTA_STEPS)
    {
        case 1:
            for (int i = N_GHOST; i < NX + N_GHOST; i++)
            {
                for (int j = N_GHOST; j < NY + N_GHOST; j++)
                {


                    u_array[i][j][rkStep + 1] = u_array[i][j][rkStep] + dt / area * total_flux[i][j];
                    //printf("%f\n", u_array[i][j][rkStep]);
                    //u_array[i][rkStep + 1] = 0.5 * (u_array[i - 1][rkStep] + u_array[i + 1][rkStep]) - 0.5 * (dt / dx_array[i]) * total_flux[i];
                }
            }
            break;

        case 2:
            for (int i = 0; i < NX + N_GHOST; i++)
            {


                for (int j = N_GHOST; j < NY + N_GHOST; j++)
                {
                    switch (rkStep)
                    {
                        case 0:
                            u_array[i][j][rkStep + 1] = u_array[i][j][rkStep] + dt / area * total_flux[i][j];
                            break;

                        case 1:
                            u_array[i][j][rkStep + 1] = 0.5 * (u_array[i][j][rkStep - 1] + u_array[i][j][rkStep] + dt / area * total_flux[i][j]);
                            break;
                        default:
                            printf("Timestep not available");
                    }
                }
            }
            break;

        case 3:
            for (int i = 0; i < NX + N_GHOST; i++)
            {


                for (int j = N_GHOST; j < NX + N_GHOST; j++)
                {
                    switch(rkStep)
                    {
                        case 0:
                            u_array[i][j][rkStep + 1] = u_array[i][j][rkStep] + dt / area * total_flux[i][j];
                            break;

                        case 1:
                            u_array[i][j][rkStep + 1] = 3.0 / 4.0 * u_array[i][j][rkStep - 1] + 1.0 / 4.0 * u_array[i][j][rkStep] + 1.0 / 4.0 * dt / area * total_flux[i][j];
                            break;

                        case 2:
                             u_array[i][j][rkStep + 1] = 1.0 / 3.0 * u_array[i][j][rkStep - 2] + 2.0 / 3.0 * u_array[i][j][rkStep] + 2.0 / 3.0 * dt / area * total_flux[i][j];
                            break;

                    }
                }
            }
        break;
    }

    return;
}



void copySolution(double*** u_array, int NX, int NY, int N_GHOST, int N_RUNGE_KUTTA_STEPS, int rkStep)
{
    for (int i = N_GHOST; i < NX + N_GHOST; i++)
    {
        for (int j = 0; j < NY + N_GHOST; j++)
        {


            u_array[i][j][0] = u_array[i][j][N_RUNGE_KUTTA_STEPS];
            //printf(" %14f", u_array[i][0]);
        }
    }

    return;
}

//function that writes each timestep to a new file
void writeSolutionToFile(double*** u_array, int N_GHOST, int NX, int NY)
 {
    //printf("\nWriting to file... ");
    static int fileIndex = 0;
    char fileName[100];


    sprintf(fileName, "./out/solution_%d.txt", fileIndex);
    FILE* file = fopen(fileName, "w");


    // Write cell data
    const int rkStep = 0;
    int i;

    for (i = N_GHOST; i < NX + N_GHOST; i++)
    {
        for (int j = 0; j < NY + N_GHOST; j++)
        {
            fprintf(file, "%f\n", u_array[i][j][rkStep]);
        }

    }

    fclose(file);

    fileIndex++;
    //printf("done.\n");
}

void ErrorCaclulator(int ARRAY_SIZE_X)
{
    double errorL1;
    double errorL2;
    double errorLinf;

    double** initCells = (double**)calloc(ARRAY_SIZE_X, sizeof(double*));
    for (int i = 0; i < ARRAY_SIZE_X; i++)
    {
        initCells[i] = (double*)calloc(ARRAY_SIZE_X, sizeof(double));
    }

    for (int i = 0; i < ARRAY_SIZE_X; i++)
    {
        for (int j = 0; j < ARRAY_SIZE_X; j++)
        {

        }
    }
}

