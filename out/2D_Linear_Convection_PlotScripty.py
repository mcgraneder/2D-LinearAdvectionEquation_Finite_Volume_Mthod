
from mpl_toolkits.mplot3d import Axes3D    ##New Library required for projected 3d plots

import numpy
from matplotlib import pyplot, cm

import numpy as np
import matplotlib.pyplot as plt
'''
#vars
nx = 80
ny = 80

#initialise arrays for 3D plot
x = numpy.linspace(0, 2, nx)
y = numpy.linspace(0, 2, ny)
X, Y = numpy.meshgrid(x, y)   


#dinitialise empty list to append data files to
filelist=["Square_Initial_conditions.txt"]

#data = numpy.loadtxt("Linear_Convection2D_final_data.txt")
#print(data)

#data1 = numpy.array(data).reshape(49,49)

#x = np.linspace(0, 1, 81)
#x = np.linspace(0, 1, 202)
# j = 1
# for i in range(1,201):
#     filelist.append("Linear_Convection_2D_final_data_{}.txt".format(i))

for fname in filelist:
	fig = pyplot.figure(figsize=(10, 9), dpi=100)
	data=np.loadtxt(fname)
	data1 = numpy.array(data).reshape(80,80)
	ax = fig.gca(projection='3d')                                              
	surf = ax.plot_surface(X, Y, data1, cmap=cm.viridis)
	ax.set_zlim(1,2)
	ax.set_ylim(0,2)
	ax.set_xlim(0,2)
	plt.savefig('Linear_Convection_2D_{}.png')
	plt.clf()
	plt.close()
	#j += 1

plt.legend()
plt.show()


'''
from mpl_toolkits.mplot3d import Axes3D    ##New Library required for projected 3d plots

import numpy
from matplotlib import pyplot, cm
#%matplotlib inline
data = numpy.loadtxt("solution_250.txt")
#print(data)
nx = 100
ny = 100
data1 = numpy.array(data).reshape(nx,ny)

print(data1)
###variable declarations

nt = 100
c = 1
dx = 2 / (nx - 1)
dy = 2 / (ny - 1)
sigma = .2
dt = sigma * dx

x = numpy.linspace(-1, 1, nx)
y = numpy.linspace(-1, 1, ny)

u = numpy.ones((ny, nx)) ##create a 1xn vector of 1's
un = numpy.ones((ny, nx)) ##

###Assign initial condition

##set hat function I.C. : u(.5<=x<=1 && .5<=y<=1 ) is 2


fig = pyplot.figure(figsize=(7, 7), dpi=100)
#ax = fig.gca(projection='3d')                      
#X, Y = numpy.meshgrid(x, y)                            
#surf = ax.plot_surface(X, Y, data1, cmap=cm.jet)
plt.contourf(x, y, data1, 90,cmap = plt.cm.jet)
#ax.set_zlim(0,2)
#ax.set_ylim(0,10)
#ax.set_xlim(0,10)
pyplot.show()
'''
filelist = []

nx = 100
ny = 100

x = numpy.linspace(-1, 1, nx)
y = numpy.linspace(-1, 1, ny)
X, Y = numpy.meshgrid(x, y)   

j = 1
for i in range(1,112):
     filelist.append("Solution_{}.txt".format(i))

for fname in filelist:
	fig = plt.figure(figsize=(6, 6), dpi=100)
	data=np.loadtxt(fname)
	data1 = numpy.array(data).reshape(nx,ny)

	#ax = fig.gca(projection='3d')                                              
	#surf = ax.plot_surface(X, Y, data1, cmap=cm.viridis)

	plt.contourf(x, y, data1, 90,cmap = plt.cm.jet)
	plt.xticks([-1, -0.5, 0, 0.5, 1])
	plt.yticks([-1, -0.5, 0, 0.5, 1])
	plt.title("2D Advection SO BW Contour Plot Limited ", fontsize=15)
	#ax.set_xticks([-1, -0.5, 0, 0.5, 1])
	#ax.set_yticks([-1, -0.5, 0, 0.5, 1])
	#ax.set_title("2D Advection SO B&W 3D Plot ", fontsize=13)
	#ax.set_zlim(0,0.18)
	#ax.set_ylim(0,2)
	#ax.set_xlim(0,2)
	plt.savefig('./Animations/BEAM_AND_WARMING_LIMITED_CONTOUR/2D_Advection_S0_BW_LIMITED_CONTOUR_{}.png'.format(j))
	plt.clf()
	plt.close()
	j += 1

plt.legend()
plt.show()

'''
'''
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

                switch (SLOPE_LIMITER)
                {
                    case 1:

                        for (int i = N_GHOST - 1; i < NX + N_GHOST + 1; i++)
                        {
                            for (int j = N_GHOST - 1; j < NY + N_GHOST + 1; j++)
                            {



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

                }
                break;

                case 3:

                    switch (SLOPE_LIMITER)
                    {
                        case 1:

                            for (int i = N_GHOST - 1; i < NX + N_GHOST + 1; i++)
                            {
                                for (int j = N_GHOST - 1; j < NY + N_GHOST + 1; j++)
                                {
                                    double r_x = (u_array[i][j][rkStep] - u_array[i - 1][j][rkStep] + 1e-5) / (u_array[i + 1][j][rkStep] - u_array[i][j][rkStep] + 1e+5);
                                    r_x = max(0, r_x);
                                    double phi_x = (r_x * r_x) / (r_x * r_x + 1.0);

                                    double r_y = (u_array[i][j][rkStep] - u_array[i][j - 1][rkStep] + 1e-5) / (u_array[i][j + 1][rkStep] - u_array[i][j][rkStep] + 1e+5);
                                    r_y = max(0, r_y);
                                    double phi_y = (r_y * r_y) / (r_y * r_y + 1.0);

                                    double du_dx = (u_array[i + 1][j][rkStep] - u_array[i][j][rkStep]) / dx_array[i][j] + 10;
                                    double du_dy = (u_array[i][j + 1][rkStep] - u_array[i][j][rkStep]) / dy_array[i][j];

                                    uWest[i][j] = u_array[i][j][rkStep] - phi_x * du_dx * dx_array[i][j] / 2.0;
                                    uEast[i][j] = u_array[i][j][rkStep] + phi_x * du_dx * dx_array[i][j] / 2.0;

                                    uNorth[i][j] = u_array[i][j][rkStep] - phi_y * du_dy * dy_array[i][j] / 2.0;
                                    uSouth[i][j] = u_array[i][j][rkStep] + phi_y * du_dy * dy_array[i][j] / 2.0;
                                }
                            }
                            break;

                        case 2:
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

                    }
                    break;

                case 4:

                switch (SLOPE_LIMITER)
                {
                    case 1:

                        for (int i = N_GHOST - 1; i < NX + N_GHOST + 1; i++)
                        {
                            for (int j = N_GHOST - 1; j < NY + N_GHOST + 1; j++)
                            {



                            }
                        }
                        break;

                    case 2:
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

                }

                break;


         }


    return;

}
'''