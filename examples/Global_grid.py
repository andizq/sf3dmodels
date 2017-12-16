import BuildGlobalGrid as BGG
import Model
import Plot_model as Pm
import Utils as U

sizex = sizey = sizez = 1000 * U.AU
Nx = Ny = Nz = 120
GRID = Model.grid([sizex, sizey, sizez], [Nx, Ny, Nz])
list_sub = ['datatab_Main.dat', 'datatab_Burger.dat']
global_prop = BGG.overlap(GRID, submodels = list_sub)
