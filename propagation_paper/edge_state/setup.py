import math
import sys
import os
import meep as mp
import matplotlib.pyplot as plt
from meep import mpb

a0 = 1
fcen = 0.46
fwidth = 0.01
sx = 24
sy = 20
Nx = int((sx-4)/(2*a0))-2
Ny = round(sy/(2*3**0.5*a0))
runTime = int(100)
plmWidth = 1

eps_si_5200 = 3.422 ** 2
eps_amor_5200 = 3.1978 ** 2
eps_crys_5200 = 4.63 ** 2

eps_si_1515 = 3.479 ** 2
eps_amor_1515 = 3.334 ** 2
eps_crys_1515 = 5.105 ** 2

r_nontri = 0.12
radius_percentage_nontri = 1.08

r_tri = 0.12
radius_percentage_tri = 0.95

cell = mp.Vector3(sx,sy,0)
geometry = []
for iterx in range(2*Nx):
    for itery in range(Ny):
        #trivial
        geometry += [mp.Cylinder(radius = r_tri, center = mp.Vector3((-Nx+iterx+7/2)*a0+1/3*radius_percentage_tri, (-Ny+itery+1/2)*3**0.5*a0+0*radius_percentage_tri),        material=mp.Medium(epsilon=eps_si_5200)),
                     mp.Cylinder(radius = r_tri, center = mp.Vector3((-Nx+iterx+7/2)*a0+1/6*radius_percentage_tri, (-Ny+itery+1/2)*3**0.5*a0+3**0.5/6*radius_percentage_tri), material=mp.Medium(epsilon=eps_si_5200)),
                     mp.Cylinder(radius = r_tri, center = mp.Vector3((-Nx+iterx+7/2)*a0-1/6*radius_percentage_tri, (-Ny+itery+1/2)*3**0.5*a0+3**0.5/6*radius_percentage_tri), material=mp.Medium(epsilon=eps_si_5200)),
                     mp.Cylinder(radius = r_tri, center = mp.Vector3((-Nx+iterx+7/2)*a0-1/3*radius_percentage_tri, (-Ny+itery+1/2)*3**0.5*a0+0*radius_percentage_tri),        material=mp.Medium(epsilon=eps_si_5200)),
                     mp.Cylinder(radius = r_tri, center = mp.Vector3((-Nx+iterx+7/2)*a0-1/6*radius_percentage_tri, (-Ny+itery+1/2)*3**0.5*a0-3**0.5/6*radius_percentage_tri), material=mp.Medium(epsilon=eps_si_5200)),
                     mp.Cylinder(radius = r_tri, center = mp.Vector3((-Nx+iterx+7/2)*a0+1/6*radius_percentage_tri, (-Ny+itery+1/2)*3**0.5*a0-3**0.5/6*radius_percentage_tri), material=mp.Medium(epsilon=eps_si_5200))]
        geometry += [mp.Cylinder(radius = r_tri, center = mp.Vector3((-Nx+iterx+4)*a0+1/3*radius_percentage_tri, (-Ny+itery+1)*3**0.5*a0+0*radius_percentage_tri),        material=mp.Medium(epsilon=eps_si_5200)),
                     mp.Cylinder(radius = r_tri, center = mp.Vector3((-Nx+iterx+4)*a0+1/6*radius_percentage_tri, (-Ny+itery+1)*3**0.5*a0+3**0.5/6*radius_percentage_tri), material=mp.Medium(epsilon=eps_si_5200)),
                     mp.Cylinder(radius = r_tri, center = mp.Vector3((-Nx+iterx+4)*a0-1/6*radius_percentage_tri, (-Ny+itery+1)*3**0.5*a0+3**0.5/6*radius_percentage_tri), material=mp.Medium(epsilon=eps_si_5200)),
                     mp.Cylinder(radius = r_tri, center = mp.Vector3((-Nx+iterx+4)*a0-1/3*radius_percentage_tri, (-Ny+itery+1)*3**0.5*a0+0*radius_percentage_tri),        material=mp.Medium(epsilon=eps_si_5200)),
                     mp.Cylinder(radius = r_tri, center = mp.Vector3((-Nx+iterx+4)*a0-1/6*radius_percentage_tri, (-Ny+itery+1)*3**0.5*a0-3**0.5/6*radius_percentage_tri), material=mp.Medium(epsilon=eps_si_5200)),
                     mp.Cylinder(radius = r_tri, center = mp.Vector3((-Nx+iterx+4)*a0+1/6*radius_percentage_tri, (-Ny+itery+1)*3**0.5*a0-3**0.5/6*radius_percentage_tri), material=mp.Medium(epsilon=eps_si_5200))]

        #non-trivial
        geometry += [mp.Cylinder(radius = r_nontri, center = mp.Vector3((-Nx+iterx+7/2)*a0+1/3*radius_percentage_nontri, (itery+1/2)*3**0.5*a0+0*radius_percentage_nontri),        material=mp.Medium(epsilon=eps_si_5200)),
                     mp.Cylinder(radius = r_nontri, center = mp.Vector3((-Nx+iterx+7/2)*a0+1/6*radius_percentage_nontri, (itery+1/2)*3**0.5*a0+3**0.5/6*radius_percentage_nontri), material=mp.Medium(epsilon=eps_si_5200)),
                     mp.Cylinder(radius = r_nontri, center = mp.Vector3((-Nx+iterx+7/2)*a0-1/6*radius_percentage_nontri, (itery+1/2)*3**0.5*a0+3**0.5/6*radius_percentage_nontri), material=mp.Medium(epsilon=eps_si_5200)),
                     mp.Cylinder(radius = r_nontri, center = mp.Vector3((-Nx+iterx+7/2)*a0-1/3*radius_percentage_nontri, (itery+1/2)*3**0.5*a0+0*radius_percentage_nontri),        material=mp.Medium(epsilon=eps_si_5200)),
                     mp.Cylinder(radius = r_nontri, center = mp.Vector3((-Nx+iterx+7/2)*a0-1/6*radius_percentage_nontri, (itery+1/2)*3**0.5*a0-3**0.5/6*radius_percentage_nontri), material=mp.Medium(epsilon=eps_si_5200)),
                     mp.Cylinder(radius = r_nontri, center = mp.Vector3((-Nx+iterx+7/2)*a0+1/6*radius_percentage_nontri, (itery+1/2)*3**0.5*a0-3**0.5/6*radius_percentage_nontri), material=mp.Medium(epsilon=eps_si_5200))]
        geometry += [mp.Cylinder(radius = r_nontri, center = mp.Vector3((-Nx+iterx+4)*a0+1/3*radius_percentage_nontri, (itery+1)*3**0.5*a0+0*radius_percentage_nontri),        material=mp.Medium(epsilon=eps_si_5200)),
                     mp.Cylinder(radius = r_nontri, center = mp.Vector3((-Nx+iterx+4)*a0+1/6*radius_percentage_nontri, (itery+1)*3**0.5*a0+3**0.5/6*radius_percentage_nontri), material=mp.Medium(epsilon=eps_si_5200)),
                     mp.Cylinder(radius = r_nontri, center = mp.Vector3((-Nx+iterx+4)*a0-1/6*radius_percentage_nontri, (itery+1)*3**0.5*a0+3**0.5/6*radius_percentage_nontri), material=mp.Medium(epsilon=eps_si_5200)),
                     mp.Cylinder(radius = r_nontri, center = mp.Vector3((-Nx+iterx+4)*a0-1/3*radius_percentage_nontri, (itery+1)*3**0.5*a0+0*radius_percentage_nontri),        material=mp.Medium(epsilon=eps_si_5200)),
                     mp.Cylinder(radius = r_nontri, center = mp.Vector3((-Nx+iterx+4)*a0-1/6*radius_percentage_nontri, (itery+1)*3**0.5*a0-3**0.5/6*radius_percentage_nontri), material=mp.Medium(epsilon=eps_si_5200)),
                     mp.Cylinder(radius = r_nontri, center = mp.Vector3((-Nx+iterx+4)*a0+1/6*radius_percentage_nontri, (itery+1)*3**0.5*a0-3**0.5/6*radius_percentage_nontri), material=mp.Medium(epsilon=eps_si_5200))]


sources = [mp.Source(mp.GaussianSource(frequency=fcen, fwidth = fwidth),
                     mp.Hz,
                     mp.Vector3(-10.5,0))
                     ]

pml_layers = [mp.PML(plmWidth)]

resolution = 32

sim = mp.Simulation(cell_size=cell,
                    boundary_layers=pml_layers,
                    geometry=geometry,
                    default_material = mp.Medium(index=1),
                    sources=sources,
                    resolution=resolution)

sim.use_output_directory("wave_in_band_gap")

def myhello(thisSim):
  # print(-0.4*sx+sim.meep_time()*v)
    if(thisSim.meep_time() > 190):
        movingSources = []
        thisSim.change_sources(movingSources)
    return 

sim.run(
        mp.at_beginning(mp.output_epsilon),
        mp.to_appended("ez",mp.at_every(1,mp.output_efield_z)),
        myhello,
        until=runTime)

# eps_data = sim.get_array(center=mp.Vector3(), size=cell, component=mp.Dielectric)
# plt.figure(dpi=100)
# plt.imshow(eps_data.transpose(), interpolation='spline36', cmap='binary')
# plt.axis('off')
# plt.savefig('./wave_in_band_gap/eps.png')

