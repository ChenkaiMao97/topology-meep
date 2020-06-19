import math
import sys
import os
import meep as mp
import matplotlib.pyplot as plt
from meep import mpb

a0 = 1
fcen = 0.6
fwidth = 0.03
sx = 20
sy = 20
Nx = int(sx/(a0))
Ny = round(sy/(3**0.5*a0))
runTime = int(100)
plmWidth = 1

eps_si_5200 = 3.422 ** 2
eps_si_1515 = 3.479 ** 2

eps_GST326_amor = 12.8
eps_GST326_crys = 40.0

eps_InSb_amor = 24.8
eps_InSb_crys = 15.1

cell = mp.Vector3(sx,sy,0)

geometry = []
geometry += [mp.Block(center=mp.Vector3(0,-Ny*3**0.5/4), material=mp.Medium(epsilon=eps_InSb_amor+(eps_InSb_crys-eps_InSb_amor)*0),
                     size=mp.Vector3(Nx,Ny*3**0.5/2)),
            mp.Block(center=mp.Vector3(0, Ny*3**0.5/4), material=mp.Medium(epsilon=eps_InSb_amor+(eps_InSb_crys-eps_InSb_amor)*1),
                     size=mp.Vector3(Nx,Ny*3**0.5/2))]


for i in range(Nx+1):
    for j in range(Ny):
        center = (-Nx/2+i+1/2*(j%2),(-Ny+j)*3**0.5/2)
        geometry += [mp.Cylinder(radius = 0.001*int(sys.argv[1]), center = mp.Vector3(center[0],center[1]), material=mp.Medium(epsilon=eps_GST326_amor+(eps_GST326_crys-eps_GST326_amor)*0)),
                     mp.Cylinder(radius = 0.001*int(sys.argv[2]), center = mp.Vector3(center[0],center[1]), material=mp.Medium(epsilon=1))
                    ]
for i in range(Nx+1):
    for j in range(Ny+1):
        center = (-Nx/2+i+1/2*(j%2),(j)*3**0.5/2)
        geometry += [mp.Cylinder(radius = 0.001*int(sys.argv[1]), center = mp.Vector3(center[0],center[1]), material=mp.Medium(epsilon=eps_GST326_amor+(eps_GST326_crys-eps_GST326_amor)*1)),
                     mp.Cylinder(radius = 0.001*int(sys.argv[2]), center = mp.Vector3(center[0],center[1]), material=mp.Medium(epsilon=1))
                    ]

sources = [mp.Source(mp.GaussianSource(frequency=fcen, fwidth = fwidth, width = 20),
                     mp.Hz,
                     mp.Vector3(-5,0))
                     ]

pml_layers = [mp.PML(plmWidth)]

resolution = 48

sim = mp.Simulation(cell_size=cell,
                    boundary_layers=pml_layers,
                    geometry=geometry,
                    default_material = mp.Medium(index=1),
                    sources=sources,
                    resolution=resolution)

sim.use_output_directory("wave_in_bulk")

def myhello(thisSim):
  # print(-0.4*sx+sim.meep_time()*v)
    if(thisSim.meep_time() > 190):
        movingSources = []
        thisSim.change_sources(movingSources)
    return 

sim.run(
        mp.at_beginning(mp.output_epsilon),
        mp.to_appended("Hz",mp.at_every(1,mp.output_hfield_z)),
        myhello,
        until=runTime)

# eps_data = sim.get_array(center=mp.Vector3(), size=cell, component=mp.Dielectric)
# plt.figure(dpi=100)
# plt.imshow(eps_data.transpose(), interpolation='spline36', cmap='binary')
# plt.axis('off')
# plt.savefig('./wave_in_band_gap/eps.png')

