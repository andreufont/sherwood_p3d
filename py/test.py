import sherwood_simulation as she_sim
import measured_power
import model_density

# specify simulation
L_hMpc=160
n_part=2048
sim=she_sim.SherwoodSimulation(L_hMpc=L_hMpc,n_part=n_part)

# specify snapshot and 3D grids
n_xy=1024
n_z=2048
snap_num=9
skewers=she_sim.Grid(simulation=sim,snapshot_num=snap_num,
            n_xy=n_xy,n_z=n_z,axis=0)

# read flux 1D power spectrum
flux_p1d=measured_power.get_measured_power(grid=skewers,power_type="p1d")

# read flux 3D power spectrum
flux_p3d=measured_power.get_measured_power(grid=skewers,power_type="p3d")
flux_p3d.plot_p3d(downsample_mu=5)

# read (real space) halo power spectrum
logMh_min=11.50
logMh_max=14.0
halos_no_rsd=she_sim.HaloGrid(simulation=sim,snapshot_num=snap_num,
            logMh_min=logMh_min,logMh_max=logMh_max,
            n_xy=n_xy,n_z=n_z,axis=1,add_rsd=False)
halo_no_rsd_p3d=measured_power.get_measured_power(grid=halos_no_rsd,
            power_type="halo")
halo_no_rsd_p3d.plot_p3d(subtract_shot_noise=False)

# read (redshift space) halo power spectrum
halos_rsd=she_sim.HaloGrid(simulation=sim,snapshot_num=snap_num,
            logMh_min=logMh_min,logMh_max=logMh_max,
            n_xy=n_xy,n_z=n_z,axis=0,add_rsd=True)
halo_rsd_p3d=measured_power.get_measured_power(grid=halos_rsd,
            power_type="halo")
halo_rsd_p3d.plot_p3d(downsample_mu=5)

# read cross power spectrum
cross_p3d=measured_power.get_measured_power(grid=halos_rsd,
            power_type="cross")
cross_p3d.plot_p3d_no_grid(halodata=halo_rsd_p3d,lyadata=flux_p3d,
            downsample_mu=5)

# play with linear power 
linP_model = model_density.LinearDensityModel()
linP_model.linP_hMpc(z=2.4,k_hMpc=1.0)
