#------------------------------------------------------------------------------------------------------------------------------------------

#CTA200H Assignment - Python Script
#Ayush Pandhi (1003227457)
#May 19, 2019

#------------------------------------------------------------------------------------------------------------------------------------------

#Importing required modules
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from mpl_toolkits.mplot3d import Axes3D

#Loading in ATNF catalog from a .txt file
d3_num, d3_name, d3_jname, d3_gl, d3_gb, d3_dm, d3_rm, d3_dist, d3_zz, d3_xx, d3_yy = np.genfromtxt('ATNF_Data.txt', unpack=True)
d3_name, d3_jname = np.genfromtxt('ATNF_Data.txt', usecols=(1,2), unpack=True, dtype=str)

#Loading a second set which will be used later
d3_num2, d3_name2, d3_jname2, d3_gl2, d3_gb2, d3_dm2, d3_rm2, d3_dist2, d3_zz2, d3_xx2, d3_yy2 = np.genfromtxt('ATNF_Data.txt', unpack=True)
d3_name2, d3_jname2 = np.genfromtxt('ATNF_Data.txt', usecols=(1,2), unpack=True, dtype=str)

#------------------------------------------------------------------------------------------------------------------------------------------

#Finding the indicies for which RM is not given (I've set no RM values to 0 in the .txt file)
rm_list = np.where(d3_rm == 0)
rm_idx = rm_list[0]

#Using np.delete to remove the indicies from the arrays
d3_num = np.delete(d3_num, rm_idx)
d3_name = np.delete(d3_name, rm_idx)
d3_jname = np.delete(d3_jname, rm_idx)
d3_gl = np.delete(d3_gl, rm_idx)
d3_gb = np.delete(d3_gb, rm_idx)
d3_dm = np.delete(d3_dm, rm_idx)
d3_rm = np.delete(d3_rm, rm_idx)
d3_dist = np.delete(d3_dist, rm_idx)
d3_xx = np.delete(d3_xx, rm_idx)
d3_yy = np.delete(d3_yy, rm_idx)
d3_zz = np.delete(d3_zz, rm_idx)

#------------------------------------------------------------------------------------------------------------------------------------------

#Some plots just to get an idea of what's going on, these aren't required for the assignment and thus are not saved as pdfs

#Galactic Projection in x-y coordinate frame with height restrictions
plt.figure(figsize=(10,8))
plt.plot(d3_xx, d3_yy, '.')
plt.xlabel('x [kpc]', fontsize=10)
plt.ylabel('y [kpc]', fontsize=10)
plt.title('Galactic Coordinate Projection of x-y axis with RM', fontsize=15)
plt.show()

#Galactic Projection in y-z coordinate frame with height restrictions
plt.figure(figsize=(10,8))
plt.plot(d3_yy, d3_zz, '.')
plt.xlabel('y [kpc]', fontsize=10)
plt.ylabel('z [kpc]', fontsize=10)
plt.title('Galactic Coordinate Projection of y-z axis with RM', fontsize=15)
plt.show()

#Galactic Projection in x-z coordinate frame with height restrictions
plt.figure(figsize=(10,8))
plt.plot(d3_xx, d3_zz, '.')
plt.xlabel('x [kpc]', fontsize=10)
plt.ylabel('z [kpc]', fontsize=10)
plt.title('Galactic Coordinate Projection of x-z axis with RM', fontsize=15)
plt.show()

#Plotting in 3 dimensions with height restrictions
fig = plt.figure(figsize=(10,8))
ax = Axes3D(fig)
ax.scatter(xs=d3_xx, ys=d3_yy, zs=d3_zz, c = 'b', marker='.', linewidth=3)
ax.set_xlabel('x [kpc]', fontsize=10)
ax.set_ylabel('y [kpc]', fontsize=10)
ax.set_zlabel('z [kpc]', fontsize=10)
ax.set_title('Galactic Coordinate Projection in 3 Dimensions with RM', fontsize=15)
plt.show()

#------------------------------------------------------------------------------------------------------------------------------------------

#Same plots as above but with a colour map for RM, again these are not required for the assignment and thus are not saved as pdfs

#Galactic Projection in x-y coordinate frame with RM colormap
plt.figure(figsize=(10,8))
plt.scatter(d3_xx, d3_yy, c=d3_rm, cmap='magma', vmin=-300, vmax=300, marker='.')
plt.xlabel('x [kpc]', fontsize=10)
plt.ylabel('y [kpc]', fontsize=10)
plt.title('Galactic Coordinate Projection of x-y axis with RM', fontsize=15)
plt.colorbar()
plt.show()

#Galactic Projection in y-z coordinate frame with RM colormap
plt.figure(figsize=(10,8))
plt.scatter(d3_yy, d3_zz, c=d3_rm, cmap='magma', vmin=-300, vmax=300, marker='.')
plt.xlabel('y [kpc]', fontsize=10)
plt.ylabel('z [kpc]', fontsize=10)
plt.title('Galactic Coordinate Projection of y-z axis with RM', fontsize=15)
plt.colorbar()
plt.show()

#Galactic Projection in x-z coordinate frame with RM colormap
plt.figure(figsize=(10,8))
plt.scatter(d3_xx, d3_zz, c=d3_rm, cmap='magma', vmin=-300, vmax=300, marker='.')
plt.xlabel('x [kpc]', fontsize=10)
plt.ylabel('z [kpc]', fontsize=10)
plt.title('Galactic Coordinate Projection of x-z axis with RM', fontsize=15)
plt.colorbar()
plt.show()

#Plotting in 3 dimensions with height restrictions
fig = plt.figure(figsize=(10,8))
ax = Axes3D(fig)
ax.scatter(xs=d3_xx, ys=d3_yy, zs=d3_zz, c=d3_rm, cmap='magma', vmin=-300, vmax=300, marker='.', linewidth=3)
ax.set_xlabel('x [kpc]', fontsize=10)
ax.set_ylabel('y [kpc]', fontsize=10)
ax.set_zlabel('z [kpc]', fontsize=10)
ax.set_title('Galactic Coordinate Projection in 3 Dimensions with RM Colormap', fontsize=15)
plt.show()

#------------------------------------------------------------------------------------------------------------------------------------------

#Adding a condition for pulsars to be in the galactic plane
h_list = np.where(np.abs(d3_zz) > 0.5)
h_idx = h_list[0]

#Using np.delete to remove the indicies from the arrays
d3h_num = np.delete(d3_num, h_idx)
d3h_name = np.delete(d3_name, h_idx)
d3h_jname = np.delete(d3_jname, h_idx)
d3h_gl = np.delete(d3_gl, rm_idx)
d3h_gb = np.delete(d3_gb, rm_idx)
d3h_dm = np.delete(d3_dm, h_idx)
d3h_rm = np.delete(d3_rm, h_idx)
d3h_dist = np.delete(d3_dist, h_idx)
d3h_xx = np.delete(d3_xx, h_idx)
d3h_yy = np.delete(d3_yy, h_idx)
d3h_zz = np.delete(d3_zz, h_idx)

#------------------------------------------------------------------------------------------------------------------------------------------

#Same as the last set of plots but with height restrcition, again these are not required for the assignment and thus are not saved as pdfs

#Galactic Projection in x-y coordinate frame with RM colormap and z < 0.5kpc
plt.figure(figsize=(10,8))
plt.scatter(d3h_xx, d3h_yy, c=d3h_rm, cmap='magma', vmin=-300, vmax=300, marker='.')
plt.xlabel('x [kpc]', fontsize=10)
plt.ylabel('y [kpc]', fontsize=10)
plt.title('Galactic Coordinate Projection of x-y axis with RM and z < 0.5kpc', fontsize=15)
plt.colorbar()
plt.show()

#Galactic Projection in y-z coordinate frame with RM colormap and z < 0.5kpc
plt.figure(figsize=(10,8))
plt.scatter(d3h_yy, d3h_zz, c=d3h_rm, cmap='magma', vmin=-300, vmax=300, marker='.')
plt.xlabel('y [kpc]', fontsize=10)
plt.ylabel('z [kpc]', fontsize=10)
plt.title('Galactic Coordinate Projection of y-z axis with RM and z < 0.5kpc', fontsize=15)
plt.colorbar()
plt.show()

#Galactic Projection in x-z coordinate frame with RM colormap and z < 0.5kpc
plt.figure(figsize=(10,8))
plt.scatter(d3h_xx, d3h_zz, c=d3h_rm, cmap='magma', vmin=-300, vmax=300, marker='.')
plt.xlabel('x [kpc]', fontsize=10)
plt.ylabel('z [kpc]', fontsize=10)
plt.title('Galactic Coordinate Projection of x-z axis with RM and z < 0.5kpc', fontsize=15)
plt.colorbar()
plt.show()

#Plotting in 3 dimensions with height restrictions and RM colormap
fig = plt.figure(figsize=(10,8))
ax = Axes3D(fig)
ax.scatter(xs=d3h_xx, ys=d3h_yy, zs=d3h_zz, c=d3h_rm, cmap='magma', vmin=-300, vmax=300, marker='.', linewidth=3)
ax.set_xlabel('x [kpc]', fontsize=10)
ax.set_ylabel('y [kpc]', fontsize=10)
ax.set_zlabel('z [kpc]', fontsize=10)
ax.set_title('Galactic Coordinate Projection in 3 Dimensions with RM Colormap and z > 0.5kpc', fontsize=15)
plt.show()

#------------------------------------------------------------------------------------------------------------------------------------------

#Plotting rotation measure and dispersion measure relationship
plt.figure(figsize=(10,8))
plt.plot(d3_dm, d3_rm, '.')
plt.xlabel('Dispersion Measure [$cm^{-3}$ $pc$]', fontsize=10)
plt.ylabel('Rotation Measure [$rad$ $m^{-2}$]', fontsize=10)
plt.title('DM vs. RM for Milky Way Pulsars', fontsize=15)
plt.xlim(-30, 1200)
plt.ylim(-4000, 4000)
plt.savefig('DM_RM_plot.pdf')
plt.show()

#Plotting it again but with the height restriction applied
plt.figure(figsize=(10,8))
plt.plot(d3h_dm, d3h_rm, '.')
plt.xlabel('Dispersion Measure [$cm^{-3}$ $pc$]', fontsize=10)
plt.ylabel('Rotation Measure [$rad$ $m^{-2}$]', fontsize=10)
plt.title('DM vs. RM for Milky Way Pulsars z < 0.5kpc', fontsize=15)
plt.xlim(-30, 1200)
plt.ylim(-4000, 4000)
plt.savefig('DM_RM_hredc_plot.pdf')
plt.show()

#------------------------------------------------------------------------------------------------------------------------------------------

#Computing parallel component of magnetic field both with and without height restriction
B_par = 1.232*(d3_rm/d3_dm)
B_par2 = 1.232*(d3h_rm/d3h_dm)

#Plotting parallel magnetic field with galactic longitude/latitude
plt.figure(figsize=(10,8))
plt.plot(d3_gl, B_par, '.')
plt.xlabel('Galactic Longitude [degrees]', fontsize=10)
plt.ylabel('Parallel Magnetic Field [$\mu G$]', fontsize=10)
plt.title('Parallel Magnetic Field vs. Galactic Longitude', fontsize=15)
plt.ylim(-15, 15)
plt.savefig('GMF_Long_plot.pdf')
plt.show()

#Plotting parallel magnetic field with galactic longitude/latitude
plt.figure(figsize=(10,8))
plt.plot(d3_gb, B_par, '.')
plt.xlabel('Galactic Latitude [degrees]', fontsize=10)
plt.ylabel('Parallel Magnetic Field [$\mu G$]', fontsize=10)
plt.title('Parallel Magnetic Field vs. Galactic Latitude', fontsize=15)
plt.ylim(-15, 15)
plt.savefig('GMF_Lat_plot.pdf')
plt.show()

#Printing average galactic magnetic field value for this data set
print('Average GMF:', np.mean(B_par), '[$\mu G$]')

#------------------------------------------------------------------------------------------------------------------------------------------

#For these plots I had to flip the y axis and move it up by 0.5kpc as ATNF uses 8.5kpc as the solar circle radius while the Spitzer image uses -8.0kpc

#Loading a Milky Way map as an overlay for data
milkyway = mpimg.imread("MW_grid_crop.jpg")

#Overplotting data with known Milky Way structure
plt.figure(figsize=(15,15))
plt.imshow(milkyway, aspect='auto', extent=[-20, 20, -20, 18])
plt.scatter(d3_xx, (-d3_yy+0.5), c=d3_rm, cmap='Greys', vmin=-0.1, vmax=0.1, marker='o', s=12)
plt.plot(0, -8.0, 'r+', markersize=20, mew=1.5)
plt.xlabel('x [kpc]', fontsize=10)
plt.ylabel('y [kpc]', fontsize=10)
plt.title('Data Set with RM Overplotted on Milky Way Structure', fontsize=15)
plt.xlim(-20, 20)
plt.ylim(-20, 18)
plt.savefig('Spitzer_binary_plot.pdf')
plt.show()

#Overplotting data with known Milky Way structure (full spectrum)
plt.figure(figsize=(20,15))
plt.imshow(milkyway, aspect='auto', extent=[-20, 20, -20, 18])
plt.scatter(d3_xx, (-d3_yy+0.5), c=d3_rm, cmap='hot', vmin=-500, vmax=500, marker='o', s=12)
plt.plot(0, -8.0, 'w+', markersize=20, mew=1.5)
plt.xlabel('x [kpc]', fontsize=10)
plt.ylabel('y [kpc]', fontsize=10)
plt.title('Set with RM Overplotted on Milky Way Structure (full colourmap)', fontsize=15)
plt.xlim(-20, 20)
plt.ylim(-20, 18)
plt.colorbar(aspect=30, label='Rotation Measure [rad m^-2]')
plt.savefig('Spitzer_cmap_plot.pdf')
plt.show()

#Overplotting data with known Milky Way structure for parallel magnetic field
plt.figure(figsize=(20,15))
plt.imshow(milkyway, aspect='auto', extent=[-20, 20, -20, 18])
plt.scatter(d3_xx, (-d3_yy+0.5), c=B_par, cmap='seismic', vmin=-3, vmax=3, marker='o', s=12)
plt.plot(0, -8.0, 'w+', markersize=20, mew=1.5)
plt.xlabel('x [kpc]', fontsize=10)
plt.ylabel('y [kpc]', fontsize=10)
plt.title('Data Set with GMF Overplotted on Milky Way Structure', fontsize=15)
plt.xlim(-20, 20)
plt.ylim(-20, 18)
plt.colorbar(aspect=30, label='Parallel Magnetic Field [$\mu G$]')
plt.savefig('Spitzer_GMF_plot.pdf')
plt.show()

#Overplotting reduced data with known Milky Way structure
plt.figure(figsize=(15,15))
plt.imshow(milkyway, aspect='auto', extent=[-20, 20, -20, 18])
plt.scatter(d3h_xx, (-d3h_yy+0.5), c=d3h_rm, cmap='Greys', vmin=-0.1, vmax=0.1, marker='o', s=12)
plt.plot(0, -8.0, 'r+', markersize=20, mew=1.5)
plt.xlabel('x [kpc]', fontsize=10)
plt.ylabel('y [kpc]', fontsize=10)
plt.title('REDUCED Data Set with RM Overplotted on Milky Way Structure', fontsize=15)
plt.xlim(-20, 20)
plt.ylim(-20, 18)
plt.savefig('Spitzer_binary_hreduc_plot.pdf')
plt.show()

#Overplotting reduced data with known Milky Way structure (full spectrum)
plt.figure(figsize=(20,15))
plt.imshow(milkyway, aspect='auto', extent=[-20, 20, -20, 18])
plt.scatter(d3h_xx, (-d3h_yy+0.5), c=d3h_rm, cmap='hot', vmin=-500, vmax=500, marker='o', s=12)
plt.plot(0, -8.0, 'w+', markersize=20, mew=1.5)
plt.xlabel('x [kpc]', fontsize=10)
plt.ylabel('y [kpc]', fontsize=10)
plt.title('REDUCED Data Set with RM Overplotted on Milky Way Structure (full colourmap)', fontsize=15)
plt.xlim(-20, 20)
plt.ylim(-20, 18)
plt.colorbar(aspect=30, label='Rotation Measure [rad m^-2]')
plt.savefig('Spitzer_cmap_hreduc_plot.pdf')
plt.show()

#Overplotting reduced data with known Milky Way structure for parallel magnetic field
plt.figure(figsize=(20,15))
plt.imshow(milkyway, aspect='auto', extent=[-20, 20, -20, 18])
plt.scatter(d3h_xx, (-d3h_yy+0.5), c=B_par2, cmap='seismic', vmin=-3, vmax=3, marker='o', s=12)
plt.plot(0, -8.0, 'w+', markersize=20, mew=1.5)
plt.xlabel('x [kpc]', fontsize=10)
plt.ylabel('y [kpc]', fontsize=10)
plt.title('REDUCED Data Set with GMF Overplotted on Milky Way Structure', fontsize=15)
plt.xlim(-20, 20)
plt.ylim(-20, 18)
plt.colorbar(aspect=30, label='Parallel Magnetic Field [$\mu G$]')
plt.savefig('Spitzer_GMF_hreduc_plot.pdf')
plt.show()

#------------------------------------------------------------------------------------------------------------------------------------------

#Individual LOS plots; these plots are not saved as pdfs and do not shed a lot of light on the overall results on their own

#Finding the indicies for longitude 30 degrees
l27_list = np.where(d3_gl < 27)
l27_idx = l27_list[0]
l33_list = np.where(d3_gl > 33)
l33_idx = l33_list[0]
l30_idx = np.hstack((l27_idx, l33_idx))

#Using np.delete to remove the indicies from the arrays
d3l30_num = np.delete(d3_num, l30_idx)
d3l30_name = np.delete(d3_name, l30_idx)
d3l30_jname = np.delete(d3_jname, l30_idx)
d3l30_gl = np.delete(d3_gl, l30_idx)
d3l30_gb = np.delete(d3_gb, l30_idx)
d3l30_dm = np.delete(d3_dm, l30_idx)
d3l30_rm = np.delete(d3_rm, l30_idx)
d3l30_dist = np.delete(d3_dist, l30_idx)
d3l30_xx = np.delete(d3_xx, l30_idx)
d3l30_yy = np.delete(d3_yy, l30_idx)
d3l30_zz = np.delete(d3_zz, l30_idx)

#Computing parallel magnetic field for these pulsars
Bl30_par = 1.232*(d3l30_rm/d3l30_dm)

#Plotting
plt.figure(figsize=(10,8))
plt.plot(d3l30_dist, Bl30_par, '.')
plt.xlabel('Distance [kpc]', fontsize=10)
plt.ylabel('Parallel Magnetic Field []', fontsize=10)
plt.title('$l$=30 Magnetic Field vs. Distance', fontsize=15)
plt.show()

#Finding the indicies for longitude 330 degrees
l327_list = np.where(d3_gl < 327)
l327_idx = l327_list[0]
l333_list = np.where(d3_gl > 333)
l333_idx = l333_list[0]
l330_idx = np.hstack((l327_idx, l333_idx))

#Using np.delete to remove the indicies from the arrays
d3l330_num = np.delete(d3_num, l330_idx)
d3l330_name = np.delete(d3_name, l330_idx)
d3l330_jname = np.delete(d3_jname, l330_idx)
d3l330_gl = np.delete(d3_gl, l330_idx)
d3l330_gb = np.delete(d3_gb, l330_idx)
d3l330_dm = np.delete(d3_dm, l330_idx)
d3l330_rm = np.delete(d3_rm, l330_idx)
d3l330_dist = np.delete(d3_dist, l330_idx)
d3l330_xx = np.delete(d3_xx, l330_idx)
d3l330_yy = np.delete(d3_yy, l330_idx)
d3l330_zz = np.delete(d3_zz, l330_idx)

#Computing parallel magnetic field for these pulsars
Bl330_par = 1.232*(d3l330_rm/d3l330_dm)

#Plotting
plt.figure(figsize=(10,8))
plt.plot(d3l330_dist, Bl330_par, '.')
plt.xlabel('Distance [kpc]', fontsize=10)
plt.ylabel('Parallel Magnetic Field []', fontsize=10)
plt.title('$l$=330 Magnetic Field vs. Distance', fontsize=15)
plt.show()

#Finding the indicies for longitude 310 degrees
l307_list = np.where(d3_gl < 307)
l307_idx = l307_list[0]
l313_list = np.where(d3_gl > 313)
l313_idx = l313_list[0]
l310_idx = np.hstack((l307_idx, l313_idx))

#Using np.delete to remove the indicies from the arrays
d3l310_num = np.delete(d3_num, l310_idx)
d3l310_name = np.delete(d3_name, l310_idx)
d3l310_jname = np.delete(d3_jname, l310_idx)
d3l310_gl = np.delete(d3_gl, l310_idx)
d3l310_gb = np.delete(d3_gb, l310_idx)
d3l310_dm = np.delete(d3_dm, l310_idx)
d3l310_rm = np.delete(d3_rm, l310_idx)
d3l310_dist = np.delete(d3_dist, l310_idx)
d3l310_xx = np.delete(d3_xx, l310_idx)
d3l310_yy = np.delete(d3_yy, l310_idx)
d3l310_zz = np.delete(d3_zz, l310_idx)

#Computing parallel magnetic field for these pulsars
Bl310_par = 1.232*(d3l310_rm/d3l310_dm)

#Plotting
plt.figure(figsize=(10,8))
plt.plot(d3l310_dist, Bl310_par, '.')
plt.xlabel('Distance [kpc]', fontsize=10)
plt.ylabel('Parallel Magnetic Field []', fontsize=10)
plt.title('$l$=310 Magnetic Field vs. Distance', fontsize=15)
plt.show()

#Finding the indicies for longitude 50 degrees
l47_list = np.where(d3_gl < 47)
l47_idx = l47_list[0]
l53_list = np.where(d3_gl > 53)
l53_idx = l53_list[0]
l50_idx = np.hstack((l47_idx, l53_idx))

#Using np.delete to remove the indicies from the arrays
d3l50_num = np.delete(d3_num, l50_idx)
d3l50_name = np.delete(d3_name, l50_idx)
d3l50_jname = np.delete(d3_jname, l50_idx)
d3l50_gl = np.delete(d3_gl, l50_idx)
d3l50_gb = np.delete(d3_gb, l50_idx)
d3l50_dm = np.delete(d3_dm, l50_idx)
d3l50_rm = np.delete(d3_rm, l50_idx)
d3l50_dist = np.delete(d3_dist, l50_idx)
d3l50_xx = np.delete(d3_xx, l50_idx)
d3l50_yy = np.delete(d3_yy, l50_idx)
d3l50_zz = np.delete(d3_zz, l50_idx)

#Computing parallel magnetic field for these pulsars
Bl50_par = 1.232*(d3l50_rm/d3l50_dm)

#Plotting
plt.figure(figsize=(10,8))
plt.plot(d3l50_dist, Bl50_par, '.')
plt.xlabel('Distance [kpc]', fontsize=10)
plt.ylabel('Parallel Magnetic Field []', fontsize=10)
plt.title('$l$=50 Magnetic Field vs. Distance', fontsize=15)
plt.xlim(-0.5, 12)
plt.show()

#Finding the indicies for longitude 10 degrees
l7_list = np.where(d3_gl < 7)
l7_idx = l7_list[0]
l13_list = np.where(d3_gl > 13)
l13_idx = l13_list[0]
l10_idx = np.hstack((l7_idx, l13_idx))

#Using np.delete to remove the indicies from the arrays
d3l10_num = np.delete(d3_num, l10_idx)
d3l10_name = np.delete(d3_name, l10_idx)
d3l10_jname = np.delete(d3_jname, l10_idx)
d3l10_gl = np.delete(d3_gl, l10_idx)
d3l10_gb = np.delete(d3_gb, l10_idx)
d3l10_dm = np.delete(d3_dm, l10_idx)
d3l10_rm = np.delete(d3_rm, l10_idx)
d3l10_dist = np.delete(d3_dist, l10_idx)
d3l10_xx = np.delete(d3_xx, l10_idx)
d3l10_yy = np.delete(d3_yy, l10_idx)
d3l10_zz = np.delete(d3_zz, l10_idx)

#Computing parallel magnetic field for these pulsars
Bl10_par = 1.232*(d3l10_rm/d3l10_dm)

#Plotting
plt.figure(figsize=(10,8))
plt.plot(d3l10_dist, Bl10_par, '.')
plt.xlabel('Distance [kpc]', fontsize=10)
plt.ylabel('Parallel Magnetic Field []', fontsize=10)
plt.title('$l$=10 Magnetic Field vs. Distance', fontsize=15)
plt.xlim(-0.5, 12)
plt.show()

#Finding the indicies for longitude 350 degrees
l347_list = np.where(d3_gl < 347)
l347_idx = l347_list[0]
l353_list = np.where(d3_gl > 353)
l353_idx = l353_list[0]
l350_idx = np.hstack((l347_idx, l353_idx))

#Using np.delete to remove the indicies from the arrays
d3l350_num = np.delete(d3_num, l350_idx)
d3l350_name = np.delete(d3_name, l350_idx)
d3l350_jname = np.delete(d3_jname, l350_idx)
d3l350_gl = np.delete(d3_gl, l350_idx)
d3l350_gb = np.delete(d3_gb, l350_idx)
d3l350_dm = np.delete(d3_dm, l350_idx)
d3l350_rm = np.delete(d3_rm, l350_idx)
d3l350_dist = np.delete(d3_dist, l350_idx)
d3l350_xx = np.delete(d3_xx, l350_idx)
d3l350_yy = np.delete(d3_yy, l350_idx)
d3l350_zz = np.delete(d3_zz, l350_idx)

#Computing parallel magnetic field for these pulsars
Bl350_par = 1.232*(d3l350_rm/d3l350_dm)

#Plotting
plt.figure(figsize=(10,8))
plt.plot(d3l350_dist, Bl350_par, '.')
plt.xlabel('Distance [kpc]', fontsize=10)
plt.ylabel('Parallel Magnetic Field []', fontsize=10)
plt.title('$l$=350 Magnetic Field vs. Distance', fontsize=15)
plt.show()

#Finding the indicies for longitude 70 degrees
l67_list = np.where(d3_gl < 67)
l67_idx = l67_list[0]
l73_list = np.where(d3_gl > 73)
l73_idx = l73_list[0]
l70_idx = np.hstack((l67_idx, l73_idx))

#Using np.delete to remove the indicies from the arrays
d3l70_num = np.delete(d3_num, l70_idx)
d3l70_name = np.delete(d3_name, l70_idx)
d3l70_jname = np.delete(d3_jname, l70_idx)
d3l70_gl = np.delete(d3_gl, l70_idx)
d3l70_gb = np.delete(d3_gb, l70_idx)
d3l70_dm = np.delete(d3_dm, l70_idx)
d3l70_rm = np.delete(d3_rm, l70_idx)
d3l70_dist = np.delete(d3_dist, l70_idx)
d3l70_xx = np.delete(d3_xx, l70_idx)
d3l70_yy = np.delete(d3_yy, l70_idx)
d3l70_zz = np.delete(d3_zz, l70_idx)

#Computing parallel magnetic field for these pulsars
Bl70_par = 1.232*(d3l70_rm/d3l70_dm)

#Plotting
plt.figure(figsize=(10,8))
plt.plot(d3l70_dist, Bl70_par, '.')
plt.xlabel('Distance [kpc]', fontsize=10)
plt.ylabel('Parallel Magnetic Field []', fontsize=10)
plt.title('$l$=70 Magnetic Field vs. Distance', fontsize=15)
plt.show()

#Finding the indicies for longitude 290 degrees
l287_list = np.where(d3_gl < 287)
l287_idx = l287_list[0]
l293_list = np.where(d3_gl > 293)
l293_idx = l293_list[0]
l290_idx = np.hstack((l287_idx, l293_idx))

#Using np.delete to remove the indicies from the arrays
d3l290_num = np.delete(d3_num, l290_idx)
d3l290_name = np.delete(d3_name, l290_idx)
d3l290_jname = np.delete(d3_jname, l290_idx)
d3l290_gl = np.delete(d3_gl, l290_idx)
d3l290_gb = np.delete(d3_gb, l290_idx)
d3l290_dm = np.delete(d3_dm, l290_idx)
d3l290_rm = np.delete(d3_rm, l290_idx)
d3l290_dist = np.delete(d3_dist, l290_idx)
d3l290_xx = np.delete(d3_xx, l290_idx)
d3l290_yy = np.delete(d3_yy, l290_idx)
d3l290_zz = np.delete(d3_zz, l290_idx)

#Computing parallel magnetic field for these pulsars
Bl290_par = 1.232*(d3l290_rm/d3l290_dm)

#Plotting
plt.figure(figsize=(10,8))
plt.plot(d3l290_dist, Bl290_par, '.')
plt.xlabel('Distance [kpc]', fontsize=10)
plt.ylabel('Parallel Magnetic Field []', fontsize=10)
plt.title('$l$=290 Magnetic Field vs. Distance', fontsize=15)
plt.show()

#Finding the indicies for longitude 130 degrees
l127_list = np.where(d3_gl < 127)
l127_idx = l127_list[0]
l133_list = np.where(d3_gl > 133)
l133_idx = l133_list[0]
l130_idx = np.hstack((l127_idx, l133_idx))

#Using np.delete to remove the indicies from the arrays
d3l130_num = np.delete(d3_num, l130_idx)
d3l130_name = np.delete(d3_name, l130_idx)
d3l130_jname = np.delete(d3_jname, l130_idx)
d3l130_gl = np.delete(d3_gl, l130_idx)
d3l130_gb = np.delete(d3_gb, l130_idx)
d3l130_dm = np.delete(d3_dm, l130_idx)
d3l130_rm = np.delete(d3_rm, l130_idx)
d3l130_dist = np.delete(d3_dist, l130_idx)
d3l130_xx = np.delete(d3_xx, l130_idx)
d3l130_yy = np.delete(d3_yy, l130_idx)
d3l130_zz = np.delete(d3_zz, l130_idx)

#Computing parallel magnetic field for these pulsars
Bl130_par = 1.232*(d3l130_rm/d3l130_dm)

#Plotting
plt.figure(figsize=(10,8))
plt.plot(d3l130_dist, Bl130_par, '.')
plt.xlabel('Distance [kpc]', fontsize=10)
plt.ylabel('Parallel Magnetic Field []', fontsize=10)
plt.title('$l$=130 Magnetic Field vs. Distance', fontsize=15)
plt.show()

#Finding the indicies for longitude 150 degrees
l147_list = np.where(d3_gl < 147)
l147_idx = l147_list[0]
l153_list = np.where(d3_gl > 153)
l153_idx = l153_list[0]
l150_idx = np.hstack((l147_idx, l153_idx))

#Using np.delete to remove the indicies from the arrays
d3l150_num = np.delete(d3_num, l150_idx)
d3l150_name = np.delete(d3_name, l150_idx)
d3l150_jname = np.delete(d3_jname, l150_idx)
d3l150_gl = np.delete(d3_gl, l150_idx)
d3l150_gb = np.delete(d3_gb, l150_idx)
d3l150_dm = np.delete(d3_dm, l150_idx)
d3l150_rm = np.delete(d3_rm, l150_idx)
d3l150_dist = np.delete(d3_dist, l150_idx)
d3l150_xx = np.delete(d3_xx, l150_idx)
d3l150_yy = np.delete(d3_yy, l150_idx)
d3l150_zz = np.delete(d3_zz, l150_idx)

#Computing parallel magnetic field for these pulsars
Bl150_par = 1.232*(d3l150_rm/d3l150_dm)

#Plotting
plt.figure(figsize=(10,8))
plt.plot(d3l150_dist, Bl150_par, '.')
plt.xlabel('Distance [kpc]', fontsize=10)
plt.ylabel('Parallel Magnetic Field []', fontsize=10)
plt.title('$l$=150 Magnetic Field vs. Distance', fontsize=15)
plt.show()

#Finding the indicies for longitude 230 degrees
l227_list = np.where(d3_gl < 227)
l227_idx = l227_list[0]
l233_list = np.where(d3_gl > 233)
l233_idx = l233_list[0]
l230_idx = np.hstack((l227_idx, l233_idx))

#Using np.delete to remove the indicies from the arrays
d3l230_num = np.delete(d3_num, l230_idx)
d3l230_name = np.delete(d3_name, l230_idx)
d3l230_jname = np.delete(d3_jname, l230_idx)
d3l230_gl = np.delete(d3_gl, l230_idx)
d3l230_gb = np.delete(d3_gb, l230_idx)
d3l230_dm = np.delete(d3_dm, l230_idx)
d3l230_rm = np.delete(d3_rm, l230_idx)
d3l230_dist = np.delete(d3_dist, l230_idx)
d3l230_xx = np.delete(d3_xx, l230_idx)
d3l230_yy = np.delete(d3_yy, l230_idx)
d3l230_zz = np.delete(d3_zz, l230_idx)

#Computing parallel magnetic field for these pulsars
Bl230_par = 1.232*(d3l230_rm/d3l230_dm)

#Plotting
plt.figure(figsize=(10,8))
plt.plot(d3l230_dist, Bl230_par, '.')
plt.xlabel('Distance [kpc]', fontsize=10)
plt.ylabel('Parallel Magnetic Field []', fontsize=10)
plt.title('$l$=230 Magnetic Field vs. Distance', fontsize=15)
plt.show()

#------------------------------------------------------------------------------------------------------------------------------------------

#Overplotting various LOS with the Spitzer image to see if there are any field reversals
plt.figure(figsize=(20,15))
plt.imshow(milkyway, aspect='auto', extent=[-20, 20, -20, 18])
plt.scatter(d3l10_xx, (-d3l10_yy+0.5), c=Bl10_par, cmap='seismic', vmin=-3, vmax=3, marker='o', s=12)
plt.scatter(d3l30_xx, (-d3l30_yy+0.5), c=Bl30_par, cmap='seismic', vmin=-3, vmax=3, marker='o', s=12)
plt.scatter(d3l50_xx, (-d3l50_yy+0.5), c=Bl50_par, cmap='seismic', vmin=-3, vmax=3, marker='o', s=12)
plt.scatter(d3l70_xx, (-d3l70_yy+0.5), c=Bl70_par, cmap='seismic', vmin=-3, vmax=3, marker='o', s=12)
plt.scatter(d3l130_xx, (-d3l130_yy+0.5), c=Bl130_par, cmap='seismic', vmin=-3, vmax=3, marker='o', s=12)
plt.scatter(d3l150_xx, (-d3l150_yy+0.5), c=Bl150_par, cmap='seismic', vmin=-3, vmax=3, marker='o', s=12)
plt.scatter(d3l230_xx, (-d3l230_yy+0.5), c=Bl230_par, cmap='seismic', vmin=-3, vmax=3, marker='o', s=12)
plt.scatter(d3l290_xx, (-d3l290_yy+0.5), c=Bl290_par, cmap='seismic', vmin=-3, vmax=3, marker='o', s=12)
plt.scatter(d3l310_xx, (-d3l310_yy+0.5), c=Bl310_par, cmap='seismic', vmin=-3, vmax=3, marker='o', s=12)
plt.scatter(d3l330_xx, (-d3l330_yy+0.5), c=Bl330_par, cmap='seismic', vmin=-3, vmax=3, marker='o', s=12)
plt.scatter(d3l350_xx, (-d3l350_yy+0.5), c=Bl350_par, cmap='seismic', vmin=-3, vmax=3, marker='o', s=12)
plt.plot(0, -8.0, 'w+', markersize=20, mew=1.5)
plt.xlabel('x [kpc]', fontsize=10)
plt.ylabel('y [kpc]', fontsize=10)
plt.title('REDUCED Data Set with GMF Overplotted on Milky Way Structure (Longitude Analysis)', fontsize=15)
plt.xlim(-20, 20)
plt.ylim(-20, 18)
plt.colorbar(aspect=30, label='Parallel Magnetic Field [$\mu G$]')
plt.savefig('Spitzer_LOS_plot.pdf')
plt.show()

#Same plot but without the Spitzer image to more clearly see the LOS pulsars
plt.figure(figsize=(20,15))
plt.scatter(d3l10_xx, (-d3l10_yy+0.5), c=Bl10_par, cmap='seismic', vmin=-3, vmax=3, marker='o', s=12)
plt.scatter(d3l30_xx, (-d3l30_yy+0.5), c=Bl30_par, cmap='seismic', vmin=-3, vmax=3, marker='o', s=12)
plt.scatter(d3l50_xx, (-d3l50_yy+0.5), c=Bl50_par, cmap='seismic', vmin=-3, vmax=3, marker='o', s=12)
plt.scatter(d3l70_xx, (-d3l70_yy+0.5), c=Bl70_par, cmap='seismic', vmin=-3, vmax=3, marker='o', s=12)
plt.scatter(d3l130_xx, (-d3l130_yy+0.5), c=Bl130_par, cmap='seismic', vmin=-3, vmax=3, marker='o', s=12)
plt.scatter(d3l150_xx, (-d3l150_yy+0.5), c=Bl150_par, cmap='seismic', vmin=-3, vmax=3, marker='o', s=12)
plt.scatter(d3l230_xx, (-d3l230_yy+0.5), c=Bl230_par, cmap='seismic', vmin=-3, vmax=3, marker='o', s=12)
plt.scatter(d3l290_xx, (-d3l290_yy+0.5), c=Bl290_par, cmap='seismic', vmin=-3, vmax=3, marker='o', s=12)
plt.scatter(d3l310_xx, (-d3l310_yy+0.5), c=Bl310_par, cmap='seismic', vmin=-3, vmax=3, marker='o', s=12)
plt.scatter(d3l330_xx, (-d3l330_yy+0.5), c=Bl330_par, cmap='seismic', vmin=-3, vmax=3, marker='o', s=12)
plt.scatter(d3l350_xx, (-d3l350_yy+0.5), c=Bl350_par, cmap='seismic', vmin=-3, vmax=3, marker='o', s=12)
plt.plot(0, -8.0, 'w+', markersize=20, mew=1.5)
plt.xlabel('x [kpc]', fontsize=10)
plt.ylabel('y [kpc]', fontsize=10)
plt.title('REDUCED Data Set with GMF Overplotted on Milky Way Structure (Longitude Analysis)', fontsize=15)
plt.xlim(-20, 20)
plt.ylim(-20, 18)
plt.colorbar(aspect=30, label='Parallel Magnetic Field [$\mu G$]')
plt.savefig('LOS_plot.pdf')
plt.show()

#------------------------------------------------------------------------------------------------------------------------------------------

#Loading CHIME pulsar data
d4_name, d4_ra, d4_dec, d4_p0, d4_dm = np.genfromtxt('CHIME_Data.txt', unpack=True)
d4_name = np.genfromtxt('CHIME_Data.txt', usecols=0, dtype=str)

#Using a loop to find the matches between ATNF (w/ RM) data and CHIME data
print('CHIME Pulsars not in ATNF Catalog:')
matches = []
matches2 = []
matches_name = []
matches_name2 = []
for i in range(446):
    i_name = d4_name[i]
    i_match = np.where(d3_name == i_name)[0]
    i_match2 = np.where(d3_name2 == i_name)[0]
    if i_match.shape[0] == 1:
        matches.append(i_match[0])
        matches_name.append(i_name)
    if i_match2.shape[0] == 1:
        matches2.append(i_match2[0])
        matches_name2.append(i_name) 
    if i_match2.shape[0] == 0:
        print(i_name)
matches = np.array(matches)
matches_name = np.array(matches_name)
matches2 = np.array(matches2)
matches_name2 = np.array(matches_name2)

#------------------------------------------------------------------------------------------------------------------------------------------

#Plotting CHIME Pulsars with known/unknown RM values
plt.figure(figsize=(15,15))
plt.imshow(milkyway, aspect='auto', extent=[-20, 20, -20, 18])
for i in matches2:
    plt.scatter(d3_xx2[i], (-d3_yy2[i] + 0.5), marker='o', c='k', s=12)
for i in matches:
    plt.scatter(d3_xx[i], (-d3_yy[i] + 0.5), marker='o', c='r', s=12)
plt.plot(0, -8.0, 'w+', markersize=20, mew=1.5)
plt.xlabel('x [kpc]', fontsize=10)
plt.ylabel('y [kpc]', fontsize=10)
plt.title('CHIME Pulsars with Known RM (Red) amd Unknown RM (Black)', fontsize=15)
plt.xlim(-20, 20)
plt.ylim(-20, 18)
plt.savefig('CHIME_RM_plot.pdf')
plt.show()

#------------------------------------------------------------------------------------------------------------------------------------------
