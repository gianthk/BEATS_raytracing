import numpy as np
import matplotlib
import matplotlib.pyplot as plt

def deg2rad(angle_deg):
    return angle_deg*(np.pi)/180

def rad2deg(angle_rad):
    return angle_rad*180/np.pi

stripe1_edge_theta_deg = np.array([-0.73166, -0.60108, -0.35582]) # deg; found with energy scans and metal foils (Mo, Pd, In). in our lab reference grazing angles are negative
stripe1_edge_theta = deg2rad(stripe1_edge_theta_deg)
stripe1_edge_en = np.array([20000., 24350., 33169.]) # eV - Mo, Pd, In
stripe1_edge_lambda = 1239.842/stripe1_edge_en # nm
stripe1_d_nominal = 2.5 # nm - d-spacing
n = np.arange(0.9999999, 1.00003, (1-0.9999999)) # average refractive index of the multilayer
d = np.arange(2.0, 2.6, 0.001) # effective d-spacing

theta_refraction = np.zeros([n.size, d.size, stripe1_edge_lambda.size], dtype='float')
theta_refraction_diff = np.zeros(theta_refraction.shape, dtype='float')
for ni in range(0, n.size):
    for di in range(0, d.size):
        for li in range(0, stripe1_edge_lambda.size):
            # compute refraction-corrected Bragg angle for a n-d space
            theta_refraction[ni, di, li] = -np.arccos(np.sqrt(n[ni]**2-(stripe1_edge_lambda[li]/(2*d[di]))**2))

            # difference with experimental values
            theta_refraction_diff[ni, di, li] = np.abs(theta_refraction[ni, di, li] - stripe1_edge_theta[li])

# we want to minimize the distance between predicted and measured values (for all 3 angles)
theta_refraction_std = np.std(theta_refraction_diff, 2)
theta_refraction_var = np.var(theta_refraction_diff, 2)
theta_refraction_sum = np.sum(theta_refraction_diff, 2)

print(d[221])
print(n[164])

plt.xlabel("d-spacing [deg]")
plt.ylabel("n")

im = plt.imshow(theta_refraction_sum, cmap=plt.cm.RdBu) # extent=[-3, 3, -3, 3]
plt.colorbar(im)
plt.show()
