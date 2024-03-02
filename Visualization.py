# -*- coding: utf-8 -*-
"""
Created on Sat Feb  3 19:21:49 2024

@author: Morteza
"""

import spectral
import numpy as np


# Specify the folder path containing the ENVI files
input_folder = "C:/Users/Morteza/OneDrive/Desktop/PhD/dataN"

# Specify the base name of the ENVI header file (without extension)
base_name = "Seurat_BEFORE"

# Construct the full paths for the header and data files
header_file = f"{input_folder}/{base_name}.hdr"
data_file = f"{input_folder}/{base_name}.dat"  # Or use the correct extension based on your data format

# Read the hyperspectral image using spectral
spectral_image = spectral.open_image(header_file)

# Access the data as a NumPy array
hyperspectral_data = spectral_image.load()
#hyperspectral_data = hyperspectral_data[:, :, 0:120]
# Print the shape of the data array
print("Shape of hyperspectral data:", hyperspectral_data.shape)

m = hyperspectral_data.shape;
#Simple visualization
viz1 = np.zeros((m[0], m[1],3))
viz1[:,:,0] = hyperspectral_data[:,:,100].squeeze();
viz1[:,:,1] = hyperspectral_data[:,:,70].squeeze();
viz1[:,:,2] = hyperspectral_data[:,:,40].squeeze();
import matplotlib.pyplot as plt
plt.imshow(viz1)
plt.show()
#colorimetric visualization


#Interpolating the standard data of standard illuminant and 
#standard observer to coincide with the wavelengths that
#our hyperspectral image has
from scipy.interpolate import interp1d

wavelengths = np.array(hyperspectral_data.bands.centers)
wavelengths = wavelengths[0:120]
# Wavelengths from 400 to 700 nm with 5 nm interval
x = np.arange(400, 701, 5)
x = x.reshape((61, 1))
x = x.ravel()
#D65 standard illuminant
from Accessories import D65
y = D65;
y = y.ravel()
# Create an interpolation function
interp_function = interp1d(x, y, kind='linear', fill_value="extrapolate")

# Define points where you want to interpolate
x_new = wavelengths

# Perform interpolation
D65n = interp_function(x_new)

#Standard oberver
from Accessories import xyz
y = xyz;


# Create interpolation functions for each column
interp_func_0 = interp1d(x, xyz[:, 0], kind='linear', fill_value='extrapolate')
interp_func_1 = interp1d(x, xyz[:, 1], kind='linear', fill_value='extrapolate')
interp_func_2 = interp1d(x, xyz[:, 2], kind='linear', fill_value='extrapolate')

# Use the interpolation functions to get the interpolated values
interpolated_values_0 = interp_func_0(x_new)
interpolated_values_1 = interp_func_1(x_new)
interpolated_values_2 = interp_func_2(x_new)

# Combine the interpolated values into the final matrix
xyzn = np.column_stack((interpolated_values_0, interpolated_values_1, interpolated_values_2))
RI_Slotff = hyperspectral_data[:,:,0:120]
RI_Slotff = RI_Slotff.reshape((670*1062,120))

#conversian of Reflectance to CIEXYZ tristimulus values

Mul_temp=np.matmul((np.transpose(xyzn)),(np.diagflat(D65n)))
Mul_temp2=np.matmul(Mul_temp,(np.transpose(RI_Slotff)))

XYZ = 1/100*Mul_temp2

#Now, XYZ to sRGB
import colour
#from colour.models import RGB_COLOURSPACE_sRGB
import skimage.exposure as exposure
XYZ = (np.transpose(XYZ)).reshape((670,1062,3))
SRGB = colour.XYZ_to_sRGB(XYZ)
# Apply the contrast stretch
SRGB = exposure.rescale_intensity(SRGB)
plt.imshow(SRGB)
plt.show()