# Module for yield calculations
# Yun-Hsuan Abby Lee <abby6677@berkeley.edu>

# For those who don't have curie, go to https://jtmorrell.github.io/curie/build/html/index.html
import curie as ci
import numpy as np
import scipy.constants as sc
import matplotlib.pyplot as plt

# *************************************************************************
# *************************************************************************
#                                  INPUTS
# *************************************************************************
# *************************************************************************

rx = ci.Reaction('54FE(p,a)51MNg')     # reaction of interest and xs library of choice(optional)
target_thickness = 0.5E-3              # m
d_thickness = 1E-6                     # m
beam_type = 'p'                        # options are 'p', 'd'
beam_energy = 16                       # MeV
beam_current = 40                      # uA
irradiation_length = 3600              # s  - half-lives are being pulled in in units of s

# *************************************************************************
# *************************************************************************
#                              END OF INPUTS
# *************************************************************************
# *************************************************************************

print(rx.library.name)                 # best xs library
ip = ci.Isotope(rx.target)             # target isotope
el = ci.Element(ip.element)            # target element
product = rx.product                   # product isotope
# Energies

try:
    slices = int(target_thickness/d_thickness)
except:
    print('Please adjust d_thickness :(')

thickness_array = np.linspace(0,target_thickness,slices)*1.0E+3    # mm

energy = np.zeros(slices)
energy[0] = beam_energy
for i in range(slices-1):
    energy_loss = el.S(energy[i], particle=beam_type)                     # MeV/cm
    if energy[i] >= energy_loss * d_thickness * 100:
        energy[i+1] = energy[i] - energy_loss * d_thickness * 100   # 100 for m to cm
    else:
        energy[i+1] = 0

plt.plot(thickness_array,energy)
plt.xlabel('Target Thickness (mm)')
plt.ylabel('Energy (MeV)')
plt.title('Energies through target')
plt.savefig('Energies through target.png')

# Cross Sections

xs = rx.interpolate(energy)             #mb

plt.figure()
plt.plot(energy,xs)
plt.xlabel('energy (MeV)')
plt.ylabel('xs (mb)')
plt.title('XS')
plt.savefig('XS.png')

# Activation

particle_I = beam_current*1E-6/sc.e

N = el.density * d_thickness * 100 / ip.mass * sc.Avogadro
Reaction = np.multiply(N * particle_I * 1E-27, xs)              # atom/s

ip_product = ci.Isotope(product)
product_lambda = np.log(2)/ip_product.half_life()
activity = np.multiply(Reaction, (1 - np.exp(- product_lambda * irradiation_length)))   #Bq
activity_converted = activity/1E+6  #MBq

plt.figure()
plt.plot(thickness_array,activity_converted)
plt.xlabel('Target Thickness (mm)')
plt.ylabel('Activity (MBq)')
plt.title('Thin target activity')
plt.savefig('Thin target activity.png')

cumulative_activity = np.zeros(slices)
for i in range(slices-1):
    cumulative_activity[i+1] = cumulative_activity[i] + activity_converted[i+1]

plt.figure()
plt.plot(thickness_array,cumulative_activity)
plt.xlabel('Target Thickness (mm)')
plt.ylabel('Activity (MBq)')
plt.title('Thick Target Activity')
plt.savefig('Thick Target Activity.png')

# Yield

product_yield = np.divide(np.divide(activity,3.7E+7),beam_current*irradiation_length/3600)

plt.figure()
plt.plot(thickness_array,product_yield)
plt.xlabel('Target Thickness (mm)')
plt.ylabel('Yield (mCi/uAh)')
plt.title('Thin Target Yield')
plt.savefig('Thin Target Yield.png')

cumulative_yield = np.zeros(slices)
for i in range(slices-1):
    cumulative_yield[i+1] = cumulative_yield[i] + product_yield[i+1]
    
plt.figure()
plt.plot(thickness_array,cumulative_yield)
plt.xlabel('Target Thickness (mm)')
plt.ylabel('Thick Target Yield (mCi/uAh)')
plt.title('Thick Target Yield')
plt.savefig('Thick Target Yield.png')

# Et Viola!

# save txt files for TYtool validation
np.savetxt('thickness.csv', thickness_array, delimiter=',')
np.savetxt('thick_target_activity.csv', cumulative_activity, delimiter=',')
np.savetxt('thick_target_yield.csv', cumulative_yield, delimiter=',')
