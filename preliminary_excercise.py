import numpy as np
import matplotlib.pyplot as plt


data_array = np.loadtxt('jla_mub.txt',skiprows =1)



def etaFunc(a, omegam):
    s = ((1.0-omegam)/omegam) ** (1.0/3.0)
    part1 = 2.0 * np.sqrt((s ** 3.0) + 1) 
    part2 = ((a**(-4.0)) - (0.1540*s*(a**(-3.0))) + (0.4304*((s/a)**2.0)) + (0.19097*(s**3.0)/a) + (0.066941*(s**4.0)))**(-1.0/8.0)
    return (part1*part2)

def luminosityDistance(z, omegam):
    DL = 3000.0 * (1.0 + z) *(etaFunc(1.0,omegam) - etaFunc(1.0/(1.0+z),omegam))
    return DL

def distanceModulus(z,h,omegam):
    mu = 25.0 - (5.0*np.log10(h)) + (5.0*np.log10(luminosityDistance(z,omegam)))
    return mu

z = np.linspace(0,2,100)
print(z)
omegam_02 = [distanceModulus(x,0.7,0.2) for x in z]
omegam_03 = [distanceModulus(x,0.7,0.3) for x in z]
omegam_04 = [distanceModulus(x,0.7,0.4) for x in z]
omegam_05 = [distanceModulus(x,0.7,0.5) for x in z]

plt.figure()
plt.plot(z,omegam_02,label = '$\\Omega_M = 0.2$')
plt.plot(z,omegam_03,label = '$\\Omega_M = 0.3$')
plt.plot(z,omegam_04,label = '$\\Omega_M = 0.4$')
plt.plot(z,omegam_05,label = '$\\Omega_M = 0.5$')
plt.plot(data_array[:,0],data_array[:,1],label='JLA')
plt.legend()
plt.show()

sample_z = []
sample_mu =[]
h_06 = [distanceModulus(x,0.6,0.3) for x in z]
h_07 = [distanceModulus(x,0.7,0.3) for x in z]
h_08 = [distanceModulus(x,0.8,0.3) for x in z]
sample_mu_error = []
for i in range (0,20, 1):
    genz = np.random.uniform(low=0.0,high=2.0)
    sample_z.append(genz)


sample_z = sorted(sample_z)

for i in range (0,20,1):
    sample_mu.append(distanceModulus(sample_z[i],0.7,0.3))

    sample_mu_error.append(np.random.normal(loc = 0.0, scale = 0.1))

print(sample_mu,sample_z,sample_mu_error)

plt.figure()
plt.errorbar(sample_z,sample_mu,yerr=sample_mu_error,label='sampled supernovae')
plt.plot(z,h_06,label = 'h = 0.6')
plt.plot(z,h_07,label = 'h = 0.7')
plt.plot(z,h_08,label = 'h = 0.8')
plt.legend()
plt.show()