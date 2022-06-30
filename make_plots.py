import matplotlib.pyplot as plt
import numpy as np
import pickle
import os, sys
from colour import Color

# Load files using pickle
with open('/data/p-one/dghuman/simulation/MCEq_tables/muon_table_dtheta_10.pkl', 'rb') as input1:
    table = pickle.load(input1)

with open('/data/p-one/dghuman/simulation/MCEq_tables/rate_muon_table_dtheta_10.pkl', 'rb') as input2:
    rate_table = pickle.load(input2)

#with open('/data/p-one/dghuman/simulation/MCEq_tables/live_time_zenith.pkl', 'rb') as input3:
#    livetime = pickle.load(input3)

with open('/data/p-one/dghuman/simulation/MCEq_tables/UofA_tables/fluxtable_as_list_1460.pkl', 'rb') as input4:
    uofa_table = pickle.load(input4)

with open('/data/p-one/dghuman/simulation/MCEq_tables/new_muon_dict.pkl', 'rb') as input5:
    full_muon_info = pickle.load(input5)

with open('/data/p-one/dghuman/simulation/MCEq_tables/muon_weights.pkl', 'rb') as input6:
    muon_weights = pickle.load(input6)

with open('/data/p-one/dghuman/simulation/MCEq_tables/muon_dict_0006.pkl', 'rb') as input7:
    muon_dict = pickle.load(input7)

with open('/data/p-one/dghuman/simulation/MCEq_tables/muon_dict_0006_triggered.pkl', 'rb') as input8:    
    muon_dict_trig = pickle.load(input8)

with open('/data/p-one/dghuman/simulation/MCEq_tables/triggered_muon_weights.pkl', 'rb') as input9:
    muon_weights_trig = pickle.load(input9)
    
e_grid = table['e_grid']
e_grid2 = uofa_table['e_grid']

#plot

l = [x for x in table.keys() if not isinstance(x, str)]
l2 = [x for x in rate_table.keys() if not isinstance(x, str)]
start = Color("#38f8ff")
end = Color("#ff38f8")
colors = list(start.range_to(end, len(l)))

l3 = [x for x in uofa_table.keys() if not isinstance(x,str)]

plt.figure(figsize=(4.2, 3))

for key in l:
    plt.loglog(e_grid, table[key], color=str(colors[l.index(key)].hex), ls='-', lw=1.5, label=r'$\theta =$ ' + str(key))
    plt.loglog(e_grid2, uofa_table[key], color=str(colors[l.index(key)].hex), ls='--', lw=1.5)#, label=r'$\theta =$ ' + str(key) + ' : mute')

plt.xlabel(r"$E_{\mu}$ [GeV]")
plt.ylabel(r"$\Phi_{\mu}$(cm$^{2}$ s sr GeV)$^{-1}$")
plt.title(r"Flux vs Energy")
plt.legend(loc='upper right',frameon=False,numpoints=1,fontsize='small')
plt.tight_layout()
plt.xlim(10,1e10)
plt.ylim(1e-17,10)
plt.grid(True)
    
plt.savefig('/data/p-one/dghuman/simulation/MCEq_tables/test_plots/E_vs_Flux_zenith.png', dpi=300)

plt.clf()

plt.figure(figsize=(4.2, 3))

for key in l2:
    plt.loglog(e_grid, rate_table[key], color=str(colors[l2.index(key)].hex), ls='-', lw=1.5, label=r'$\theta =$ ' + str(key))

plt.xlabel(r"$E_{\mu}$ [GeV]")
plt.ylabel(r"$\dot{N}_{\mu}$(s GeV)$^{-1}$")
plt.title(r"Rate vs Energy")
plt.legend(loc='upper right',frameon=False,numpoints=1,fontsize='medium')
plt.tight_layout()
plt.xlim(10,1e10)
plt.ylim(1e-17,10)
    
plt.savefig('/data/p-one/dghuman/simulation/MCEq_tables/test_plots/E_vs_Rate_zenith.png', dpi=300)

plt.clf()

n, bins, patches = plt.hist(x=np.array(full_muon_info['Energy']), bins=np.array(e_grid), alpha=1.0, edgecolor='black', linestyle='solid', histtype='step', label='raw')
n2, bins2, patches2 = plt.hist(x=np.array(e_grid)[full_muon_info['E_bin']], bins=np.array(e_grid), alpha=1.0, edgecolor='cyan', linestyle='dashed', histtype='step', label='After binning')
plt.yscale('log')
plt.xscale('log')
plt.ylabel(r"Count")
plt.xlabel(r"E (GeV)")
#plt.xlim(0, 10)
#plt.ylim(0, 1000)
plt.title(r"Energy Distribution")
plt.legend(loc='best',frameon=False,numpoints=1,fontsize='small')
plt.tight_layout()
plt.grid(True)

plt.savefig('/data/p-one/dghuman/simulation/MCEq_tables/test_plots/E_dist.png', dpi=300)


'''

plt.clf()

plt.figure(figsize=(4.2, 3))

mainbin = np.arange(-5, 85, 10)
n, bins, patches = plt.hist(x=np.array(full_muon_info['theta']), bins=mainbin, label=r"MCEq", alpha=1.0, edgecolor='black', linestyle='solid', histtype='step', weights=np.array(full_muon_info['time_length']))
n2, bins2, patches2 = plt.hist(x=np.array(full_muon_info['theta']), bins=mainbin, label=r"Mute", alpha=1.0, edgecolor='cyan', linestyle='--', histtype='step', weights=np.array(full_muon_info['mute_time_length']))
#plt.yscale('log')
plt.ylabel(r"Livetime (s)")
plt.xlabel(r"Zenith ($^{o}$)")
plt.yscale('log')
plt.title(r"Livetime Contributions")
plt.legend(loc='upper left',frameon=False,numpoints=1,fontsize='medium')
plt.tight_layout()
plt.grid(True)

plt.savefig('/data/p-one/dghuman/simulation/MCEq_tables/test_plots/live_time.png', dpi=300)
'''

plt.clf()

livetimes = np.array(muon_weights['mintimes'])

plt.figure(figsize=(4.2, 3))
mainbin = np.arange(-5, 85, 10)
width = 1*(mainbin[1] - mainbin[0])
#n, bins, patches = plt.hist(x=np.array(full_muon_info['theta']), bins=mainbin, label=r"MCEq", alpha=1.0, edgecolor='black', linestyle='solid', histtype='step', weights=np.array(full_muon_info['time_length']))
plt.bar(mainbin, livetimes, align='edge', width=width)
plt.yscale('log')
plt.ylabel(r"Livetime (s)")
plt.xlabel(r"Zenith ($^{o}$)")
#plt.yscale('log')
plt.title(r"Livetime Contributions")
#plt.legend(loc='upper left',frameon=False,numpoints=1,fontsize='medium')
plt.tight_layout()
plt.grid(True)

plt.savefig('/data/p-one/dghuman/simulation/MCEq_tables/test_plots/live_time.png', dpi=300)

plt.clf()

weights = np.array(muon_weights['angle_weight'])

plt.figure(figsize=(4.2, 3))
width = 1*(mainbin[1] - mainbin[0])
#n, bins, patches = plt.hist(x=np.array(full_muon_info['theta']), bins=mainbin, label=r"MCEq", alpha=1.0, edgecolor='black', linestyle='solid', histtype='step', weights=np.array(full_muon_info['time_length']))
plt.bar(mainbin, weights, align='edge', width=width)
plt.yscale('log')
plt.ylabel(r"Weight")
plt.xlabel(r"Zenith ($^{o}$)")
#plt.yscale('log')
plt.title(r"Weights per Zenith")
#plt.legend(loc='upper left',frameon=False,numpoints=1,fontsize='medium')
plt.tight_layout()
plt.grid(True)

plt.savefig('/data/p-one/dghuman/simulation/MCEq_tables/test_plots/weights_zenith.png', dpi=300)

plt.clf()

# plot livetime against E for each theta
for theta_i in range(len(weights)):
    values = np.zeros(len(e_grid))
    values[muon_weights[theta_i]['E_index']] = muon_weights[theta_i]['livetime']
    plt.loglog(e_grid, values, color=str(colors[theta_i].hex), ls='', marker='.', label=r'$\theta =$ ' + str(theta_i*10))

plt.xlabel(r"$E_{\mu}$ [GeV]")
plt.ylabel(r"Livetime (s)")
plt.title(r"LiveTime vs Energy")
plt.legend(loc='upper right',frameon=False,numpoints=1,fontsize='small')
plt.tight_layout()
plt.xlim(10,1e10)
#plt.ylim(1e-17,10)
plt.grid(True)
    
plt.savefig('/data/p-one/dghuman/simulation/MCEq_tables/test_plots/E_vs_livetime.png', dpi=300)

plt.clf()

plt.figure(figsize=(4.2, 3))
width = 1*(mainbin[1] - mainbin[0])
n, bins, patches = plt.hist(x=np.array(full_muon_info['zeros']), bins=e_grid, label=r"zeros", alpha=1.0, edgecolor='black', linestyle='solid', histtype='step')
plt.yscale('log')
plt.ylabel(r"Count")
plt.xlabel(r"$E_{\mu}$ [GeV]")
plt.xscale('log')
plt.title(r"Energy Count (0)")
#plt.legend(loc='upper left',frameon=False,numpoints=1,fontsize='medium')
plt.tight_layout()
plt.grid(True)

plt.savefig('/data/p-one/dghuman/simulation/MCEq_tables/test_plots/zeros_E.png', dpi=300)

plt.clf()

# plot weight against E for each theta
for theta_i in range(len(weights)):
    values = np.zeros(len(e_grid))
    values[muon_weights[theta_i]['E_index']] = muon_weights[theta_i]['e_weight']
    plt.loglog(e_grid, values, color=str(colors[theta_i].hex), ls='', marker='.', label=r'$\theta =$ ' + str(theta_i*10))

plt.xlabel(r"$E_{\mu}$ [GeV]")
plt.ylabel(r"Weight (Energy*Angle)")
plt.title(r"Weight vs Energy")
plt.legend(loc='upper right',frameon=False,numpoints=1,fontsize='small')
plt.tight_layout()
plt.xlim(10,1e10)
#plt.ylim(1e-17,10)
plt.grid(True)
    
plt.savefig('/data/p-one/dghuman/simulation/MCEq_tables/test_plots/E_vs_weight.png', dpi=300)

plt.clf()

n, bins, patches = plt.hist(x=np.array(muon_dict['Energy']), bins=e_grid, label=r"Untriggered", alpha=1.0, edgecolor='black', linestyle='solid', histtype='step')
n2, bins2, patches2 = plt.hist(x=np.array(muon_dict_trig['Energy']), bins=e_grid, label=r"triggered", alpha=1.0, edgecolor='red', linestyle='dashed', histtype='step')
plt.yscale('log')
plt.ylabel(r"Count")
plt.xlabel(r"Energy")
plt.xscale('log')
plt.title(r"Energy distribution")
plt.legend(loc='best',frameon=False,numpoints=1,fontsize='medium')
plt.tight_layout()
plt.grid(True)

plt.savefig('/data/p-one/dghuman/simulation/MCEq_tables/test_plots/E-dist.png', dpi=300)

plt.clf()
zbins = np.linspace(0, 1, 100)
n, bins, patches = plt.hist(x=np.cos(np.array(muon_dict['theta'])*np.pi/180), bins=zbins, label=r"Untriggered", alpha=1.0, edgecolor='black', linestyle='solid', histtype='step')
n2, bins2, patches2 = plt.hist(x=np.cos(np.array(muon_dict_trig['theta'])*np.pi/180), bins=zbins, label=r"triggered", alpha=1.0, edgecolor='red', linestyle='dashed', histtype='step')
plt.yscale('log')
plt.ylabel(r"Count")
plt.xlabel(r"cos(Zenith)")
plt.title(r"Zenith distribution")
plt.legend(loc='best',frameon=False,numpoints=1,fontsize='medium')
plt.tight_layout()
plt.grid(True)

plt.savefig('/data/p-one/dghuman/simulation/MCEq_tables/test_plots/Zenith-dist.png', dpi=300)

plt.clf()
abins = np.linspace(0, 360, 360)
n, bins, patches = plt.hist(x=np.array(muon_dict['phi']), bins=abins, label=r"Untriggered", alpha=1.0, edgecolor='black', linestyle='solid', histtype='step')
n2, bins2, patches2 = plt.hist(x=np.array(muon_dict_trig['phi']), bins=abins, label=r"triggered", alpha=1.0, edgecolor='red', linestyle='dashed', histtype='step')
plt.yscale('log')
plt.ylabel(r"Count")
plt.xlabel(r"Azimuth")
plt.title(r"Azimuth distribution")
plt.legend(loc='best',frameon=False,numpoints=1,fontsize='medium')
plt.tight_layout()
plt.grid(True)

plt.savefig('/data/p-one/dghuman/simulation/MCEq_tables/test_plots/Azimuth-dist.png', dpi=300)

