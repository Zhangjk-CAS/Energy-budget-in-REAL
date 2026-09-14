import matplotlib.pyplot as plt
import numpy as np
import matplotlib
matplotlib.use('Agg')

n = 9
x = [1,2,3,4, 6,7,8,9, 11,12,13,14]
# means = np.array([
#     [ 1.75, 1.75, 1.64, 1.48,   0.93, 0.95, 0.86, 0.81,   2.74, 2.71, 2.64, 2.13],   # Taiwan Strait
#     [-15.62, -15.60, -13.43, -17.42,    -18.69, -17.92, -13.47, -17.07,   -12.95, -13.92, -13.88, -16.36],  # Luzon Strait
#     [ 0.69, 2.78, 3.22, 2.78,   1.55, 3.65, 4.89, 3.65,   0.76, 2.26, 2.35, 2.26],   # Mindoro Strait
#     # [-1.08, -1.07, -1.39,  -3.56, -3.65, -3.80,   1.48,  1.62,  1.29]    # Karimata Strait
#     [1.08, 1.07, 1.39, 1.35,   3.56, 3.65, 3.80, 2.97,   -1.48, -1.62, -1.29, -0.49]    # Karimata Strait
# ])

# stds = np.array([
#     [1.23, 1.29, 1.26, 0.86,   1.14, 1.26, 1.27, 0.90,   0.68, 0.70, 0.69, 0.39],   # Taiwan Strait
#     [6.51, 5.84, 5.34, 5.94,   5.43, 6.97, 5.53, 6.07,   7.21, 4.61, 4.40, 5.17],   # Luzon Strait
#     [2.60, 1.53, 2.22, 1.53,   2.05, 1.23, 1.86, 1.23,   2.31, 1.33, 1.76, 1.32],   # Mindoro Strait
#     [2.17, 2.25, 2.17, 1.51,   0.94, 0.91, 0.92, 0.63,   0.71, 0.66, 0.75, 0.72]    # Karimata Strait
# ])

# new luzon strait 
means = np.array([
    [ 1.75, 1.75, 1.64, 1.48,   0.93, 0.95, 0.86, 0.81,   2.74, 2.71, 2.64, 2.13],   # Taiwan Strait
    [-5.30, -4.91, -6.71, -7.86,    -7.70, -7.70, -10.27, -10.43,   -4.88, -4.07, -6.39, -6.88],  # Luzon Strait
    [ 0.69, 2.78, 3.22, 2.78,   1.55, 3.65, 4.89, 3.65,   0.76, 2.26, 2.35, 2.26],   # Mindoro Strait
    # [-1.08, -1.07, -1.39,  -3.56, -3.65, -3.80,   1.48,  1.62,  1.29]    # Karimata Strait
    [1.08, 1.07, 1.39, 1.35,   3.56, 3.65, 3.80, 2.97,   -1.48, -1.62, -1.29, -0.49]    # Karimata Strait
])

stds = np.array([
    [1.23, 1.29, 1.26, 0.86,   1.14, 1.26, 1.27, 0.90,   0.68, 0.70, 0.69, 0.39],   # Taiwan Strait
    [5.03, 4.59, 5.32, 4.11,   4.86, 5.02, 5.87, 3.08,   4.01, 3.79, 3.86, 3.21],   # Luzon Strait
    [2.60, 1.53, 2.22, 1.53,   2.05, 1.23, 1.86, 1.23,   2.31, 1.33, 1.76, 1.32],   # Mindoro Strait
    [2.17, 2.25, 2.17, 1.51,   0.94, 0.91, 0.92, 0.63,   0.71, 0.66, 0.75, 0.72]    # Karimata Strait
])


# labels = ['E-G (Total)', 'J-G (Total)', 'J-L (Total)',
#           'E-G (Win)', 'J-G (Win)', 'J-L (Win)',
#           'E-G (Sum)', 'J-G (Sum)', 'J-L (Sum)']
x_labels = [2.5,7.5,12.5]
labels = ['Annual', 'Winter', 'Summer']

straits = ['a) Taiwan Strait', 'b) Luzon Strait', 'c) Mindoro Strait', 'd) Karimata Strait']

fig, axes = plt.subplots(2, 2, figsize=(7, 14/3), sharex=True)

i = 0
cstr = ['blue', 'red', 'green', 'royalblue', 'tomato']
lnstr = ['o', 'o', 'o', '<', '<']
labelstr = ['Reg-EG','Reg-JG','Reg-JL','GLORYS']
for ii in range(2):
    for ij in range(2):
        for j in [0,4,8]:
            ax = axes[ii,ij]
            for nn in range(j, j+4):
                if (ii+ij+j) == 0 :
                    ax.errorbar(x[nn], means[i,nn], yerr=stds[i,nn], 
                        fmt=lnstr[nn-j], color=cstr[nn-j], label=labelstr[nn-j], capsize=4, elinewidth=1.2, markersize=5)
                else:
                    ax.errorbar(x[nn], means[i,nn], yerr=stds[i,nn], 
                        fmt=lnstr[nn-j], color=cstr[nn-j], capsize=4, elinewidth=1.2, markersize=5)
        ax.set_title(straits[i], fontsize=11, loc='left')
        # ax.axhline(0, color='k', linewidth=0.8, linestyle='--', alpha=0.6)
        ax.grid(True, linestyle='--', alpha=0.4)
        ax.legend(loc='upper left',ncol=2,frameon=False)
        if ij == 0:
            ax.set_ylabel('Transport (Sv)', fontsize=10)
        i = i + 1 

axes[1,0].set_xticks(x_labels)
axes[1,0].set_xticklabels(labels, ha='center', fontsize=11)
axes[1,1].set_xticks(x_labels)
axes[1,1].set_xticklabels(labels, ha='center', fontsize=11)

axes[0,0].set_ylim(-10,10)
axes[1,0].set_ylim(-10,10)
axes[0,1].set_ylim(-20, 5)
axes[1,1].set_ylim(-10,10)

plt.tight_layout()
pn = "figures/transport-errorbar-force0519.png"
plt.savefig(pn, dpi=600)
print(pn)