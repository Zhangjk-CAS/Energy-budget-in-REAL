import matplotlib.pyplot as plt
import numpy as np


# adv = [(-3.47, -3.58, 1.16, -4.63, -2.42),
#        (-12.8, -12.6, -10.4, -7.5, -6.7),
#        (-10.21, -8.78, -4.49, 1.19, -2.77),
#        (-13.81, -12.78, -8.56, -8.87, -2.6)]

# pwork = [(135.75, 135.39, 129.94, 142.69, 124.81),
#          (101.62, 83.64, 79.43, 21.79, 16.7)]

# eddy =[(-20.72, -21.34, -19.16, -11.52, -9.66),
#        (16.63, 13.53, 6.77, 9.31, 1.61)]

adv = [(-3.58, -3.79, 1.08, -3.03, -1.45),
       (-14.21, -14.09, -11.50, -8.11, -6.78),
       (-8.71, -10.15, -4.39, 1.27, -2.73),
       (-12.90, -13.95, -8.64, -7.94, -2.81)]

pwork = [(136.27, 136.62, 130.19, 148.30, 124.81),
         (83.63, 101.82, 79.67, 16.69, 19.05)]

eddy =[(-21.86, -20.84, -19.25, -11.45, -9.72),
       (13.51, 16.65, 6.74, 9.30, 3.56)]

colors = ['blue', 'red', 'green', 'blue', 'green'] 
ln = ['Reg-EG','Reg-JG', 'Reg-JL', 'GLORYS', 'LICOM3']
width = 0.9
group_width = 7
fig = plt.figure(figsize=(7.5, 5))
gs = fig.add_gridspec(2, 2, height_ratios=[1, 1], width_ratios=[1, 1])

print('bar of energy advection')
ax1 = fig.add_subplot(gs[0, :])
for n in range(4):
    bars = ax1.bar(np.arange(5)+group_width*n, adv[n], color=colors, alpha=0.7, width=width, linestyle='-', edgecolor='black')
    ax1.bar_label(bars, labels=[f'{val:.3g}' for val in adv[n]], fontsize=6, padding=1)
    edge_styles = ['-', '-', '-', '--', '--']  
    for bar, style in zip(bars, edge_styles):
        bar.set_linestyle(style)  
# ax1.set_title('(a)Energy advection', loc='left', pad=0.2)

bars = ax1.bar(100, 1, color=colors[0], alpha=0.7, width=width, linestyle='-', edgecolor='black', label=ln[0])
bars = ax1.bar(100, 1, color=colors[1], alpha=0.7, width=width, linestyle='-', edgecolor='black', label=ln[1])
bars = ax1.bar(100, 1, color=colors[2], alpha=0.7, width=width, linestyle='-', edgecolor='black', label=ln[2])
bars = ax1.bar(100, 1, color=colors[3], alpha=0.7, width=width, linestyle='--', edgecolor='black', label=ln[3])
bars = ax1.bar(100, 1, color=colors[4], alpha=0.7, width=width, linestyle='--', edgecolor='black', label=ln[4])
ax1.legend(loc='upper center', ncol=5, frameon=False, fontsize=7)

print('bar of energy convertion loss')
ax2 = fig.add_subplot(gs[1, 0])
for n in range(2):  
    bars = ax2.bar(np.arange(5)+group_width*n, eddy[n], color=colors, alpha=0.7, width=width, linestyle='-', edgecolor='black')
    ax2.bar_label(bars, labels=[f'{val:.3g}' for val in eddy[n]], fontsize=6, padding=1)
    edge_styles = ['-', '-', '-', '--', '--']  
    for bar, style in zip(bars, edge_styles):
        bar.set_linestyle(style)  
# ax2.set_title('(b)Energy conversion loss', loc='left', pad=0.2)

print('bar of pressure work')
ax3 = fig.add_subplot(gs[1, 1])
for n in range(2):  
    bars = ax3.bar(np.arange(5)+group_width*n, pwork[n], color=colors, alpha=0.7, width=width, linestyle='-', edgecolor='black')
    ax3.bar_label(bars, labels=[f'{val:.3g}' for val in pwork[n]], fontsize=6, padding=1)
    edge_styles = ['-', '-', '-', '--', '--']  
    for bar, style in zip(bars, edge_styles):
        bar.set_linestyle(style)  
# ax3.set_title('(c)Pressure work', loc='left', pad=0.2)

ax1.set_xticks(np.arange(28)[::group_width] + width*2 + 0.5)  
ax1.set_xticklabels(['In MKE', 'In EKE', 'In MPE', 'In EPE'])
ax2.set_xticks(np.arange(14)[::group_width] + width*2 + 0.2)
ax2.set_xticklabels(['In KE', 'In PE'])
ax3.set_xticks(np.arange(14)[::group_width] + width*2 + 0.2)
ax3.set_xticklabels(['In MKE', 'In EKE'])

ax1.set_ylim(-25,25)
ax2.set_ylim(-25,25)
ax3.set_ylim(0,159)
ax1.set_xlim(-1,26)

# plt.xlabel('terms', fontsize=12)
ax1.set_ylabel('GW', fontsize=12)
ax2.set_ylabel('GW', fontsize=12)
# ax3.set_ylabel('GW', fontsize=12)
# plt.title('Basic Bar Chart', fontsize=14)

pnstr2 = f'energyBar'
pn = 'figures/' + pnstr2 + '.png'
plt.savefig(pn, dpi=600)
plt.close()