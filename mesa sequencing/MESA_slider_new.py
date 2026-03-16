#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
from glob import glob
from matplotlib.artist import Artist
from scipy.interpolate import interp1d

import mesa_reader as mr

from matplotlib.widgets import Button, Slider

fontsize=18
plt.rcParams.update({'font.size': fontsize})

def get_data(f):            # f is the filename to read
    h = mr.MesaData(f)
    mods = h.model_number
    log_Teffs = h.log_Teff
    ages = h.star_age
    log_Ls = h.log_L
    log_LHs = h.log_LH
    center_h1 = h.center_h1
    center_he4 = h.center_he4
    center_he3 = h.center_he3
    center_c12 = h.center_c12
    center_si28 = h.center_si28
    log_LHes = h.log_LHe
    star_mass = h.star_mass
    Z_init = 1. - center_h1[0] - center_he4[0] - center_he3[0]
    return mods, log_Teffs, ages, log_Ls, log_LHs, log_LHes, center_h1, center_he4, center_c12, center_si28, star_mass, Z_init, h

files1 = glob('*/*hist*')
files2 = glob('*hist*')
files = files1 + files2

if len(files) > 0:
    print('Found the following files:')
    print('  File no.      Name')
    for ii,fi in enumerate(files):
        out = f'     {ii}          {fi}'
        print(out)
else:
    print('No history files found.')
    exit(0)

cnum = input('\nEnter the number of the file name you wish to plot:  ')
num = int(cnum)
if num < 0 or num >= len(files):
    print('Incorrect number entered. Must be in the range 0 to ' + f'{len(files)-1}')
    exit(0)
    
f1 = files[num]
print('\nReading ' + f1 + '...\n')

mods, log_Teffs, ages, log_Ls, log_LHs, log_LHes, center_h1, center_he4, center_c12,  center_si28, star_mass, Z_init, h = get_data(f1)

# Parse track to find approximate EEPs
modH0=-1
modMS=-1
modHeburn=-1
modHe0=-1
for n in range(len(mods)-1,1,-1):
    if center_h1[n]<=1.e-08:
        modH0=n
    if abs(log_LHs[n] - log_Ls[n]) <= 0.0005:
        modMS=n
    if log_LHes[n] >= (log_Ls[n]-2.0):
        modHeburn=n
    if center_he4[n]<=1.e-08:
        modHe0=n
modmap = np.array([1,  100,  250,      500,   750,    1000])
modnums= np.array([1,modMS,modH0,modHeburn,modHe0,len(center_h1)-1])
mask = modnums > 0
modnums = modnums[mask]
modmap  = modmap[mask]
fint = interp1d(modmap,modnums)


# Create the figure and the point on the track that we will manipulate
fig, ax = plt.subplots(figsize=(10,7))
ax.plot(log_Teffs,log_Ls,'-') 
point, = ax.plot(log_Teffs[0],log_Ls[0],'o',ms=8) 
tlab=ax.text(0.05, 0.95, '')
ax.set_xlabel('$\log T_{eff}$')
ax.set_ylabel('$\log L_\star/L_{\odot}$')
ax.invert_xaxis()

# adjust the main plot to make room for the sliders, buttons, and variable values
plt.subplots_adjust(right=0.76, bottom=0.32)

# Make a horizontal slider to control the model number
axfreq = fig.add_axes([0.25, 0.07, 0.65, 0.03])
freq_slider = Slider(
    ax=axfreq,
    label='Model Number',
    valmin=1,
    valmax=int(max(modmap)),
    valinit=250,
)

# --- A/B/C/D marker state ---
marker_labels = ['A', 'B', 'C', 'D']
marker_colors = ['tab:red', 'tab:green', 'tab:purple', 'tab:orange']

# Each entry holds the plot objects for that marker (dot and text), or None if not placed
placed_markers = {label: {'dot': None, 'text': None} for label in marker_labels}

def place_marker(label):
    """Place or move the named marker to the current slider position and save a PDF."""
    global cnum
    color = marker_colors[marker_labels.index(label)]
    jval = int(freq_slider.val)
    i = int(fint(jval))
    x = log_Teffs[i]
    y = log_Ls[i]

    # Remove existing marker objects if already placed
    if placed_markers[label]['dot'] is not None:
        placed_markers[label]['dot'].remove()
        placed_markers[label]['text'].remove()

    # Draw marker dot and label
    dot, = ax.plot(x, y, 'o', ms=12, color=color, zorder=5,
                   markeredgecolor='black', markeredgewidth=0.8)
    txt = ax.text(x, y, f'{label} ', color=color, fontsize=14,
                  fontweight='bold', va='bottom', ha='right', zorder=6)

    placed_markers[label]['dot'] = dot
    placed_markers[label]['text'] = txt
    fig.canvas.draw_idle()

    # Save the figure to a PDF file
    pdf_filename = f'HRD-{cnum}-{label}.pdf'
    fig.savefig(pdf_filename, bbox_inches='tight')
    print(f'   Saved {pdf_filename}...')

# Create A/B/C/D buttons in a row above the slider
button_axes = []
buttons = []
btn_width = 0.08
btn_gap = 0.02
btn_start_x = 0.42
btn_y = 0.14
btn_height = 0.05

for idx, label in enumerate(marker_labels):
    bx = btn_start_x + idx * (btn_width + btn_gap)
    bax = fig.add_axes([bx, btn_y, btn_width, btn_height])
    btn = Button(bax, label, color=marker_colors[idx], hovercolor='0.85')
    btn.label.set_fontsize(14)
    btn.label.set_fontweight('bold')
    btn.on_clicked(lambda event, lbl=label: place_marker(lbl))
    button_axes.append(bax)
    buttons.append(btn)

# Instruction label to the left of the buttons
fig.text(0.15, btn_y + btn_height / 2, 'Place marker at current position:',
         fontsize=11, color='0.3', va='center')


def update(val):
    global tlab
    jval = int(freq_slider.val)
    i = np.array([int(fint(jval))])
    age     = ages[i]
    lg_Teff = log_Teffs[i]
    lg_L    = log_Ls[i]
    lg_H    = log_LHs[i]
    lg_He   = log_LHes[i]
    H_cntr  = center_h1[i]
    He_cntr = center_he4[i]
    C_cntr  = center_c12[i]
    Si_cntr = center_si28[i]
    mass    = star_mass[i]
    point.set_xdata(lg_Teff)
    point.set_ydata(lg_L)
    labm2 = r'$Z_{\rm init}=$' + '{0:.4f}'.format(Z_init)
    labm1 = r'$M_\star/M_\odot=$' + '{0:.4f}'.format(mass[0])
    lab0  = r'Age (yr) =' + '{0:.3g}'.format(age[0])
    lab1  = r'$\log T_{\rm eff}=$'
    lab2 = '{0:.3f}'.format(lg_Teff[0])
    lab3 = r'$\log L_\star=$'
    lab4 = '{0:.3f}'.format(lg_L[0])
    lab5 = r'$\log L_{He}=$'
    lab6 = '{0:.3f}'.format(lg_He[0])
    lab7 = r'$\log L_{H}=$'
    lab8 = '{0:.3f}'.format(lg_H[0])
    lab9 = r'$H_{\rm center}=$'
    lab10 = '{0:.4f}'.format(H_cntr[0])
    lab11 = r'$He_{\rm center}=$'
    lab12 = '{0:.4f}'.format(He_cntr[0])
    lab13 = r'$C_{\rm center}=$'
    lab14 = '{0:.4f}'.format(C_cntr[0])
    lab15 = r'$Si_{\rm center}=$'
    lab16 = '{0:.4f}'.format(Si_cntr[0])
    lab = labm1 + '\n' + labm2 + '\n' + lab0 + '\n' + lab1 + lab2 + '\n' + lab3 + lab4 + '\n' + lab7 + lab8 + '\n' + lab5 + lab6 + '\n' + lab9 + lab10 + '\n' + lab11 + lab12 + '\n' + lab13 + lab14 + '\n' + lab15 + lab16
    Artist.remove(tlab)
    tlab=ax.text(1.02, 0.97, lab,
        verticalalignment='top', horizontalalignment='left',
        transform=ax.transAxes,
        color='C3', fontsize=15)
    fig.canvas.draw_idle()

update(250)
freq_slider.on_changed(update)

plt.show()
