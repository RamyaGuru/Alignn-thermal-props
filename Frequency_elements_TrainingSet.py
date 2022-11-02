# Use /home/ftavazza/Sofware/README_bokeh to activate whatever is needed to run Bokeh
# Elements never used are grayed-out in the perioedic table

from jarvis.db.jsonutils import loadjson
from jarvis.core.atoms import Atoms
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.metrics import mean_absolute_error
from scipy.stats import kde
from scipy.stats import gaussian_kde

from bokeh.io import export_svgs
from bokeh.io import export_png

from bokeh.models import (
    ColumnDataSource,
    LinearColorMapper,
    LogColorMapper,
    ColorBar,
    FixedTicker,
    BasicTicker,
)
from bokeh.plotting import figure
from bokeh.io import output_file, save, show
from bokeh.sampledata.periodic_table import elements
from bokeh.transform import dodge
from matplotlib.colors import Normalize, LogNorm, to_hex
from matplotlib.cm import plasma, viridis, magma, inferno, ScalarMappable
from bokeh.core.properties import value, field

import numpy as np

run = "run21"

f_in='../../{}/PhononData-Freq-300_1000_20.json'.format(run)

dos_data=loadjson(f_in)

freq = np.linspace(-300, 1000, len(dos_data[0]['pdos_elast']))
        
def label_dynamic_stability(freq, dos_dict):
    stable_list = []
    unstable_list = []
    zero_indx = np.where(freq > 0)[0][0] - 1
    neg_freq = np.linspace(-300, 0, zero_indx)
    for d in dos_dict:
        intdos_neg_target = np.trapz(neg_freq, d['pdos_elast'][:zero_indx])
        intdos_target = np.trapz(freq, d['pdos_elast'])
        if intdos_neg_target / intdos_target > 0.1:
            unstable_list.append(d)
        else:
            stable_list.append(d)
    return stable_list, unstable_list  
        

stable_list, unstable_list = label_dynamic_stability(freq, dos_data)

data = unstable_list

periodic=['H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe','Cs','Ba','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn','Nh','Fl','Mc','Lv','Ts','Og','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr']

#f_name2 = "PeriodicTable_frequency-Training_ingr.png"
f_name2 = "PeriodicTable_frequency-Training-no-O_unstable.png"
f_name3 = "PeriodicTable_CompoundType-Training.png"
f_name4 = "PeriodicTable_Natoms-Training.png"
f_name5 = "PeriodicTable_VolpAt-Training.png"


def main():
    print("Number of materials in dataset=",len(data))

    # el_count = projection on each element of how often that element is part of a material in the dataset
    # el_compound = projection on each element of its average compound type (how often it is part of a binary, ternary etc.)
    # el_nat = projection on each element of the average numer of atoms in the primitive cell of materials it's part of
    # el_vat = projection on each element of the average volume per atom of materials it's part of
    el_count={}
    el_compound={}
    el_nat={}
    el_vat={}
    for el in periodic:
        el_count[el]=0
        el_compound[el]=0
        el_nat[el]=0
        el_vat[el]=0

    for i in data:
        box=i['atoms']['lattice_mat']
        coords=i['atoms']['coords']
        pp = i['atoms']['elements']
        mat = Atoms(lattice_mat=box, coords=coords, elements=pp)
        v_pa=round(float(mat.volume) / float(mat.num_atoms))
        nat = len(pp)
        elems=set(pp)
        n_elem=len(elems) 
        for ik in elems:   
            el_count[ik] = el_count[ik] + 1
            el_compound[ik] = el_compound[ik] + n_elem
            el_nat[ik] = el_nat[ik] + nat
            el_vat[ik] = el_vat[ik] + v_pa
            

    for el in periodic:
        if (el_count[el] > 0):
           el_compound[el] = el_compound[el]/el_count[el]
           el_nat[el] = el_nat[el]/el_count[el]
           el_vat[el] = el_vat[el]/el_count[el]
           

    periodic_new=[]
    el_frequency=[]
    el_compound_type=[]
    el_nat1=[]
    el_vat1=[]
    for key in el_count.keys():
        #if key != 'O':
        if el_count[key] > 0:
           periodic_new.append(key)
           el_frequency.append(el_count[key])
           el_compound_type.append(el_compound[key])
           el_nat1.append(el_nat[key])
           el_vat1.append(el_vat[key])
           
    np.savez("unstable_elem_freq.npz", periodic_new, el_frequency)       

    # el_frequency = how often a specific element appeared in materials for which we 
    #                computed the VALIDATION (phonon) spectrum
    x = plot_ptable_trend(periodic_new, el_frequency)
    export_png(x, filename=f_name2)

    # x = plot_ptable_trend(periodic_new, el_compound_type)
    # export_png(x, filename=f_name3)

    # x = plot_ptable_trend(periodic_new, el_nat1)
    # export_png(x, filename=f_name4)

    # x = plot_ptable_trend(periodic_new, el_vat1)
    # export_png(x, filename=f_name5)


def plot_ptable_trend(
    data_elements=["Rb", "S", "Se"],
    data_list=[10, 20, 30],
    input_file=None,
    output_html="ptable.html",
    bokeh_palette="Inferno256",
    cmap=inferno,
    log_scale=0,
    width=1050,
    alpha=0.65,
    cbar_height=520,
    cbar_font="14pt",
    save_plot=True,
):
    """
    Generate periodic table chemical trends.
    Either provide a file or list of data_elements, &data_list.
    Note that Bokeh already provided a periodic table.
    This module will take your data to color code them.
    See an example: https://www.nature.com/articles/s41598-019-45028-y
    Fig. 3
    Forked from https://github.com/arosen93/ptable_trends
    """
    output_file(output_html)
    # Define number of and groups
    period_label = ["1", "2", "3", "4", "5", "6", "7"]
    group_range = [str(x) for x in range(1, 19)]
    if input_file is not None:
        data_elements = []
        data_list = []

        f = open(input_file, "r")
        lines = f.read().splitlines()
        f.close()
        for i in lines:
            data_elements.append(i.split()[0])
            data_list.append(i.split()[1])

    print("Data_elements:")
    print(data_elements)
    print("Data_list:")
    print(data_list)
    
    
    data = [float(i) for i in data_list]

    if len(data) != len(data_elements):
        raise ValueError("Unequal number of atomic elements and data points")

    # lanthanides = [x.lower() for x in elements["symbol"][56:70].tolist()]
    # actinides = [x.lower() for x in elements["symbol"][88:102].tolist()]
    period_label.append("blank")
    period_label.append("La")
    period_label.append("Ac")

    count = 0
    for i in range(56, 70):
        elements.period[i] = "La"
        elements.group[i] = str(count + 4)
        count += 1

    count = 0
    for i in range(88, 102):
        elements.period[i] = "Ac"
        elements.group[i] = str(count + 4)
        count += 1

    # Define matplotlib and bokeh color map
    if log_scale == 0:
        color_mapper = LinearColorMapper(
            palette=bokeh_palette, low=min(data), high=max(data)
        )
        norm = Normalize(vmin=min(data), vmax=max(data))
    elif log_scale == 1:
        for i in range(len(data)):
            if data[i] < 0:
                raise ValueError(
                    "Entry for element "
                    + data_elements[i]
                    + " is negative but"
                    " log-scale is selected"
                )
        color_mapper = LogColorMapper(
            palette=bokeh_palette, low=min(data), high=max(data)
        )
        norm = LogNorm(vmin=min(data), vmax=max(data))
    color_scale = ScalarMappable(norm=norm, cmap=cmap).to_rgba(
        data, alpha=None
    )

    # Define color for blank entries
    blank_color = "#c4c4c4"
    color_list = []
    for i in range(len(elements)):
        color_list.append(blank_color)

    # Compare elements in dataset with elements in periodic table
    for i in range(len(data)):
        element_entry = elements.symbol[
            elements.symbol.str.lower() == data_elements[i].lower()
        ]
        if not element_entry.empty:
            element_index = element_entry.index[0]
        else:
            print("WARNING: Invalid chemical symbol: " + data_elements[i])
        if color_list[element_index] != blank_color:
            print("WARNING: Multiple entries for element " + data_elements[i])
        color_list[element_index] = to_hex(color_scale[i])

    # Define figure properties for visualizing data
    source = ColumnDataSource(
        data=dict(
            group=[str(x) for x in elements["group"]],
            period=[str(y) for y in elements["period"]],
            sym=elements["symbol"],
            atomic_number=elements["atomic number"],
            type_color=color_list,
        )
    )
    # Plot the periodic table
    p = figure(
        x_range=group_range, y_range=list(reversed(period_label)),\
            title="Counts in Training Set", title_location = "right", tools="save"
    )
    p.plot_width = width
    p.outline_line_color = None
    p.toolbar_location = "above"
    p.rect(
        "group",
        "period",
        0.9,
        0.9,
        source=source,
        alpha=alpha,
        color="type_color",
    )
    p.axis.visible = False
    text_props = {
        "source": source,
        "angle": 0,
        "color": "black",
        "text_align": "left",
        "text_baseline": "middle",
    }
    x = dodge("group", -0.4, range=p.x_range)
    y = dodge("period", 0.3, range=p.y_range)
    
    y2 = dodge("period", -0.05, range=p.y_range)
    p.text(
        x=x,
        y=y2,
        text="sym",
        text_font=value("helvetica"),
        text_font_style="bold",
        text_font_size="20pt",
        **text_props
    )
    p.text(x=x, y=y, text="atomic_number", text_font_size="9pt", **text_props)
    color_bar = ColorBar(
        color_mapper=color_mapper,
        ticker=BasicTicker(desired_num_ticks=10),
        border_line_color=None,
        label_standoff=6,
        major_label_text_font_size=cbar_font,
        location=(0, 0),
        orientation="vertical",
        scale_alpha=alpha,
        width=8,
    )
    p.title.align = "center"
    p.title.text_font_size = "20px"
    if cbar_height is not None:
        color_bar.height = cbar_height

    p.add_layout(color_bar, "right")
    p.grid.grid_line_color = None
    if save_plot:
        save(p)
    else:
        show(p)
    return p

main()
