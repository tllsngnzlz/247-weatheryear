
"""
Function for plotting a nodal balance considering pypsa-eur-sec component structure
if you can ;)

"""

import pandas as pd
import matplotlib.colors as mc
import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox
import yaml
import pypsa

from pypsa.descriptors import Dict
from pypsa.plot import add_legend_patches


datacenters = {"IE5 0": "dublin", "DK1 0": "frederica"} #snakemake.config['ci']['datacenters']
locations = list(datacenters.keys())
names = list(datacenters.values())
    
rename = {}

for name in names:
    temp = {
    f'{name} H2 Electrolysis': 'hydrogen storage',
    f'{name} H2 Fuel Cell': 'hydrogen storage',
    f'{name} battery charger': 'battery storage',
    f'{name} battery discharger': 'battery storage',
    f'{name} export': 'grid',
    f'{name} import': 'grid',
    f'{name} onwind': 'wind',	
    f'{name} solar': 'solar',
    f'{name} load': 'load',
    f'{name} adv_geothermal': "clean dispatchable",
    f'{name} allam_ccs': "NG-Allam"}
    rename = rename | temp

links = {
    'vcc12': 'spatial shifting',	
    'vcc21': 'spatial shifting'}

rename = rename | links
 

preferred_order_balances = pd.Index([
    "clean dispatchable",
    "NG-Allam",
    "solar",
    "wind",
    "load",
    "battery storage",
    "hydrogen storage",
    'grid',
    'spatial shifting',
    ])


def retrieve_nb(n, node):
    '''
    Retrieve nodal energy balance per hour
    This simple function works only for the Data center nodes: 
        -> lines and links are bidirectional AND their subsets are exclusive.
        -> links include fossil gens
    NB {-1} multiplier is a nodal balance sign
    '''

    components=['Generator', 'Load', 'StorageUnit', 'Store', 'Link', 'Line']
    nodal_balance = pd.DataFrame(index=n.snapshots)
    
    for i in components:
        if i == 'Generator':
            node_generators = n.generators.query('bus==@node').index
            nodal_balance = nodal_balance.join(n.generators_t.p[node_generators])
        if i == 'Load':
            node_loads = n.loads.query('bus==@node').index
            nodal_balance = nodal_balance.join(-1*n.loads_t.p_set[node_loads])
        if i == 'Link':
            node_export_links = n.links.query('bus0==@node').index
            node_import_links = n.links.query('bus1==@node').index
            nodal_balance = nodal_balance.join(-1*n.links_t.p0[node_export_links])
            nodal_balance = nodal_balance.join(-1*n.links_t.p1[node_import_links])
            ##################
        if i == 'StorageUnit':
            #node_storage_units = n.storage_units.query('bus==@node').index
            #nodal_balance = nodal_balance.join(n.storage_units_t.p_dispatch[node_storage_units])
            #nodal_balance = nodal_balance.join(n.storage_units_t.p_store[node_storage_units]) 
            continue   
        if i == 'Line':
            continue
        if i == 'Store':
            continue

        nodal_balance = nodal_balance.rename(columns=rename).groupby(level=0, axis=1).sum()

    return nodal_balance


def plot_nb(n, node, 
            start='2013-03-01 00:00:00', 
            stop='2013-03-08 00:00:00'):

    fig, ax = plt.subplots()
    fig.set_size_inches((6,4.5))

    tech_colors = snakemake.config['tech_colors']
    df = retrieve_nb(n, node)

    #format time & set a range to display
    ldf = df.loc[start:stop,:]
    duration = (pd.to_datetime(stop)-pd.to_datetime(start)).days
    ldf.index = pd.to_datetime(ldf.index, format='%%Y-%m-%d %H:%M:%S').strftime('%m.%d %H:%M')

    #get colors
    for item in ldf.columns:
        if item not in tech_colors:
            print("Warning!",item,"not in config/tech_colors")

    #set logical order
    new_index = preferred_order_balances.intersection(ldf.columns).append(ldf.columns.difference(preferred_order_balances))
    ldf = ldf.loc[:,new_index]
    
    ldf.plot(kind="bar",stacked=True,
            color=tech_colors, 
            ax=ax, width=1, 
            edgecolor = "black", linewidth=0.05
            )

    #visually ensure net energy balance at the node
    net_balance=ldf.sum(axis=1)
    x = 0
    for i in range(len(net_balance)):
        ax.scatter(x = x, y = net_balance[i], color='black', marker="_")
        x += 1

    plt.xticks(rotation=90)
    ax.grid(alpha=0.3)
    ax.set_axisbelow(True)
    #ax.set_xlabel("Hours")
    ax.set(xlabel=None)
    ax.xaxis.set_major_locator(plt.MaxNLocator((duration*2+1)))

    ax.set_ylabel("Nodal balance [MW*h/h]")
    #ax.legend(loc="upper left", ncol = 3, prop={"size":8})

    add_legend_patches(
        ax,
        colors = [tech_colors[c] for c in ldf.columns],
        labels= ldf.columns,
        legend_kw= dict(bbox_to_anchor=(1, 1), loc = "upper left", frameon = False,
        ))

    fig.tight_layout()
    _start = ldf.index[0].split(' ')[0]
    _stop = ldf.index[-1].split(' ')[0]
    fig.savefig(snakemake.output.plot.replace("balance.pdf",f"balance_{node}_{_start}-{_stop}.pdf"), 
                bbox_inches = Bbox([[0,0],[7.7,4.5]])
                )


if __name__ == "__main__":
    # Detect running outside of snakemake and mock snakemake for testing
    if 'snakemake' not in globals():
        # from _helpers import mock_snakemake
        # snakemake = mock_snakemake('plot_balances', 
        #             year='2025', zone='IEDK', palette='p1', policy="cfe100", flexibility='0')   
        
        snakemake = Dict()
        with open(f"../config.yaml",'r') as f:
            snakemake.config = yaml.safe_load(f)
        snakemake.input = Dict()
        snakemake.output = Dict()

    folder="../results/test-IEDK-all"
    #snakemake.input.data = f"{folder}/networks/{scenario}/ref.csv"
    snakemake.output.plot = f"{folder}/plots/2025/IEDK/p1/cfe100/balance.pdf"
    solved_network = f"{folder}/networks/2025/IEDK/p1/cfe100/0.nc"
    n = pypsa.Network(solved_network)

    #Wildcards & Settings
    # year = snakemake.wildcards.year
    # area = snakemake.config['area']
    # zone = snakemake.wildcards.zone
    # tech_palette = snakemake.wildcards.palette
    # policy = snakemake.wildcards.policy[:3]
    # flexibility = snakemake.wildcards.flexibility
    # n = pypsa.Network(snakemake.input.network)

    datacenters = snakemake.config['ci']['datacenters']
    locations = list(datacenters.keys())
    names = list(datacenters.values())

    node = 'frederica'

for node in names:
    plot_nb(n, node, '2013-01-01 00:00:00', '2013-01-08 00:00:00')
    plot_nb(n, node, '2013-03-01 00:00:00', '2013-03-08 00:00:00')
    plot_nb(n, node, '2013-05-01 00:00:00', '2013-05-08 00:00:00')

