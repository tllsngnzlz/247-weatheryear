import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def concatenate_summary_files(summary_files):
    dfs = []
    for file_path in summary_files:
        year = file_path.split('/')[-2]  # Adjust this based on the actual path structure
        df = pd.read_csv(file_path, index_col=0)
        new_index = [x[:-5] if 'system_inv' in x else x for x in df.index]
        df.index = new_index
        df['year'] = year
        multi_index = pd.MultiIndex.from_arrays([df.index, df['year']], names=['variable', 'year'])
        df = df.drop(columns=['year']).set_index(multi_index)
        dfs.append(df)
    
    return pd.concat(dfs)

def plot_system_invest(df, output_path, tech_colors):
    exp_generators = ['offwind-ac', 
                    'offwind-dc', 
                    'onwind', 
                    'solar']
    exp_links = ['OCGT']
    exp_chargers = ['battery charger', 'H2 Electrolysis']
    exp_dischargers = ['battery discharger', 'H2 Fuel Cell']

    exp_generators = [item.replace(' ', '_') for item in exp_generators]
    exp_chargers = [item.replace(' ', '_') for item in exp_chargers]
    exp_dischargers = [item.replace(' ', '_') for item in exp_dischargers]
    
    preferred_order = pd.Index([
        "advanced dispatchable",
        "NG-Allam",
        'Gas OC',
        "offshore wind",
        "onshore wind",
        "solar",
        "battery",
        "hydrogen storage",
        "hydrogen electrolysis",
        "hydrogen fuel cell"])

    rename_system_simple = {
        'offwind-ac': 'offshore wind',
        'offwind-dc': 'offshore wind',
        'onwind': 'onshore wind',
        'solar': 'solar',
        'OCGT': 'Gas OC',
        'battery_discharger': 'battery',
        'H2_Fuel_Cell': 'hydrogen fuel cell',
        'H2_Electrolysis': 'hydrogen electrolysis'
    }

    rename_scen = {'ref': 'no\npolicy',
                    'res100': '100%\nRES',
                    'cfe80': '80%',
                    'cfe85':'85%',
                    'cfe90':'90%',
                    'cfe95':'95%',
                    'cfe98':'98%',
                    'cfe99':'99%',
                    'cfe100':'100%'
                    }
    
    
    # calculate the total system investment for each year in multiindex
    gens = df.loc[["system_inv_" + t for t in exp_generators]].rename({"system_inv_" + t : t for t in exp_generators})
    links = df.loc[["system_inv_" + t for t in exp_links]].rename({"system_inv_" + t : t for t in exp_links})
    dischargers = df.loc[["system_inv_" + t for t in exp_dischargers]].rename({"system_inv_" + t : t for t in exp_dischargers})
    chargers = df.loc[["system_inv_" + t for t in exp_chargers]].rename({"system_inv_" + t : t for t in exp_chargers})
    chargers = chargers.drop(['battery_charger']) # display only battery discharger capacity

    ldf = pd.concat([gens, links, dischargers, chargers])
    # ldf.drop(columns=['ref'], inplace=True)

    scenarios = ldf.columns
    years = ldf.index.get_level_values('year').unique()
    
    fig, axes = plt.subplots(len(scenarios), 1, figsize=(18, 6))
    # Ensure axes is always an iterable (list of Axes objects)
    if not isinstance(axes, np.ndarray):
        axes = [axes]

    for scenario_idx, scenario in enumerate(scenarios):
        # Sum values across all variables for each year and scenario
        ldf_sum = ldf.groupby(level='year').sum().sort_values(by=scenario, ascending=True)
        
        # Create a DataFrame to hold sorted variables for stacking
        sorted_df = pd.DataFrame()
        for year in ldf_sum.index:
            year_data = ldf.xs(year, level='year')[scenario]
            #set year data name to year
            year_data.name = year
            year_data_sorted = year_data.sort_values(ascending=True)
            sorted_df = pd.concat([sorted_df, year_data_sorted], axis=1)
        
        # rename index entries and sum index entries with the same name
        # sorted_df.index = sorted_df.index.map(rename_system_simple)
        sorted_df = sorted_df.groupby(rename_system_simple,).sum()
        # set the index to the preferred order
        sorted_df = sorted_df.reindex(preferred_order)
        # Plot years next to each other with columns touching
        ax = axes[scenario_idx]
        sorted_df.T.plot(kind='bar', stacked=True, color=tech_colors, width=1 ,ax=ax, legend=(scenario_idx==0), )
        ax.set_title(scenario)
        
        # Only add legend to the first subplot to avoid repetition
        if scenario_idx == 0:
            axes[scenario_idx].legend(title='technology', loc='upper right', bbox_to_anchor=(1.05, 1))
        
    
        # drop rows with all zeros or close to zero
        # ldf = ldf.loc[(ldf.sum(axis=1) > 1e-1)] # > 0.1 MW1

    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()   

def plot_system_cfe(df, output_path):
    
    ldf = df.loc["system_grid_cfe_wavg",:]

    # ldf.drop(columns=['ref'], inplace=True)

    scenarios = ldf.columns
    years = ldf.index.get_level_values('year').unique()
    
    fig2, axes2 = plt.subplots(len(scenarios), 1, figsize=(18, 6))
        # Ensure axes is always an iterable (list of Axes objects)
    if not isinstance(axes2, np.ndarray):
        axes2 = [axes2]

    for scenario_idx, scenario in enumerate(scenarios):
        # Sum values across all variables for each year and scenario
        sorted_df = ldf.groupby(level='year').sum().sort_values(by=scenario, ascending=True)[scenario]

        ax = axes2[scenario_idx]
        # Plot years next to each other with columns touching
        sorted_df.plot(kind='bar', color='purple',ax=ax, legend=(scenario_idx==0), )
        ax.set_title(scenario)
        
        # Only add legend to the first subplot to avoid repetition
        if scenario_idx == 0:
            axes2[scenario_idx].legend(title='technology', loc='upper right', bbox_to_anchor=(1.05, 1))
        
    
        # drop rows with all zeros or close to zero
        # ldf = ldf.loc[(ldf.sum(axis=1) > 1e-1)] # > 0.1 MW1

    plt.tight_layout()
    plt.savefig(output_path)


if __name__ == "__main__":
    # Detect running outside of snakemake and mock snakemake for testing
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('compare_all_weatheryears', policy="ref", palette='p1', zone='DE', year='2025', participation='10', weather_year='1980-1982')

    
    summary_files = snakemake.input.summary_files # this is the list of files from the expand function
    tech_colors = snakemake.config['tech_colors']
    
    concatenated_df = concatenate_summary_files(summary_files)
    
    plot_system_invest(concatenated_df, snakemake.output.plot_invest, tech_colors)
    plot_system_cfe(concatenated_df, snakemake.output.plot_cfe)
