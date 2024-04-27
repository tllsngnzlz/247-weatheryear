import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import numpy as np
from solve_network import geoscope

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

def plot_objective(df, output_path):
    
    ldf = df.loc["objective",:]

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

def plot_rldc(files):
    
    policies = snakemake.config["scenario"]["policy"]
    name = snakemake.config['ci']['name']
    node = geoscope(snakemake.config['scenario']['zone'][0], snakemake.config['area'])['node']

    fig, axes = plt.subplots(len(policies), 2, figsize=(18, 15))
    axes = np.atleast_2d(axes)

    dfs = []
    for policy_idx, policy in enumerate(policies):
        selected_files = [file for file in files if policy in file]
        for selected_file in selected_files:
            year = selected_file.split('/')[-2]
            df = pd.read_csv(selected_file, index_col=0)
            df[node].plot(ax=axes[policy_idx, 0], title=policy)
            df[f"{name} load"].plot(ax=axes[policy_idx, 1], title=policy)
            df['year'] = year
            df['policy'] = policy
            multi_index = pd.MultiIndex.from_arrays([df.index, df['year'], df['policy']], names=['variable', 'year', 'policy'])
            df = df.drop(columns=['year', 'policy']).set_index(multi_index)
            dfs.append(df)
    dfs = pd.concat(dfs)
    dfs.to_csv(snakemake.output.rldc_csv)
    agg_df = dfs.groupby(level=['policy', 'year']).mean().assign(variable='rldc_sum')
    agg_df = agg_df.rename(columns={node : 'value'})
    agg_df = agg_df.set_index('variable', append=True).reorder_levels(['variable', 'year', 'policy'])
    global pick_df 
    pick_df = pick_df.melt(ignore_index=False, var_name='policy', value_name='value')
    pick_df = pick_df.set_index('policy', append=True)
    pick_df = pd.concat([pick_df,agg_df.loc[:,'value'].to_frame()],axis=0)
    plt.tight_layout()
    plt.savefig(snakemake.output.plot_rldc)

def pick_wy_table(df):

    # Step 1: Split the DataFrame based on 'policy'
    unique_policies = df.index.get_level_values('policy').unique()
    split_dfs = {policy: df.xs(policy, level='policy') for policy in unique_policies}

    # Step 2 & 3: Unstack 'variable' level and filter columns for each split DataFrame
    unstacked_filtered_dfs = {}
    for policy, split_df in split_dfs.items():
        # Unstack the 'variable' level
        unstacked_df = split_df.reset_index()
        unstacked_df = unstacked_df.pivot(index='year', columns='variable', values='value')
        # Filter columns
        columns_of_interest = {'objective': 0.25 , 'system_grid_cfe_wavg':0.25, 'emissions':0.25, 'rldc_sum':0.25}
        # Filtering and handling cases where the columns might not exist in the unstacked DataFrame
        filtered_df = unstacked_df.loc[:, unstacked_df.columns.isin(columns_of_interest.keys())]
        # Normalizing the filtered DataFrame and calculating the weighted score
        rank_df=filtered_df.copy()
        rank_df['system_grid_cfe_wavg'] = 1-rank_df.loc[:,'system_grid_cfe_wavg'] # flip the cfe to get the right ranking behaviour
        rank_df = (rank_df - rank_df.min()) / (rank_df.max() - rank_df.min())
        rank_df['weighted_score'] = rank_df.mul(columns_of_interest, axis=1).sum(axis=1)
        rank_df.columns = [col + "_weight" for col in rank_df.columns]
        # Storing the filtered, unstacked DataFrame with weights in a dictionary keyed by 'policy'
        filtered_df = pd.concat([filtered_df, rank_df], axis=1).sort_values(by='weighted_score_weight', ascending=False)
        unstacked_filtered_dfs[policy] = filtered_df
    unstacked_filtered_dfs = pd.concat(unstacked_filtered_dfs)
    unstacked_filtered_dfs.to_csv(snakemake.output.pick_csv)
    return unstacked_filtered_dfs   

def plot_wy(df):
    df.reset_index(inplace=True)
    # Split the dataframe based on the unique policies in column 0
    groups = df.groupby(df.columns[0])

    # Prepare the colormap
    total_years = len(df['year'].unique())
    colormap = cm.get_cmap('plasma', total_years)  # Get the 'plasma' colormap
    normalize = mcolors.Normalize(vmin=0, vmax=total_years-1)

    line_styles = ['-', '--', '-.', ':']
    markers = ['x'] #, 's', '^', 'D', '*', 'p', 'h']

    fig, ax = plt.subplots(groups.ngroups,1,figsize=(12, 6*groups.ngroups), sharex=True)
    if not isinstance(ax, np.ndarray):
        ax = [ax]

    # Iterate over each group
    for j, (group_name, group_df) in enumerate(groups):
        # Identify columns to plot, excluding the last "result" column
        columns_to_plot = [column for column in group_df.columns if '_weight' in column or '_score' in column]

        # Plot the line plots for each year
        for i, (index, row) in enumerate(group_df.iterrows()):
            # Extract the year and the values for the criteria columns
            year = row['year']
            # Calculate the color
            color = colormap(normalize(i))
            linestyle = line_styles[i % len(line_styles)]
            marker = markers[i % len(markers)]
            # Plot the criteria performances for this year
            ax[j].plot(columns_to_plot, row[columns_to_plot], linestyle=linestyle, marker=marker, label=str(year), color=color)

        # Calculate specific percentiles
        percentiles = [100, 90, 50, 10, 0]
        percentile_values = np.percentile(group_df['weighted_score_weight'], percentiles)
        # Highlight and annotate the percentiles
        for percentile, value in zip(percentiles, percentile_values):
            closest_year = group_df.iloc[(group_df['weighted_score_weight']-value).abs().argsort()[:1]]['year'].values[0]
            #ax[j].axhline(y=value, linestyle='--', color='gray', alpha=0.7)
            ax[j].text(len(columns_to_plot)-0.7 , value, f'{percentile}th percentile: {value:.2f}\n{closest_year}', verticalalignment='center')
            
        # Add legend and labels
        ax[j].legend(loc='upper left', bbox_to_anchor=(1.3 ,1), title='Year', ncol=2)
        ax[j].set_title(f'Group: {group_name}')

        
        # Improve the layout
        plt.xticks(rotation=45, ha="right")  # Rotate the x-axis labels for better readability
        plt.tight_layout()  # Adjust the layout to make room for the legend and x-axis labels

        plt.savefig(snakemake.output.plot_pick)

if __name__ == "__main__":
    # Detect running outside of snakemake and mock snakemake for testing
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('compare_weatheryears', policy=["ref", "res100"], palette='p1', zone='DE', year='2025', participation='10', weather_year='1980')

    
    summary_files = snakemake.input.summary_files # this is the list of files from the expand function
    tech_colors = snakemake.config['tech_colors']
    
    concatenated_df = concatenate_summary_files(summary_files)
    concatenated_df.to_csv(snakemake.output.summary_csv)    
    pick_df = concatenated_df
    plot_system_invest(concatenated_df, snakemake.output.plot_invest, tech_colors)
    plot_system_cfe(concatenated_df, snakemake.output.plot_cfe)
    plot_objective(concatenated_df, snakemake.output.plot_objective)

    rldc_files = snakemake.input.rldc_files
    plot_rldc(rldc_files)
    pick_wy = pick_wy_table(pick_df)
    plot_wy(pick_wy)

