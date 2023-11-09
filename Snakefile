configfile: "config.yaml"

wildcard_constraints:
    policy="[\-a-zA-Z0-9\.]+",
    weather_year="([0-9]+)(-[0-9]+)?(\+([0-9]+)(-[0-9]+)?)*",


RDIR = os.path.join(config['results_dir'], config['run'])
CDIR = config['costs_dir']

def parse_year_wildcard(w):
    """
    Parse a {year} wildcard to a list of years. Parses a wildcard with years and ranges separated by '+' and '-'.
    
    Args:
    w (str): A string containing years and ranges like '1980+1990+2000-2002'.
    
    Returns:
    list: A list of years expanded from the ranges.
    """
    years = []
    parts = w.split('+')  # Split the string by '+'
    
    for part in parts:
        if '-' in part:
            # It's a range, split it into start and end, then generate the range of years.
            start, end = map(int, part.split('-'))
            # Ensure that the start year is less than the end year.
            if end < start:
                raise ValueError(f"Malformed range of years {part}.")
            years.extend(range(start, end + 1))  # Include the end year in the range.
        else:
            # It's a single year, just append it to the list.
            years.append(int(part))
    
    return sorted(years)

rule merge_all_plots:
    input: 
        expand(RDIR + "/plots/{participation}/{year}/{zone}/{palette}/{weather_year}/SUMMARY.pdf", **config['scenario'])


rule plot_summary_all_networks:
    input: 
        expand(RDIR + "/plots/{participation}/{year}/{zone}/{palette}/{weather_year}/used.pdf", **config['scenario'])


rule make_summary_all_networks:
    input: 
        expand(RDIR + "/csvs/{participation}/{year}/{zone}/{palette}/{weather_year}/summary.csv", **config['scenario'])


rule summarise_all_networks:
    input: 
        expand(RDIR + "/summaries/{participation}/{year}/{zone}/{palette}/{weather_year}/{policy}.yaml", **config['scenario'])


rule solve_all_networks:
    input: 
        expand(RDIR + "/networks/{participation}/{year}/{zone}/{palette}/{weather_year}/{policy}.nc", **config['scenario'])


rule merge_plots:
    input:
        used=RDIR + "/plots/{participation}/{year}/{zone}/{palette}/{weather_year}/used.pdf",
        config=RDIR + '/configs/config.yaml'
    output:
        final=RDIR + "/plots/{participation}/{year}/{zone}/{palette}/{weather_year}/SUMMARY.pdf"
    threads: 2
    resources: mem_mb=2000
    script:
        'scripts/merge_plots.py'


rule plot_summary:
    input:
        summary=RDIR + "/csvs/{participation}/{year}/{zone}/{palette}/{weather_year}/summary.csv",
        config=RDIR + '/configs/config.yaml'
    output:
        used=RDIR + "/plots/{participation}/{year}/{zone}/{palette}/{weather_year}/used.pdf"
    threads: 2
    resources: mem_mb=2000
    script:
        'scripts/plot_summary.py'


rule make_summary:
    input:
        expand(RDIR + "/summaries/{participation}/{year}/{zone}/{palette}/{weather_year}/{policy}.yaml",
               **config['scenario'])
    output:
        summary=RDIR + "/csvs/{participation}/{year}/{zone}/{palette}/{weather_year}/summary.csv"
    threads: 2
    resources: mem_mb=2000
    script: 'scripts/make_summary.py'


if config['solve_network'] == 'solve':
    rule solve_network:
        input:
            network2030 = 'input/v6_elec_s_37_lv1.0__3H-B-solar+p3_2030.nc',
            network2025 = 'input/elec{weather_year}_s_37_lv1.0__3H-B-solar+p3_2025.nc',
            costs2030=CDIR + "/costs_2030.csv",
            costs2025=CDIR + "/costs_2025.csv"
        output:
            network=RDIR + "/networks/{participation}/{year}/{zone}/{palette}/{weather_year}/{policy}.nc",
            grid_cfe=RDIR + "/networks/{participation}/{year}/{zone}/{palette}/{weather_year}/{policy}.csv"
        log:
            solver=RDIR + "/logs/{participation}/{year}/{zone}/{palette}/{weather_year}/{policy}_solver.log",
            python=RDIR + "/logs/{participation}/{year}/{zone}/{palette}/{weather_year}/{policy}_python.log",
            memory=RDIR + "/logs/{participation}/{year}/{zone}/{palette}/{weather_year}/{policy}_memory.log"
        threads: 12
        resources: mem=8000
        script: "scripts/solve_network.py"

rule summarise_network:
    input:
        network=RDIR + "/networks/{participation}/{year}/{zone}/{palette}/{weather_year}/{policy}.nc",
	    grid_cfe=RDIR + "/networks/{participation}/{year}/{zone}/{palette}/{weather_year}/{policy}.csv"
    output:
        yaml=RDIR + "/summaries/{participation}/{year}/{zone}/{palette}/{weather_year}/{policy}.yaml"
    threads: 2
    resources: mem_mb=2000
    script: 'scripts/summarise_network.py'


rule copy_config:
    output: RDIR + '/configs/config.yaml'
    threads: 1
    resources: mem_mb=1000
    script: "scripts/copy_config.py"


# additional rules for cluster communication -> not included into a workflow 
rule sync_solution:
    params:
        cluster="iegor.riepin@gateway.hpc.tu-berlin.de:/scratch/iegor.riepin/247-cfe/results/report"
    shell: 
        """
        rsync -uvarh --no-g {params.cluster} results/
        """

rule sync_plots:
    params:
        cluster="iegor.riepin@gateway.hpc.tu-berlin.de:/scratch/iegor.riepin/247-cfe/results/report/plots/"
    shell: 
        """
        rsync -uvarh --no-g {params.cluster} report/plots
        """


# illustrate workflow
rule dag:
     message: "Plot dependency graph of the workflow."
     output:
         dot="workflow/dag.dot",
         graph="workflow/graph.dot",
         pdf="workflow/graph.pdf"
     shell:
         """
         snakemake --rulegraph > {output.dot}
         sed -e '1,3d' < {output.dot} > {output.graph}
         dot -Tpdf -o {output.pdf} {output.graph}
         """