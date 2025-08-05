# RNA-Seq Bioinformatics Pipeline

This repository contains scripts and configurations for RNA-Seq analysis including de novo assembly, scaffolding, and visualization.

## Structure

- `scripts/`: Shell and Python scripts for each stage
- `envs/`: Conda environment files
- `notebooks/`: Jupyter Notebooks
- `results/`: Output visualizations like Circos plots

## Usage

1. Create the environment:
```bash
conda env create -f envs/trinity_env.yml
conda activate trinity_env
```

2. Run Trinity:
```bash
qsub scripts/trinity/run_trinity_M-38.pbs
```

3. Generate Circos plot:
```bash
python scripts/circos/prepare_links.py scaffolds_vs_ref.paf > scripts/circos/links.txt
circos -conf scripts/circos/circos.conf
```
