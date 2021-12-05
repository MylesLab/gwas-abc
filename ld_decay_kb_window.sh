#!/bin/bash
#SBATCH --time=3-00:00:00
#SBATCH --mem=64G
#SBATCH --output=ld_decay_kb_window.out
#SBATCH --job-name=ld_decay_kb_window

Rscript ld_decay_kb_window.R