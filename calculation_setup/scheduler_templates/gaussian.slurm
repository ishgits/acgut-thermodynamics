#!/usr/bin/env bash
# Customize resources and your site's Gaussian module/account before submission.
#SBATCH --job-name=atcgu_example
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=24:00:00
set -euo pipefail
# Example module command; replace with your installation's setup.
# module load gaussian16
: "${ATCGU_INPUT:?Set ATCGU_INPUT to the reviewed absolute .com input path}"
: "${ATCGU_CHECKOUT:?Set ATCGU_CHECKOUT to this repository checkout}"
bash "$ATCGU_CHECKOUT/calculation_setup/scheduler_templates/run_gaussian.sh" "$ATCGU_INPUT"
