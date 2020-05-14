#!/bin/sh
# note - there can be no line spaces between # SBATCH directives .
# SBATCH -- job - name = heat2D
# SBATCH -- output = heat.out
# SBATCH -- error = gpuTest_ % j . err
# SBATCH -- ntasks =1
# SBATCH -- cpus - per - task =1
# SBATCH -- time =12:00:00
# SBATCH --partition=am-148-s20
# SBATCH--qos=am-148-s20
# SBATCH--account=am-148-s20

module load cuda10.0/10.0

./transport.exe
