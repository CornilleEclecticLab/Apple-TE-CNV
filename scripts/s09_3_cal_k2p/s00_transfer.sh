#!/bin/bash
#SBATCH --job-name=ATIP_BAM_NeoTransfer
#SBATCH --output=../../logs/ATIP_BAM_NeoTransfer-err2.txt
#SBATCH --error=../../logs/ATIP_BAM_NeoTransfer-out2.txt
#SBATCH --time=7-00:00:00
#SBATCH -p compute
##SBATCH -c 1
##SBATCH --mem=2G
#SBATCH -c 20
#SBATCH --mem=45G

N=15

# Remote host and base directory holding the BAM files produced by the
# SortMeRNA + Stampy mapping step. Set them for your own system, either here
# or in the environment before running this script.
REMOTE="${REMOTE:-user@cluster}"
REMOTE_STAGE="${REMOTE_STAGE:-/path/to/mapping/stage}"

# There are 94/110 samples from Duan and from 134 ATIP

#rsync -avh --partial --progress ${REMOTE}:${REMOTE_STAGE}/input_fastq/FILE_FOR_SORTMERNA/done/Sie*bam . #ERROR

############
# from Duan
#rsync -avh --partial --progress --ignore-existing ${REMOTE}:${REMOTE_STAGE}/Database_ATIP/Mapping/run_sortmerna_stampy_Duan/s02_run_sortmerna_stampy_output/*/Dom_29*_sorting.bam ../../data/BAM/originel/Duan/

#Cer Cor Dom16-35 Hyb_02 Roc

############
# from ATIP
#rsync -avh --partial --progress --ignore-existing ${REMOTE}:${REMOTE_STAGE}/Database_ATIP/Mapping/s02_run_sortmerna_stampy_output/*/F2*_sorting.bam ../../data/BAM/originel/ATIP
COMMANDS=(
# "rsync -avh --partial --progress --ignore-existing ${REMOTE}:${REMOTE_STAGE}/Database_ATIP/Mapping/s02_run_sortmerna_stampy_output/*/A022*_sorting.bam ../../data/BAM/originel/ATIP/"
# "rsync -avh --partial --progress --ignore-existing ${REMOTE}:${REMOTE_STAGE}/Database_ATIP/Mapping/s02_run_sortmerna_stampy_output/*/A023*_sorting.bam ../../data/BAM/originel/ATIP/"
# "rsync -avh --partial --progress --ignore-existing ${REMOTE}:${REMOTE_STAGE}/Database_ATIP/Mapping/s02_run_sortmerna_stampy_output/*/A024*_sorting.bam ../../data/BAM/originel/ATIP/"
# "rsync -avh --partial --progress --ignore-existing ${REMOTE}:${REMOTE_STAGE}/Database_ATIP/Mapping/s02_run_sortmerna_stampy_output/*/A025*_sorting.bam ../../data/BAM/originel/ATIP/"
# "rsync -avh --partial --progress --ignore-existing ${REMOTE}:${REMOTE_STAGE}/Database_ATIP/Mapping/s02_run_sortmerna_stampy_output/*/A026*_sorting.bam ../../data/BAM/originel/ATIP/"
# "rsync -avh --partial --progress --ignore-existing ${REMOTE}:${REMOTE_STAGE}/Database_ATIP/Mapping/s02_run_sortmerna_stampy_output/*/A027*_sorting.bam ../../data/BAM/originel/ATIP/"
# "rsync -avh --partial --progress --ignore-existing ${REMOTE}:${REMOTE_STAGE}/Database_ATIP/Mapping/s02_run_sortmerna_stampy_output/*/A028*_sorting.bam ../../data/BAM/originel/ATIP/"
# "rsync -avh --partial --progress --ignore-existing ${REMOTE}:${REMOTE_STAGE}/Database_ATIP/Mapping/s02_run_sortmerna_stampy_output/*/A029*_sorting.bam ../../data/BAM/originel/ATIP/"
# "rsync -avh --partial --progress --ignore-existing ${REMOTE}:${REMOTE_STAGE}/Database_ATIP/Mapping/s02_run_sortmerna_stampy_output/*/A031*_sorting.bam ../../data/BAM/originel/ATIP/"
# "rsync -avh --partial --progress --ignore-existing ${REMOTE}:${REMOTE_STAGE}/Database_ATIP/Mapping/s02_run_sortmerna_stampy_output/*/A032*_sorting.bam ../../data/BAM/originel/ATIP/"
# "rsync -avh --partial --progress --ignore-existing ${REMOTE}:${REMOTE_STAGE}/Database_ATIP/Mapping/s02_run_sortmerna_stampy_output/*/A033*_sorting.bam ../../data/BAM/originel/ATIP/"
# "rsync -avh --partial --progress --ignore-existing ${REMOTE}:${REMOTE_STAGE}/Database_ATIP/Mapping/s02_run_sortmerna_stampy_output/*/A034*_sorting.bam ../../data/BAM/originel/ATIP/"
# "rsync -avh --partial --progress --ignore-existing ${REMOTE}:${REMOTE_STAGE}/Database_ATIP/Mapping/s02_run_sortmerna_stampy_output/*/A035*_sorting.bam ../../data/BAM/originel/ATIP/"
# "rsync -avh --partial --progress --ignore-existing ${REMOTE}:${REMOTE_STAGE}/Database_ATIP/Mapping/s02_run_sortmerna_stampy_output/*/A036*_sorting.bam ../../data/BAM/originel/ATIP/"
# "rsync -avh --partial --progress --ignore-existing ${REMOTE}:${REMOTE_STAGE}/Database_ATIP/Mapping/s02_run_sortmerna_stampy_output/*/A037*_sorting.bam ../../data/BAM/originel/ATIP/"

# "rsync -avh --partial --progress --ignore-existing ${REMOTE}:${REMOTE_STAGE}/Database_ATIP/Mapping/s02_run_sortmerna_stampy_output/*/A038*_sorting.bam ../../data/BAM/originel/ATIP/"
# "rsync -avh --partial --progress --ignore-existing ${REMOTE}:${REMOTE_STAGE}/Database_ATIP/Mapping/s02_run_sortmerna_stampy_output/*/A039*_sorting.bam ../../data/BAM/originel/ATIP/"
# "rsync -avh --partial --progress --ignore-existing ${REMOTE}:${REMOTE_STAGE}/Database_ATIP/Mapping/s02_run_sortmerna_stampy_output/*/BACC02*_sorting.bam ../../data/BAM/originel/ATIP/"
# "rsync -avh --partial --progress --ignore-existing ${REMOTE}:${REMOTE_STAGE}/Database_ATIP/Mapping/s02_run_sortmerna_stampy_output/*/BACC03*_sorting.bam ../../data/BAM/originel/ATIP/"
# "rsync -avh --partial --progress --ignore-existing ${REMOTE}:${REMOTE_STAGE}/Database_ATIP/Mapping/s02_run_sortmerna_stampy_output/*/BACC04*_sorting.bam ../../data/BAM/originel/ATIP/"
# "rsync -avh --partial --progress --ignore-existing ${REMOTE}:${REMOTE_STAGE}/Database_ATIP/Mapping/s02_run_sortmerna_stampy_output/*/BACC05*_sorting.bam ../../data/BAM/originel/ATIP/"
# "rsync -avh --partial --progress --ignore-existing ${REMOTE}:${REMOTE_STAGE}/Database_ATIP/Mapping/s02_run_sortmerna_stampy_output/*/BACC06*_sorting.bam ../../data/BAM/originel/ATIP/"
# "rsync -avh --partial --progress --ignore-existing ${REMOTE}:${REMOTE_STAGE}/Database_ATIP/Mapping/s02_run_sortmerna_stampy_output/*/BACC07*_sorting.bam ../../data/BAM/originel/ATIP/"
# "rsync -avh --partial --progress --ignore-existing ${REMOTE}:${REMOTE_STAGE}/Database_ATIP/Mapping/s02_run_sortmerna_stampy_output/*/BACC08*_sorting.bam ../../data/BAM/originel/ATIP/"
# "rsync -avh --partial --progress --ignore-existing ${REMOTE}:${REMOTE_STAGE}/Database_ATIP/Mapping/s02_run_sortmerna_stampy_output/*/BACC09*_sorting.bam ../../data/BAM/originel/ATIP/"
# "rsync -avh --partial --progress --ignore-existing ${REMOTE}:${REMOTE_STAGE}/Database_ATIP/Mapping/s02_run_sortmerna_stampy_output/*/BACC10*_sorting.bam ../../data/BAM/originel/ATIP/"
# "rsync -avh --partial --progress --ignore-existing ${REMOTE}:${REMOTE_STAGE}/Database_ATIP/Mapping/s02_run_sortmerna_stampy_output/*/BACC11*_sorting.bam ../../data/BAM/originel/ATIP/"
# "rsync -avh --partial --progress --ignore-existing ${REMOTE}:${REMOTE_STAGE}/Database_ATIP/Mapping/s02_run_sortmerna_stampy_output/*/C10*_sorting.bam ../../data/BAM/originel/ATIP/"
# "rsync -avh --partial --progress --ignore-existing ${REMOTE}:${REMOTE_STAGE}/Database_ATIP/Mapping/s02_run_sortmerna_stampy_output/*/D03*_sorting.bam ../../data/BAM/originel/ATIP/"
# "rsync -avh --partial --progress --ignore-existing ${REMOTE}:${REMOTE_STAGE}/Database_ATIP/Mapping/s02_run_sortmerna_stampy_output/*/D04*_sorting.bam ../../data/BAM/originel/ATIP/"

"rsync -avh --partial --progress --ignore-existing ${REMOTE}:${REMOTE_STAGE}/Database_ATIP/Mapping/s02_run_sortmerna_stampy_output/*/D05*_sorting.bam ../../data/BAM/originel/ATIP/"
"rsync -avh --partial --progress --ignore-existing ${REMOTE}:${REMOTE_STAGE}/Database_ATIP/Mapping/s02_run_sortmerna_stampy_output/*/D06*_sorting.bam ../../data/BAM/originel/ATIP/"
"rsync -avh --partial --progress --ignore-existing ${REMOTE}:${REMOTE_STAGE}/Database_ATIP/Mapping/s02_run_sortmerna_stampy_output/*/D07*_sorting.bam ../../data/BAM/originel/ATIP/"
"rsync -avh --partial --progress --ignore-existing ${REMOTE}:${REMOTE_STAGE}/Database_ATIP/Mapping/s02_run_sortmerna_stampy_output/*/D08*_sorting.bam ../../data/BAM/originel/ATIP/"
"rsync -avh --partial --progress --ignore-existing ${REMOTE}:${REMOTE_STAGE}/Database_ATIP/Mapping/s02_run_sortmerna_stampy_output/*/D09*_sorting.bam ../../data/BAM/originel/ATIP/"
"rsync -avh --partial --progress --ignore-existing ${REMOTE}:${REMOTE_STAGE}/Database_ATIP/Mapping/s02_run_sortmerna_stampy_output/*/F02*_sorting.bam ../../data/BAM/originel/ATIP/"
"rsync -avh --partial --progress --ignore-existing ${REMOTE}:${REMOTE_STAGE}/Database_ATIP/Mapping/s02_run_sortmerna_stampy_output/*/F03*_sorting.bam ../../data/BAM/originel/ATIP/"
"rsync -avh --partial --progress --ignore-existing ${REMOTE}:${REMOTE_STAGE}/Database_ATIP/Mapping/s02_run_sortmerna_stampy_output/*/F04*_sorting.bam ../../data/BAM/originel/ATIP/"
"rsync -avh --partial --progress --ignore-existing ${REMOTE}:${REMOTE_STAGE}/Database_ATIP/Mapping/s02_run_sortmerna_stampy_output/*/F05*_sorting.bam ../../data/BAM/originel/ATIP/"
"rsync -avh --partial --progress --ignore-existing ${REMOTE}:${REMOTE_STAGE}/Database_ATIP/Mapping/s02_run_sortmerna_stampy_output/*/F06*_sorting.bam ../../data/BAM/originel/ATIP/"
"rsync -avh --partial --progress --ignore-existing ${REMOTE}:${REMOTE_STAGE}/Database_ATIP/Mapping/s02_run_sortmerna_stampy_output/*/F07*_sorting.bam ../../data/BAM/originel/ATIP/"
"rsync -avh --partial --progress --ignore-existing ${REMOTE}:${REMOTE_STAGE}/Database_ATIP/Mapping/s02_run_sortmerna_stampy_output/*/F08*_sorting.bam ../../data/BAM/originel/ATIP/"
"rsync -avh --partial --progress --ignore-existing ${REMOTE}:${REMOTE_STAGE}/Database_ATIP/Mapping/s02_run_sortmerna_stampy_output/*/F09*_sorting.bam ../../data/BAM/originel/ATIP/"
"rsync -avh --partial --progress --ignore-existing ${REMOTE}:${REMOTE_STAGE}/Database_ATIP/Mapping/s02_run_sortmerna_stampy_output/*/SP05*_sorting.bam ../../data/BAM/originel/ATIP/"
"rsync -avh --partial --progress --ignore-existing ${REMOTE}:${REMOTE_STAGE}/Database_ATIP/Mapping/s02_run_sortmerna_stampy_output/*/SP06*_sorting.bam ../../data/BAM/originel/ATIP/"
"rsync -avh --partial --progress --ignore-existing ${REMOTE}:${REMOTE_STAGE}/Database_ATIP/Mapping/s02_run_sortmerna_stampy_output/*/SP07*_sorting.bam ../../data/BAM/originel/ATIP/"
"rsync -avh --partial --progress --ignore-existing ${REMOTE}:${REMOTE_STAGE}/Database_ATIP/Mapping/s02_run_sortmerna_stampy_output/*/SP08*_sorting.bam ../../data/BAM/originel/ATIP/"
"rsync -avh --partial --progress --ignore-existing ${REMOTE}:${REMOTE_STAGE}/Database_ATIP/Mapping/s02_run_sortmerna_stampy_output/*/SP09*_sorting.bam ../../data/BAM/originel/ATIP/"
)

# --- Controlled parallel execution ---
# Run the commands in parallel, at most $N at a time
echo "Starting parallel rsync transfers..."
for cmd in "${COMMANDS[@]}"; do
    echo "Launching: $cmd"
    bash -c "$cmd" &
    
    # Cap the number of concurrent jobs
    while (( $(jobs -r | wc -l) >= N )); do
        sleep 5
    done
done

# Wait for every transfer to finish
wait
echo "✅ All transfers completed."