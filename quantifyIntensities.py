import os
import argparse
import tifffile
import numpy as np
import pandas as pd
import subprocess
import sys
from tqdm.auto import tqdm
import shutil

def run_command(command):
    process = subprocess.Popen(command, shell=True, bufsize=1,
                stdout=subprocess.PIPE, stderr = subprocess.STDOUT,encoding='utf-8', errors = 'replace' )
    while True:
        realtime_output = process.stdout.readline()
        if realtime_output == '' and process.poll() is not None:
            break
        if realtime_output:
            print(realtime_output.strip(), flush=False)
            sys.stdout.flush()
    if process.returncode:
        sys.exit()

parser = argparse.ArgumentParser()

parser.add_argument('--sample_names', type=str, nargs='+', required=True)
parser.add_argument('--dir', type=str, required=True)
parser.add_argument('--channel_info_path', type=str, required=True)
parser.add_argument('--nextflow_dir', type=str, default='/mnt/data2/shared/Lymphomoid-IF-software/nextflow')

args = parser.parse_args()

working_dir = os.getcwd()
print(working_dir)

channels_info = pd.read_csv(args.channel_info_path, delimiter='\t')
channels_info.index = channels_info['Channel_name']
nuclear_markers = channels_info[channels_info['Cellular_location'].str.lower() == 'nucleus']['Channel_name']
membrane_markers = channels_info[channels_info['Cellular_location'].str.lower() == 'cytoplasm']['Channel_name']


for name in tqdm(args.sample_names):
    print(name)
    output_dir = os.path.join(args.dir, name)
    
    masks_dir = os.path.join(output_dir, 'segmentation')

    # if not os.path.exists(os.path.join(masks_dir, f'mesmer-{name}', f'cytoplasm.tif')) or \
    #     not os.path.exists(os.path.join(masks_dir, f'mesmer-{name}', f'nuclei.tif')):

    #     print('Loading segmentation masks...')
    #     mask_nuclei = tifffile.imread(os.path.join(masks_dir, 'tmp', f'{name}_nuclear_mask.tif'))
    #     mask_whole_cell = tifffile.imread(os.path.join(masks_dir, 'tmp', f'{name}_whole_cell_mask.tif'))

    #     # Remove the nucleus mask from the whole-cell mask to obtain the cytoplasm
    #     print('Generating cytoplasm mask...')
    #     mask_cytoplasm = mask_whole_cell*~(mask_nuclei > 0)

    #     # Select nucleus masks where there is a match with the whole cell masks
    #     mask_nuclei = mask_whole_cell*(mask_nuclei > 0)
        
    #     print('Saving final masks...')
    #     os.makedirs(os.path.join(masks_dir, f'mesmer-{name}'), exist_ok=True)
    #     tifffile.imsave(os.path.join(masks_dir, f'mesmer-{name}', f'cytoplasm.tif'), mask_cytoplasm)
    #     tifffile.imsave(os.path.join(masks_dir, f'mesmer-{name}', f'nuclei.tif'), mask_nuclei)

    if os.path.exists(os.path.join(masks_dir, 'tmp')):
        shutil.rmtree(os.path.join(masks_dir, 'tmp'))
    
    # Extract the mean intensity of each marker from the mask of each cell with MCMICRO
    if not os.path.exists(os.path.join(output_dir, 'quantification')):
        os.chdir(args.dir)
        print('Start quantification...')
        run_command(f"{args.nextflow_dir} run labsyspharm/mcmicro --profile singularity --in {name} --start-at quantification --stop-at quantification --probability-maps mesmer --quant-opts '--masks nuclei.tif cytoplasm.tif'")

        os.chdir(working_dir)

        sample_nuclei = pd.read_csv(os.path.join(output_dir, 'quantification', f'mesmer-{name}_nuclei.csv'))
        sample_nuclei = sample_nuclei.loc[:, ~sample_nuclei.columns.isin(membrane_markers)]
        
        sample_cytoplasm = pd.read_csv(os.path.join(output_dir, 'quantification', f'mesmer-{name}_cytoplasm.csv'))
        sample_cytoplasm = sample_cytoplasm.loc[:, ~sample_cytoplasm.columns.isin(nuclear_markers)]

        sample = pd.merge(sample_nuclei, sample_cytoplasm, on='CellID', how='left', suffixes=("_nucleus", "_cytoplasm"))
        sample = sample.rename(columns={row['Channel_name']: row['Antibody'] for _, row in channels_info.iterrows()})
        sample.index = sample.index.astype(str)
        sample.to_csv(os.path.join(output_dir, 'quantification', f'mesmer-{name}_merged.csv'), index=False)



