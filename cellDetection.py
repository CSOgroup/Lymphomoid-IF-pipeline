import os
import argparse
import tifffile
import numpy as np
import pandas as pd
import xmltodict
import shutil
from tqdm.auto import tqdm
import subprocess

from utils import run_command

parser = argparse.ArgumentParser()

parser.add_argument('--sample_names', type=str, nargs='+', required=True)
parser.add_argument('--dir', type=str, required=True)
parser.add_argument('--channel_info_path', type=str, required=True)
parser.add_argument('--deepcell_path', type=str, default='/mnt/data2/shared/Lymphomoid-IF-software/deepcell.sif')
parser.add_argument('--nucleus_channel', type=str, default='DAPI')

args = parser.parse_args()

channels_info = pd.read_csv(args.channel_info_path, delimiter='\t')
channels_info.index = channels_info['Channel_name']
channels_info['Order'] = -np.ones(channels_info.shape[0], dtype=int)

for filename in tqdm(args.sample_names):
    name = filename.replace('.ome.tif', '')
    print(f'Processing sample {name}...')

    output_dir = os.path.join(args.dir, name)

    image_path = os.path.join(output_dir, 'registration', f'{name}.ome.tif')
    channels_dir = os.path.join(output_dir, 'channels')
    masks_dir = os.path.join(output_dir, 'segmentation', 'tmp')

    with tifffile.TiffFile(image_path) as tif:
        
        channels_ordered = [c['@Name'] for c in xmltodict.parse(tif.ome_metadata)['OME']['Image']['Pixels']['Channel']]
        
        markers = pd.DataFrame(np.column_stack((np.zeros(len(channels_ordered), dtype=int), channels_ordered)), columns=['cycle', 'marker_name'])
        markers.to_csv(os.path.join(output_dir, 'markers.csv'), index=False)
        
        for i, channel in enumerate(channels_ordered):
            channels_info.loc[channel, 'Order'] = i
        
        nucleus_channel_index = channels_info.loc[args.nucleus_channel]['Order']
        membrane_channel_indexes = channels_info.loc[channels_info['Cellular_location'].str.lower() == 'cytoplasm']['Order'].values
    
        os.makedirs(masks_dir, exist_ok=True)
        os.makedirs(channels_dir, exist_ok=True)

        image = tif.asarray()

        # Select nucleus and membrane channels from the full .tif image
        tifffile.imsave(os.path.join(channels_dir, f'{name}_nucleus.tif'), image[nucleus_channel_index])
        tifffile.imsave(os.path.join(channels_dir, f'{name}_membrane.tif'), np.max(image[membrane_channel_indexes], axis=0))

        # Perform nucleus and whole-cell segmentation using DeepCell's Mesmer
        run_command(f"singularity exec --env TF_CPP_MIN_LOG_LEVEL=2 --cleanenv --no-home --bind {output_dir}:/data --nv {args.deepcell_path} /usr/local/bin/python /usr/src/app/run_app.py mesmer --nuclear-image {os.path.join('/data', 'channels', f'{name}_nucleus.tif')} --output-directory {os.path.join('/data', 'segmentation', 'tmp')} --output-name {name}_nuclear_mask.tif --compartment nuclear --squeeze ")
        run_command(f"singularity exec --env TF_CPP_MIN_LOG_LEVEL=2 --cleanenv --no-home --bind {output_dir}:/data --nv {args.deepcell_path} /usr/local/bin/python /usr/src/app/run_app.py mesmer --nuclear-image {os.path.join('/data', 'channels', f'{name}_nucleus.tif')} --membrane-image {os.path.join('/data', 'channels', f'{name}_membrane.tif')} --output-directory {os.path.join('/data', 'segmentation', 'tmp')} --output-name {name}_whole_cell_mask.tif --compartment whole-cell --squeeze")

        print(f'Segmentation completed for sample {name}')

        shutil.rmtree(channels_dir)
print()
print('Channels extracted successfully.')