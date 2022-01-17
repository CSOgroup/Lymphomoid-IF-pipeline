import os
import argparse
import tifffile
import numpy as np
from tqdm.auto import tqdm

parser = argparse.ArgumentParser()

parser.add_argument('--sample_names', type=str, nargs='+', default=None)
parser.add_argument('--dir', type=str, default='/mnt/ndata/varrone/Data/graviton/HLS')
parser.add_argument('--deepcell_path', type=str, default='/mnt/data2/shared/Lymphomoid-IF-software/deepcell.sif')
parser.add_argument('--nucleus_channel', type=int, default=0)
parser.add_argument('--membrane_channels', type=int, nargs='+', default=[1,3,4,5])

args = parser.parse_args()


for filename in tqdm(args.sample_names):
    name = filename.replace('.tif', '')
    print(f'Processing sample {name}...')

    output_dir = os.path.join(args.dir, name)
    os.makedirs(output_dir, exist_ok=True)

    image_path = os.path.join(output_dir, 'registration', f'{name}.tif')
    channels_dir = os.path.join(output_dir, 'channels')
    masks_dir = os.path.join(output_dir, 'segmentation', 'tmp')
    
    all_chs = tifffile.imread(image_path)
    os.makedirs(masks_dir, exist_ok=True)
    os.makedirs(channels_dir, exist_ok=True)

    # Select nucleus and membrane channels from the full .tif image
    tifffile.imsave(os.path.join(channels_dir, f'{name}_nucleus.tif'), all_chs[args.nucleus_channel])
    tifffile.imsave(os.path.join(channels_dir, f'{name}_membrane.tif'), np.max(all_chs[args.membrane_channels], axis=0))

    # Perform nucleus and whole-cell segmentation using DeepCell's Mesmer
    result = os.popen(f"singularity exec --env TF_CPP_MIN_LOG_LEVEL=2 --cleanenv --no-home --bind {output_dir}:/data --nv {args.deepcell_dir} /usr/local/bin/python /usr/src/app/run_app.py mesmer --nuclear-image {os.path.join('/data', 'channels', f'{name}_nucleus.tif')} --output-directory {os.path.join('/data', 'segmentation', 'tmp')} --output-name {name}_nuclear_mask.tif --compartment nuclear --squeeze ").read()
    result = os.popen(f"singularity exec --env TF_CPP_MIN_LOG_LEVEL=2 --cleanenv --no-home --bind {output_dir}:/data --nv {args.deepcell_dir} /usr/local/bin/python /usr/src/app/run_app.py mesmer --nuclear-image {os.path.join('/data', 'channels', f'{name}_nucleus.tif')} --membrane-image {os.path.join('/data', 'channels', f'{name}_membrane.tif')} --output-directory {os.path.join('/data', 'segmentation', 'tmp')} --output-name {name}_whole_cell_mask.tif --compartment whole-cell --squeeze").read()

    print(f'Segmentation completed for sample {name}')
print()
print('Channels extracted successfully.')