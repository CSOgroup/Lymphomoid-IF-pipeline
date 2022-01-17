import os
import argparse
import tifffile
import numpy as np
import pandas as pd
from tqdm.auto import tqdm
import shutil

parser = argparse.ArgumentParser()

parser.add_argument('--sample_names', type=str, default=None)
parser.add_argument('--dir', type=str, default='/mnt/ndata/varrone/Data/graviton/HLS')
parser.add_argument('--nextflow_dir', type=str, default='/mnt/data2/shared/Lymphomoid-IF-software/nextflow')
parser.add_argument('--markers', type=str, nargs='*', default=['DAPI', 'CD20', 'Ki67', 'CD4', 'CD8', 'CD68'])

args = parser.parse_args()

filenames = [args.sample_name] if args.sample_name else os.listdir(args.masks_dir)

for name in tqdm(filenames):
    print(name)
    output_dir = os.path.join(args.dir, name)
    
    masks_dir = os.path.join(output_dir, 'segmentation')

    if not os.path.exists(os.path.join(masks_dir, f'mesmer-{name}', f'cytoplasm.tif')) or \
        not os.path.exists(os.path.join(masks_dir, f'mesmer-{name}', f'nuclei.tif')):

        mask_nuclei = tifffile.imread(os.path.join(masks_dir, 'tmp', f'{name}_nuclear_mask.tif'))
        mask_whole_cell = tifffile.imread(os.path.join(masks_dir, 'tmp', f'{name}_whole_cell_mask.tif'))

        # Remove the nucleus mask from the whole-cell mask to obtain the cytoplasm
        mask_cytoplasm = mask_whole_cell*~(mask_nuclei > 0)

        # Select nucleus masks where there is a match with the whole cell masks
        mask_nuclei = mask_whole_cell*(mask_nuclei > 0)
        
        os.makedirs(os.path.join(masks_dir, f'mesmer-{name}'), exist_ok=True)
        tifffile.imsave(os.path.join(masks_dir, f'mesmer-{name}', f'cytoplasm.tif'), mask_cytoplasm)
        tifffile.imsave(os.path.join(masks_dir, f'mesmer-{name}', f'nuclei.tif'), mask_nuclei)

    if os.path.exists(os.path.join(masks_dir, 'tmp')):
        shutil.rmtree(os.path.join(masks_dir, 'tmp'))
    
    # Extract the mean intensity of each marker from the mask of each cell with MCMICRO
    if not os.path.exists(os.path.join(output_dir, 'quantification')):
        markers = pd.DataFrame(np.column_stack((np.zeros(len(args.markers), dtype=int), args.markers)), columns=['cycle', 'marker_name'])
        markers.to_csv(os.path.join(output_dir, 'markers.csv'), index=False)

        os.chdir(args.dir)
        result = os.popen(f"{args.nextflow_dir} run labsyspharm/mcmicro --profile O2 --in {name} --start-at quantification --stop-at quantification --probability-maps mesmer --quant-opts '--masks nuclei.tif cytoplasm.tif'").read()
        print(result)
