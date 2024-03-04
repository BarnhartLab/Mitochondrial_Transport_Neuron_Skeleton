
import argparse
import numpy as np
import h5py
import os
import pandas as pd
from typing import Dict
from skimage import io
"""
RUNME_tif_to_h5.py
Performs a bulk conversion of tif image stacks in a directory to hdf5
- Only works with tifs for now
- Each stack is converted to an h5 file with a 'dataset' for each of its channels
- tif filenames (and their new names) must be unique to prevent any mix-up 
- Files are renamed and channels are labelled according to a csv file with columns:
 ['old_name', 'new_name', 'ignore', 'ch0', ..., 'chn']. 
 TODO: make spec file optional
"""


def main():

    cfg = parse_these_args()
    ### Find all tifs located within the parent directory ###
    full_dir = os.path.expanduser(cfg.dir)

    if not os.path.isdir(full_dir):
        raise FileNotFoundError(f"Directory not found: {full_dir}")
    else:
        dir_files = []
        for root, subdirs, files in os.walk(full_dir):  # recursive walk from the parent dir
            for f in files:
                if f[-4:] == '.tif':
                    #print(os.path.join(root, f))
                    dir_files.append(os.path.join(root, f))

    if cfg.spec != '':
        with open(os.path.expanduser(cfg.spec)) as csv_file:
            spec = pd.read_csv(csv_file)
    else:
        raise FileNotFoundError("A specification file is needed")
        # TODO: a default option where no csv is passed, all tifs converted to h5 with the same file name

    ### CHECK SPEC FILE ###
    if len(spec) != spec['old_name'].nunique():
        print(f"The following files share the same name: "
              f"{[f for f, n in spec['old_name'].value_counts().items() if n > 1]}")
        raise Exception("Spec csv: Input filenames are not unique. Find and change duplicate entries in old_name")
    elif len(spec) != spec['new_name'].nunique():
        print(f"The following files share the same name: "
              f"{[f for f, n in spec['new_name'].values_counts().items() if n > 1]}")
        raise Exception("Spec csv: Output filenames are not unique. Find and change duplicate entries in new_name")
    else:
        print('Old and new filenames are unique')

    # Make a new directory for the h5s
    output_dir = os.path.join(full_dir, 'h5_files')
    os.mkdir(output_dir)
    print(f"Converted files will be written to {output_dir}")

    path_list = []
    for i, row in spec.iterrows():
        im_path = [p for p in dir_files if row['old_name'] == os.path.split(p)[-1]]

        if row['ignore']:  # spec file says to ignore this
            path_list.append('ignore')
            continue
        else:
            if len(im_path) > 1:
                raise Exception(f"Multiple files named {row['old_name']} "
                                f"were found under the parent directory: {print(im_path)}")
            elif len(im_path) == 0:
                raise FileNotFoundError(f"{row['old_name']} not found under parent directory")
            else:
                this_path = im_path[0]
                path_list.append(this_path)

                stack = io.imread(this_path)
                color_axis = np.argmin(stack.shape)  # assume that the color axis is shorter than x, y, z
                stack = np.moveaxis(stack, color_axis, -1)  # Make sure colors are the last dim

                chan_cols = [str(k) for k in row.axes[0] if str(k[0]) == 'c']
                channel_labels = [label for label in row[chan_cols] if str(label) != 'nan']
                print(channel_labels)
                print(stack.shape)
                write_stack_to_h5(stack, os.path.join(output_dir, row['new_name']), channel_labels)


def write_stack_to_h5(im, save_path, channel_labels=None):

    if channel_labels is None:
        channel_labels = [i for i in range(0, im.shape[-1])]
    else:
        if len(channel_labels) != im.shape[-1]:
            print(channel_labels)
            print(im.shape)
            raise Exception(f"Number of channel_labels for {save_path} does not"
                            f"match the number of channels in the tiff")

    f = h5py.File(save_path, "w")
    for i, ch in enumerate(channel_labels):
        dset = f.create_dataset(ch, data=im[:, :, :, i], dtype='uint8')
    f.close()

def parse_these_args() -> Dict:

    parser = argparse.ArgumentParser(description='Reformat tiff stacks based on a specification file')
    parser.add_argument('dir', type=str, help='Parent directory containing tifs or subdirectories containing tifs')
    parser.add_argument('--spec', '-s', type=str, help='CSV specification file where the first column are old file'
                                                       'names, followed by columns describing how channels should be '
                                                       'labelled, and then the new names these files should be '
                                                       'given', default='')
    return parser.parse_args()


main()

# data_dir = os.path.expanduser('~/barnhart-lab/data/MCFO-morphology/image-data/multi-channel/')
# dsets = [d for d in os.listdir(data_dir) if d[-3:] == 'tif']
# print(dsets)
# for this_dset in dsets:
#     im_path = data_dir + this_dset
#     im = io.imread(im_path)
#     stack_path = data_dir + this_dset.split('.')[0] + '_pngs/'
#     os.mkdir(stack_path)
#     print(f'PNGs written to: {stack_path}')
#     for z in range(0, im.shape[0]):
#         fname = stack_path + this_dset.split('.')[0] + f'_{z}.png'
#         io.imsave(fname, im[z])