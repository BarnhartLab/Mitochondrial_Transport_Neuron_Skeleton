
import matplotlib.pyplot as plt
import mpl_interactions.ipyplot as iplt
import h5py
import os.path
from typing import Tuple, Dict, Union
import numpy as np
import pandas as pd
from skimage import io
import skimage.morphology
from typing import List

from fig_utils import rand_cmap
plt.style.use('mystyle.mplstyle')


def main():

    data_dir = os.path.expanduser('~/Barnhart/data/mito-seg/211206/')

    ### Make paired lists of paths to image data and pmaps ###
    im_data = []
    prob_maps = []
    for f in os.listdir(os.path.join(data_dir, 'h5_files')):
        if 'Probabilities.h5' in str(f):
            prob_maps.append(f)
        else:
            im_data.append(f)

    if len(im_data) != len(prob_maps):
        raise Exception(f'Number of image files does not match number of probability maps. Please check {data_dir}')

    file_pairs = zip(sorted(im_data), sorted(prob_maps))
    #########################################################
    i = 0
    for im_f, p_f in file_pairs:
        dset_num = str(im_f).split('.')[0]
        #assert(im_f.split('.')[0] == p_f.split('.')[0])
        print(im_f)
        print(p_f)
        im = read_h5_multi_channel(os.path.join(data_dir, 'h5_files', im_f), channels='tdTomato')
        p = read_h5_pmap(os.path.join(data_dir, 'h5_files', p_f))

        print(im.shape)
        print(p.shape)

        fig, axes = plt.subplots(2, 2, figsize=[3.8, 3.8], sharex=True, sharey=True)
        ax = axes.flatten()

        ax[0].set_title(f"im - {dset_num}")
        ax[0].imshow(np.max(im, axis=0), cmap='hot')
        ax[0].axes.xaxis.set_visible(False)
        ax[0].axes.yaxis.set_visible(False)

        ax[1].set_title(f"pmap")
        ax[1].imshow(np.max(p, axis=0), cmap='cool')
        ax[1].axes.xaxis.set_visible(False)
        ax[1].axes.yaxis.set_visible(False)

        # closing_size = 5
        # footprint = skimage.morphology.ball(closing_size)
        # p_morphed = skimage.morphology.closing(p, footprint)

        # ax[2].set_title("Closed")
        # ax[2].imshow(np.max(p_morphed, axis=0), cmap='binary')

        ### Binarize ###
        bin_thresh = 0.70
        bin_mask = binarize(p, bin_thresh)
        ax[2].set_title(f"bin_mask, p>={bin_thresh}")
        ax[2].imshow(np.max(bin_mask, axis=0), cmap='binary')
        ax[2].axes.xaxis.set_visible(False)
        ax[2].axes.yaxis.set_visible(False)
        ### Binary closing ###
        closing_size = 10
        footprint = skimage.morphology.ball(closing_size)
        bin_morphed = skimage.morphology.binary_closing(bin_mask, footprint)
        ### Label connected components ###
        label_mat, n_labels = skimage.morphology.label(bin_morphed, background=0, return_num=True)
        labels, n_px = np.unique(label_mat, return_counts=True)
        labels_npx = pd.Series(n_px, index=labels).sort_values(ascending=False)
        ordered_labels = labels_npx.index

        largest_label = ordered_labels[1]
        neurite_mask = label_mat == largest_label

        ax[3].set_title(f"neurite_mask, struct={closing_size}")
        ax[3].imshow(np.max(neurite_mask, axis=0), cmap='binary')
        ax[3].axes.xaxis.set_visible(False)
        ax[3].axes.yaxis.set_visible(False)

        plt.savefig(os.path.join(data_dir, '211207_z-proj', f'{str(dset_num)}.png'))
        plt.close(fig)

        ### SAVE new H5 containing all the processing results ###
        results = {'im': im, 'pmap': p, 'bin_mask': bin_mask, 'neurite_mask': neurite_mask}
        save_data_h5(os.path.join(data_dir, '211207_processed-h5'), file_name=dset_num, dsets=results)


def binarize(im: np.array, min_thresh: float) -> np.array:
    binary_map = np.zeros(im.shape, dtype=int)
    binary_map[im >= min_thresh] = 1
    return binary_map


def read_h5_pmap(fp: str) -> np.ndarray:
    return read_h5_multi_channel(fp, channels='pmap')[:, :, :, 0]


def read_h5_multi_channel(fp: str, channels: Union[str, List] = None) -> np.ndarray:

    f = h5py.File(fp, "r")
    arr = []

    if not channels:
        arr = np.array([f[ch] for ch in f.keys()])
        print(f'Channels not specified. All are included: {list(f.keys())}')
    elif isinstance(channels, str):
        arr = np.array(f[channels])
    elif isinstance(channels, list):
        arr = np.array(f[[ch for ch in channels]])
    else:
        print("~")
        raise ValueError("Expecting a string channel label or list of strings")

    f.close()
    return np.squeeze(np.array(arr))  # singleton dimension removed if len(channels) == 1


def save_data_h5(save_dir, file_name: str, dsets: Dict) -> None:
    """
    TODO: save h5 chunked
    :param save_dir: Directory to save h5 file
    :param file_name:
    :param dsets:
    :return:
    """
    output_file = os.path.join(save_dir, f'{file_name}.h5')

    f = h5py.File(output_file, 'w')
    for label, d in dsets.items():
        f.create_dataset(label, data=d)
    f.close()
    print(f'File written to: {output_file}')


main()

