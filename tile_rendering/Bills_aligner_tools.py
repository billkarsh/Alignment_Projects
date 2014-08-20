from __future__ import print_function

__author__ = 'Joergen Kornfeld'

import socket
import subprocess
import random
import os
import re
import shutil
import glob
import matplotlib.pyplot as mplot
import numpy as np
import scipy
import scipy.stats
import io
import StringIO
import Image
import PIL.ImageOps
import imreg
from libtiff import TIFF
from scipy.signal.signaltools import correlate2d as c2d
from scipy.misc import imread
import scipy as sp
import sys
import cv2
import time
try:
    import stackUtilities as su
except ImportError:
    stack_utilities = False
    print("No stack utilities found, using slower image loading.")

import cPickle
from skimage import exposure


class render_job_info(object):
    """
    Class solely made for pickling of single render job info - we use the
    pickle mechanism here as a convenient way of file-system based temporary
    information.
    """
    def __init__(self):
        self.num_tiles = 0
        self.out_path = ''
        self.out_format = ''
        self.crop_px = 0
        self.invert_images = False
        self.gaussian_sigma = 0.0  # 0.0 indicates no filtering
        self.mean_normalization = True,
        self.clahe = False
        self.clahe_clip_lim = 1.0, # value hint: use between 1. and 2.
        self.clahe_size = 10,
        self.hist_stretch = False,
        self.hist_stretch_percentile = 0.5,
        self.global_xy_dims = ()
        self.affines = []
        self.tile_paths = []
        self.scale = 1.0


def single_render_job(job_file_path):
    """
    Function to be called by each qsub-ed render job. The job information is
    stored in a pickled render_job_info object, at job_file_path.

    :param job_file_path: str
        Path to pickled render_job_info instance.
    :return:
    """

    print("Render job execution started on host {0} at {1}.".format(
        socket.gethostname(), time.strftime("%Y-%m-%d %H:%M:%S")))

    with open(job_file_path, 'rb') as job_file:
        job_info = cPickle.load(job_file)

    print("Render job parameters:")
    print('===================================================================')
    print("Affines and input paths omitted for readability.")
    print("Output path: {0}".format(job_info.out_path))
    print("Output format: {0}".format(job_info.out_format))
    print("X output size: {0}".format(job_info.global_xy_dims[0]))
    print("Y output size: {0}".format(job_info.global_xy_dims[1]))
    print("Crop pixels: {0}".format(job_info.crop_px))
    print("Scaling: {0}".format(job_info.scale))
    print("Invert images: {0}".format(job_info.invert_images))
    print("Mean normalization: {0}".format(job_info.mean_normalization))
    print("Histogram stretch: {0}".format(job_info.hist_stretch))
    print("Hist stretch percentile: {0}".format(job_info.hist_stretch_percentile))
    print("CLAHE: {0}".format(job_info.clahe))
    print("CLAHE clip lim: {0}".format(job_info.clahe_clip_lim))
    print("Gaussian sigma: {0}".format(job_info.gaussian_sigma))
    print('===================================================================')


    render_tiles_in_layer(affines = job_info.affines,
                          global_xy_layer_dimensions = job_info.global_xy_dims,
                          tile_paths = job_info.tile_paths,
                          out_path = job_info.out_path,
                          out_format = job_info.out_format,
                          scale = job_info.scale,
                          crop_px = job_info.crop_px,
                          invert_images = job_info.invert_images,
                          hist_stretch = job_info.hist_stretch,
                          hist_stretch_percentile =
                          job_info.hist_stretch_percentile,
                          clahe = job_info.clahe,
                          mean_normalization=job_info.mean_normalization,
                          clahe_clip_lim = job_info.clahe_clip_lim,
                          gaussian_sigma = job_info.gaussian_sigma)


    # We remove the job file - this is an indicator that everything went
    # fine
    os.remove(job_file_path)

    return


def qsub_render_jobs():
    """
    Function that mass-submits tile rendering jobs into the SGE queue with
    qsub. Everything is specified in the code, see below.

    :return:
    """

    # Specify parameters here:
    # =========================================================================

    # Be aware of some python module dependencies
    python_exe_path = '/home/kornfeld/bin/anaconda/bin/python'
    single_render_job_script_path = \
        '/home/kornfeld/repos/skeleton-analysis/EMUtilities/'\
        '/Bills_aligner_tools.py'

    render_layer_range = [0, 15493]

    # skip layers can be useful for debugging: 100 means that only every
    # 100th layer is rendered from the total layer range; use 1 for no skipping
    skip_layers = 1

    tile_path = '/big-data/kornfeld/rawStacks/j0251/j0251/'
    affine_TXT_path = '/big-data/kornfeld/rawStacks/j0251/align_06_full/temp0' \
                      '/stack/X_A_TXT/'


    # this folder will be used for the temporary storage of the job info
    # files - each job deletes its info file after it is done
    job_info_base_path = '/home/kornfeld/tmp/'

    # specify the sun grid engine queue that should be used on the cluster
    sge_queue = 'somahead01q'

    # additional flags that get added to the qsub command
    sge_additional_flags = ''

    num_tiles = 16
    # linear interpolation between start and end layer of gaussian sigma to
    # be applied to the layers - this can be used if you have a SNR gradient
    # from the beginning to the end of your stack, e.g. due to
    # decreasing detector sensitivity over time
    gauss_sigma_gradient = [0.45, 0.60]

    # linear interpolation between start and end layer of clahe_clip_lim to
    # be applied to the layers - this can be used if you have a SNR gradient
    # from the beginning to the end of your stack, e.g. due to
    # decreasing detector sensitivity over time
    clahe_clip_lim_gradient = [1.0, 2.4]

    # This option enables the removal of linear drifts of the affines. It can
    # be useful, but you should know what you're doing when using it - it can
    # mess up the entire stack by correcting a linear that drift that is
    # actually a right one.
    remove_slow_drifts = True

    # used as a workaround to manually fix the global coordinate system after
    # slow drift removal. This can and should be automated at some point.
    global_shift = (-750,-3800)

    out_path = '/big-data/kornfeld/rawStacks/j0251/raw_registered_j0251/'
    out_filename_prepend = 'j0251_'
    out_format = 'tif'
    crop_px = 25
    invert_images = True

    # can be used for downscaling for debugging, eg. 0.5;
    scale_layer = 1.0

    clahe = True
    hist_stretch = False
    mean_normalization = True # always adisable, even with clahe enabled.
    global_xy_dims = (27000, 27000)

    # No need to change code below
    # =========================================================================









    if gauss_sigma_gradient[0] > gauss_sigma_gradient[1]:
        raise Exception('Only increasing sigmas for gaussian filtering '
                        'gradient supported.')

    # Read in all affine transformations that were calculated by Bill's
    # alignment tool before.
    affines = parse_affines_from_bills_x_a_txt(affine_TXT_path, num_tiles)

    if remove_slow_drifts:
        affines = remove_slow_drifts_in_affines(affines, global_shift)

    # get required global max dimensions; not yet implemented
    # here, we need to analyze all the affines
    # for layer_affines in affines:
    #     for tile_affine in layer_affines:
    #         tile_affine[3]
    #         tile_affine[5]

    all_tifs = glob.glob(tile_path + '*.tif')
    all_tifs.sort()

    num_layers = len(all_tifs) / num_tiles # + 1?

    if render_layer_range[1] > num_layers:
        raise Exception('Requested range not possible.')

    if render_layer_range[1] > len(affines):
        raise Exception('Not enough affines for requested rendering jobs.')

    for layer_cnt in range(render_layer_range[0], render_layer_range[1],
                           skip_layers):
        # we generate a job-info file for each job that contains all relevant
        # information

        this_job_file_info = render_job_info()

        this_job_file_info.num_tiles = num_tiles
        this_job_file_info.out_path = out_path + out_filename_prepend + \
                                                 "_0:05d}.".format(layer_cnt) \
                                      + out_format

        this_job_file_info.out_format = out_format
        this_job_file_info.crop_px = crop_px
        this_job_file_info.invert_images = invert_images
        this_job_file_info.gaussian_sigma = gauss_sigma_gradient[0] +\
            float(layer_cnt) / float(num_layers) *\
            (gauss_sigma_gradient[1]-gauss_sigma_gradient[0])

        this_job_file_info.clahe_clip_lim = clahe_clip_lim_gradient[0] +\
            float(layer_cnt) / float(num_layers) *\
            (clahe_clip_lim_gradient[1]-clahe_clip_lim_gradient[0])
        this_job_file_info.scale = scale_layer
        this_job_file_info.clahe = clahe
        this_job_file_info.global_xy_dims = global_xy_dims
        this_job_file_info.affines = affines[layer_cnt]
        this_job_file_info.tile_paths = all_tifs[layer_cnt*num_tiles:\
                                        (layer_cnt+1)*num_tiles]

        this_job_file_info.mean_normalization = mean_normalization

        job_file_path = job_info_base_path + 'l_' + str(layer_cnt)\
                        + '.jobfile'

        job_script_path = job_info_base_path + 'l_' + str(layer_cnt)\
                        + '.sh'


        job_log_path = job_info_base_path + 'l_' + str(layer_cnt)\
                        + '.log'

        job_err_path = job_info_base_path + 'l_' + str(layer_cnt)\
                        + '.err'

        with open(job_file_path, 'wb') as job_file:
            cPickle.dump(this_job_file_info, job_file)

        # creating call script
        with open(job_script_path, 'w') as call_script:
            call_script.write('#!/bin/bash\n')
            call_script.write("{0} {1} {2}".format(python_exe_path,
                                               single_render_job_script_path,
                                               job_file_path))

        # set executable bit
        os.chmod(job_script_path, 0744)

        print("Submitting job for layer {0}".format(layer_cnt))

        # run qsub - shell = True is important for the proper executation
        # environment
        subprocess.call("qsub -q {0} -o {1} -e {2} {3} {4}".format(
            sge_queue,
            job_log_path,
            job_err_path,
            sge_additional_flags,
            job_script_path), shell=True)


    print('All jobs submitted.')

    return

def render_tiles_in_layer(affines = [],
                          global_xy_layer_dimensions=(),
                          tile_paths=[],
                          out_path='',
                          out_format='tif',
                          crop_px = 25,
                          scale = 1.0,
                          invert_images = True,
                          mean_normalization = True,
                          hist_stretch = False,
                          hist_stretch_percentile = 0.5,
                          clahe = True,
                          clahe_clip_lim = 1.5,
                          clahe_size = 10,
                          gaussian_sigma = 0.5):
    """
    Reads tif files, applies given affine transformations and combines the
    tiles into a single layer file. Saves either raw or tif. Be aware of the
    2GB or 4GB tif limitations for libtiff and some tif libraries (checkout
    BigTiff, but I didn't find a python wrapper).

    :param affines: list of list of float
        affine transformation, eg read in by parse_affines_from_bills_x_a_txt
    :param global_xy_layer_dimensions: tuple of int
    :param tile_paths:
    :param out_path:
    :param out_format:
    :return:
    """
    tot_ref_time = time.time()

    # create a numpy array that holds the final results
    final_layer = np.zeros(global_xy_layer_dimensions, dtype=np.uint8)

    for tile_path, affine in zip(tile_paths, affines):
        print("loading: {0}".format(tile_path))
        # create a temp numpy array that is then merged into final_layer
        #tile_layer = np.zeros(global_xy_layer_dimensions, dtype=np.uint8)

        # fast image loader with properly buffered reading, replace with
        # anything that can read images into numpy array

        #ref_time = time.time()
        this_tile = su.load_image(tile_path, out_format='numpy')

        if this_tile.dtype != np.uint8:
            raise Exception('Data type not supported, uint8 required')

        # prepare affine, for some reason the data type gets messed when not
        # explicitly specified
        this_affine = np.array(affine).reshape([2,3]).astype(np.float32)


        # perform cropping on the image tile
        if crop_px > 0:
            this_tile = this_tile[crop_px:this_tile.shape[0]-crop_px,
                                  crop_px:this_tile.shape[1]-crop_px]

        # the opencv affine implementation performs reasonably well (compared
        # to other libraries...)
        cv2.warpAffine(this_tile, this_affine, global_xy_layer_dimensions,
                       dst=final_layer, flags=cv2.INTER_NEAREST,
                       borderMode=cv2.BORDER_TRANSPARENT)


    if gaussian_sigma > 0.0:
        # hard coded kernel size here, should change with sigma actually..:(
        # gaussian sigmas from 0.5 to 0.6 turned out to be good for denoising
        # EM data
        cv2.GaussianBlur(final_layer, (5,5), gaussian_sigma, final_layer)

    if hist_stretch == True:
        # calculating the percentile on the whole image is slow, we sample a
        # few times and take the median
        lower_percentile = hist_stretch_percentile
        upper_percentile = 100. - hist_stretch_percentile

        p1 = []
        p2 = []
        max_start = np.min(final_layer.shape) - 101

        for cnt in range(0, 1000):
            start = random.randint(0, max_start)
            end = start + 100
            this_p1, this_p2 = np.percentile(final_layer[start:end, start:end],
                                        (lower_percentile, upper_percentile))
            p1.append(this_p1)
            p2.append(this_p2)

        p1 = np.median(np.array(p1))
        p2 = np.median(np.array(p2))

        scale_factor = np.array((255.-0.)/(p2-p1)+0.)

        # basic opencv operations are crazy fast compared to numpy and
        # everything else - skimage is slowest, as always (over 20x slower..)
        cv2.subtract(final_layer, p1, final_layer)
        cv2.multiply(final_layer, scale_factor, final_layer)

    if clahe == True:
        # again, only the opencv clahe implementation performed well enough
        # to handle large datasets, 900 MB can be clahe-enhanced in about 30s
        # on a single core
        cv_clahe = cv2.createCLAHE(clahe_clip_lim,(clahe_size,clahe_size))
        final_layer = cv_clahe.apply(final_layer)

    if mean_normalization:
        # it is assumed that no more information is below a gray value of 4
        # the mask itself is necessary to calculate the mean only to the
        # actual images, and not the surrounding black parts.
        ret, thres_mask = cv2.threshold(final_layer, 4, 255, cv2.THRESH_BINARY)

        # we use final_layer as src AND as mask, because only non-zero
        # elements should be taken into account for the mean calculation
        mean = cv2.mean(final_layer, mask=thres_mask)[0]

        # again, we use final_layer as a mask as well, to not change the
        # values of the unused background
        final_layer = cv2.add(final_layer, 127 - mean, mask=thres_mask)

    if invert_images:

        # using a mask for the inversion is possible, but dangerous: you do
        # not want to exclude some voxels accidentally from the inversion...
        # this is less problematic for the normalization, were, in the worst
        # case, only a few pixels become shift a bit.
        #if not mean_normalization:
        #    ret, thres_mask = cv2.threshold(final_layer,
        #                                    3, 255, cv2.THRESH_BINARY)
        cv2.subtract(255, final_layer, final_layer)


    if scale != 1.0:
        final_layer = cv2.resize(final_layer, (0,0), fx=scale, fy=scale)

    if out_format == 'tif':
        #cv2.imwrite(out_path, final_layer)
        try:
            final_tif = Image.fromarray(final_layer)
            final_tif.save(out_path)
        except:
            raise Exception('Tif cannot be written, too big? Try raw instead.')

    elif out_format == 'raw':
        final_layer.tofile(out_path)

    print("Finished, took: {0}s".format(time.time()-tot_ref_time))

    return


def parse_affines_from_bills_x_a_txt(path_to_X_A_TXT_folder, num_tiles):
    """
    Parses affine transformations and returns them.

    path_to_X_A_TXT_folder: string
        Full path to the folder generated by Bill's xview

    num_tiles: int
        Number of tiles per layer

    returns: list of list of floats, list of paths of all affine files
        Affine transformations of all layers of all tiles

    """

    all_layer_files = [file for file in os.listdir(path_to_X_A_TXT_folder)
                       if file.lower().endswith('.txt')]

    digits = []
    for file_name in all_layer_files:
        file_digit = ''.join(x for x in file_name if x.isdigit())
        digits.append(int(file_digit))

    all_layer_files_proper_order = [x for (y,x) in sorted(zip(digits,
                                    all_layer_files),
                                    key=lambda pair: pair[0])]

    affines = []

    cnt = 0.0
    total_cnt = len(all_layer_files_proper_order)
    for layer_file in all_layer_files_proper_order:

        print('{0:3.2f}% done'.format(cnt / total_cnt * 100.))
        cnt += 1.

        with open(path_to_X_A_TXT_folder + layer_file) as f:
            tile_lines = f.readlines()

        this_layer_affines = []

        tile_cnt = 0

        # extract affines for each tile
        for tile_line in tile_lines:
            tile_cnt += 1
            # each affine is specified in column 3 to 8
            this_tile_affines = tile_line.rstrip().split('\t')[2:8]
            this_layer_affines.append(this_tile_affines)


        # Bill's aligner leaves out tiles that had problems
        if tile_cnt < num_tiles:
            print("Warning: only {0} tiles found in layer {1}, "
                  "this certainly requires further action.".format(
                    tile_cnt, layer_file))
            raise Exception('Incomplete affine specifications.')

        affines.append(this_layer_affines)

    return affines

def remove_slow_drifts_in_affines(affines, global_shift=(0,0)):
    """
    Function that removes const drifts of the x and y parameters of the
    given affine transformations (affines should come from
    parse_affines_from_bills_x_a_txt).

    :param affines:
    :return: affines
    """

    # copied from public domain:
    # http://stackoverflow.com/questions/11686720/
    # is-there-a-numpy-builtin-to-reject-outliers-from-a-list
    def reject_outliers(data, m = 20.):
        d = np.abs(data - np.median(data))
        mdev = np.median(d)
        s = d/mdev if mdev else 0
        return data[s<m]

    num_tiles = len(affines[0])

    # Extract x,y affine params from all tiles and all layers for fit
    np_affines = np.array(affines).astype(np.float32)

    x_affines = []
    y_affines = []

    # Calc slope
    for this_tile in range(0, num_tiles):
        x_affines.append(np.mean(reject_outliers(np.diff(np_affines[:,
                                                           this_tile,2]))))
        y_affines.append(np.mean(reject_outliers(np.diff(np_affines[:,
                                                         this_tile,5]))))

    x_correction = -np.mean(x_affines)
    y_correction = -np.mean(y_affines)

    x_correction_line = np.tile(np.arange(0, len(affines)) * x_correction,
                                (num_tiles,1))
    y_correction_line = np.tile(np.arange(0, len(affines)) * y_correction,
                                (num_tiles,1))

    x_correction_line = x_correction_line + global_shift[0]
    y_correction_line = y_correction_line + global_shift[1]

    np_affines[:,:,2] = np_affines[:,:,2] + np.swapaxes(x_correction_line,0,1)
    np_affines[:,:,5] = np_affines[:,:,5] + np.swapaxes(y_correction_line,0,1)

    print("const x correction applied: {0}".format(x_correction))
    print("const y correction applied: {0}".format(y_correction))

    print("global shift x applied: {0}".format(global_shift[0]))
    print("global shift y applied: {0}".format(global_shift[1]))

    #mplot.figure()
    #mplot.plot(y_correction_line[0,:])

    #mplot.figure()
    #mplot.plot(np_affines[:,0,2])

    #mplot.figure()
    #mplot.plot(np_affines[:,0,5])

    affines = np_affines.tolist()

    return affines


def plot_affines(affines, tile):
    """
    Plots the course of the affine values over z (layers). This is handy to
    diagnose alignment problems (sudden jumps, etc).

    :param affines: list of list of floats
    :param tile: int
        Specify the linear tile number
    :return:
    """

    # didn't get array slicing working here
    plt_affines = np.zeros([len(affines), 6])
    for i in range(0, len(affines)-1):
        plt_affines[i] = np.array(affines[i][tile])


    mplot.figure()
    mplot.title('tile # %d' % tile)
    mplot.subplot(6,1,1)
    mplot.plot(np.arange(0,len(affines)).T, plt_affines[:,0])

    mplot.subplot(6,1,2)
    mplot.plot(np.arange(0,len(affines)).T, plt_affines[:,1])

    mplot.subplot(6,1,3)
    mplot.plot(np.arange(0,len(affines)).T, plt_affines[:,2])

    mplot.subplot(6,1,4)
    mplot.plot(np.arange(0,len(affines)).T, plt_affines[:,3])

    mplot.subplot(6,1,5)
    mplot.plot(np.arange(0,len(affines)).T, plt_affines[:,4])

    mplot.subplot(6,1,6)
    mplot.plot(np.arange(0,len(affines)).T, plt_affines[:,5])
    mplot.show()

    return

def gen_bills_aligner_j0251_stack_import_file():
    base_file_name = '/big-data/kornfeld/rawStacks/j0251/j0251/j0251_'

    tile_coord_map = {}

    tile_coord_map[0] = [63, 29]
    tile_coord_map[1] = [6491, 0]
    tile_coord_map[2] = [12934, 18]
    tile_coord_map[3] = [19372, 45]
    tile_coord_map[4] = [19319, 6465]
    tile_coord_map[5] = [12881, 6455]
    tile_coord_map[6] = [6442, 6452]
    tile_coord_map[7] = [15, 6495]
    tile_coord_map[8] = [7, 12968]
    tile_coord_map[9] = [6436, 12908]
    tile_coord_map[10] = [12874, 12898]
    tile_coord_map[11] = [19313, 12896]
    tile_coord_map[12] = [19305, 19330]
    tile_coord_map[13] = [12865, 19349]
    tile_coord_map[14] = [6426, 19378]
    tile_coord_map[15] = [0, 19451]

    slices_z = range(0, 15494)

    tiles = range(1, 17)

    with open('/big-data/kornfeld/rawStacks/j0251/bills_aligner_j0251_full_layout.txt', 'w') as text_file:

        for curr_slice in slices_z:
            for curr_tile in tiles:
                text_file.write(str(curr_slice - slices_z[0]) + '\t' +
                                str(curr_tile) + '\t' +
                                str(1.0) + '\t' +
                                str(0.0) + '\t' +
                                str(tile_coord_map[curr_tile-1][0]) + '\t' +
                                str(0.0) + '\t' +
                                str(1.0) + '\t' +
                                str(tile_coord_map[curr_tile-1][1]) + '\t' +
                                str(-999) + '\t' +
                                str(-999) + '\t' +
                                str(0) + '\t' +
                                base_file_name +
                                "{0:06d}".format(curr_slice) + '_' +
                                "{0:03d}".format(curr_tile) + '.tif ' +
                                '\n')


    return

def gen_bills_aligner_stack_import_file():
    base_file_name = '/nobackup/hackathon/data/denklab/songbird_stack/songbird_stack_'
    ovlap = 150 # use about half of the overlap value in pixels
    tile_len = 6600


    # only rectangular tiling patterns supported currently, in "snake" order
    slices_z = range(0, 1)
    tiles_x = range(0, 4)
    tiles_y = range(0, 4)

    # tile_coord_map is a hash: tile_num -> xy upper-left coord of image
    tile_coord_map = dict()

    # generate tile coord map
    linear_tile_cnt = 0
    for curr_y in tiles_y:
        for curr_x in tiles_x:
            if curr_y%2:
                tile_coord_map[linear_tile_cnt]  = \
                    [(max(tiles_x) - curr_x) * (tile_len - ovlap),
                     curr_y * (tile_len - ovlap)]
            else:
                tile_coord_map[linear_tile_cnt]  = \
                    [curr_x * (tile_len - ovlap),
                     curr_y * (tile_len - ovlap)]

            linear_tile_cnt += 1

    tiles = range(1, linear_tile_cnt+1)

    with open("/home/jk/single_layer_coords.txt", "w") as text_file:

        for curr_slice in slices_z:
            for curr_tile in tiles:
                text_file.write(str(curr_slice - slices_z[0]) + '\t' +
                                str(curr_tile) + '\t' +
                                str(1.0) + '\t' +
                                str(0.0) + '\t' +
                                str(tile_coord_map[curr_tile-1][0]) + '\t' +
                                str(0.0) + '\t' +
                                str(1.0) + '\t' +
                                str(tile_coord_map[curr_tile-1][1]) + '\t' +
                                str(-99) + '\t' +
                                str(-99) + '\t' +
                                str(0.0) + '\t' +
                                base_file_name +
                                "{0:06d}".format(curr_slice) + '_' +
                                "{0:03d}".format(curr_tile) + '.tif ' +
                                '\n')

    return


if __name__ == '__main__':
    single_render_job(sys.argv[1])
	