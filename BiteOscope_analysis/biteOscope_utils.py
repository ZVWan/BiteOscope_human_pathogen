'''
author: Felix Hol
date: 2022  Sept
content: create test videos of single DLC generated tracks (creates one video per tracklet). Crop to trajectory
'''

import cv2
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter
import glob
import pandas as pd
import numpy as np
import warnings
from tqdm.notebook import tqdm
import os
import re
warnings.filterwarnings("ignore")

def save_individual_tracklets(pkl, min_length, **kwargs):
    saveDir = kwargs.get('saveDir', None)
    allTracklets = pd.read_pickle(pkl)
    header = allTracklets.pop('header')
    scorer = header.get_level_values("scorer").unique().to_list()
    bpts = header.get_level_values("bodyparts").unique().to_list()
    coords = header.get_level_values("coords").unique().to_list()
    animal_names = ['ind1']

    columns = pd.MultiIndex.from_product(
                [scorer, animal_names, bpts, coords],
                names=["scorer", "individuals", "bodyparts", "coords"],
            )

    get_frame_ind = lambda s: int(re.findall(r"\d+", s)[0])

    for k, tracklet in allTracklets.items():
        if len(tracklet) > min_length - 1:
            inds, temp = zip(*[(get_frame_ind(k), v) for k, v in tracklet.items()])
            data = np.asarray(temp, dtype=np.float16)
            flatData = data[:, :, 0:3].reshape((np.shape(data)[0], -1))
            df = pd.DataFrame(flatData, columns=columns, index=inds)
            if saveDir:
                output_name = saveDir + os.path.split(pkl)[1][:-7] + '_' + str('%05d' % k) + '.h5'
            else:
                output_name = pkl[:-7] + '_' + str('%05d' % k) + '.h5'
            df.to_hdf(output_name, "df_with_missing", format="table", mode="w")



def createVideoDLCTrackCrop(hdf, videoIn, videoOut, minDuration, maxDuration, everyNthFrame, pCutOff):
    cropExtra = 5          ### the full frame is cropped to the extreme values of the body parts + cropExtra ##Previously 35
    df = pd.read_hdf(hdf)
    scorer = df.columns.get_level_values(0)[0]
    for ind in df.columns.get_level_values(1).unique():
        for bp in df.columns.get_level_values(2).unique():
            df[scorer, ind, bp, 'x'].mask(df[scorer, ind, bp, 'likelihood'] < pCutOff,inplace=True)
            df[scorer, ind, bp, 'y'].mask(df[scorer, ind, bp, 'likelihood'] < pCutOff,inplace=True)
    if len(df) > minDuration - 1:
        scorer = df.columns.get_level_values(0)[0]
        bptsToPlot = ['proboscis1', 'proboscis2', 'head', 'thorax', 'abdomenC', 'abdomenR', 'abdomenL', 'bottom']
            #    'rightForeLeg1', 'rightForeLeg2', 'rightForeLeg3',
            #    'leftForeLeg1', 'leftForeLeg2', 'leftForeLeg3',
            #    'rightMidleg1', 'rightMidleg2', 'rightMidleg3',
            #    'leftMidleg1', 'leftMidleg2', 'leftMidleg3',
            #    'rightHindleg1', 'rightHindleg2', 'rightHindleg3',
            #    'leftHindleg1', 'leftHindleg2', 'leftHindleg3', 'leftMPalp','rightMPalp','antennaL','antennaR']
        cmap = mpl.cm.get_cmap('plasma', len(bptsToPlot) + 2)
        start = df.first_valid_index() - 1
        stop = df.last_valid_index()
        if stop - start > maxDuration:
            stop = start + maxDuration
        if start < 1:
            start = 1

        cap = cv2.VideoCapture(videoIn)
        fps = cap.get(cv2.CAP_PROP_FPS)
        cap.set(1, start)
        ret, frame = cap.read()
        (h, w) = frame.shape[:2]
        dpi = 100
        writer = FFMpegWriter(fps=fps, codec='h264')

        ind = df.columns.get_level_values(1).unique()[0]
        x_min = int(df[scorer, ind].loc[:,(slice(None),'x')].min().min())
        x_max = int(df[scorer, ind].loc[:,(slice(None),'x')].max().max())
        y_min = int(df[scorer, ind].loc[:,(slice(None),'y')].min().min())
        y_max = int(df[scorer, ind].loc[:,(slice(None),'y')].max().max())

        if x_min < cropExtra + 1:
            x_min = 0
        else:
            x_min = x_min - cropExtra
        if y_min < cropExtra + 1:
            y_min = 0
        else:
            y_min = y_min - cropExtra
        if x_max < w - cropExtra + 1:
            x_max = x_max + cropExtra
        else:
            x_max = w
        if y_max < h - cropExtra + 1:
            y_max = y_max + cropExtra
        else:
            y_max = h

        fig = plt.figure(frameon=False, figsize=((x_max - x_min) / dpi, (y_max - y_min) / dpi))
        plt.set_cmap('plasma')

        with writer.saving(fig, videoOut, dpi=dpi):
            for i in tqdm(range(start, stop, everyNthFrame)):
                cap.set(1, i)
                ret, frame = cap.read()
                plt.imshow(frame[y_min:y_max, x_min:x_max])
                j = 0
                try:
                    for ind in df.columns.get_level_values(1).unique():
                        for bp in bptsToPlot:
                            plt.scatter(df[scorer, ind, bp].x[i] - x_min, df[scorer, ind, bp].y[i] - y_min, color=cmap(j), alpha=0.5)
                            j += 1
                except:
                    print('no data for frame ' + str(i))
                plt.axis('off')
                plt.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
                writer.grab_frame()
                plt.clf()

        plt.close(fig)



def createVideoDLCTrackCrop1(hdf, videoIn, videoOut, minDuration, maxDuration, everyNthFrame, pCutOff):
    cropExtra = 5          ### the full frame is cropped to the extreme values of the body parts + cropExtra ##Previously 35
    df = dfAll
    scorer = df.columns.get_level_values(0)[0]
    for ind in df.columns.get_level_values(1).unique():
        for bp in df.columns.get_level_values(2).unique():
            df[scorer, ind, bp, 'x'].mask(df[scorer, ind, bp, 'likelihood'] < pCutOff,inplace=True)
            df[scorer, ind, bp, 'y'].mask(df[scorer, ind, bp, 'likelihood'] < pCutOff,inplace=True)
    if len(df) > minDuration - 1:
        scorer = df.columns.get_level_values(0)[0]
        bptsToPlot = ['proboscis1', 'proboscis2', 'head', 'thorax', 'abdomenC', 'abdomenR', 'abdomenL', 'bottom',
               'rightForeLeg1', 'rightForeLeg2', 'rightForeLeg3',
               'leftForeLeg1', 'leftForeLeg2', 'leftForeLeg3',
               'rightMidleg1', 'rightMidleg2', 'rightMidleg3',
               'leftMidleg1', 'leftMidleg2', 'leftMidleg3',
               'rightHindleg1', 'rightHindleg2', 'rightHindleg3',
               'leftHindleg1', 'leftHindleg2', 'leftHindleg3', 'leftMPalp','rightMPalp','antennaL','antennaR']
        cmap = mpl.cm.get_cmap('plasma', len(bptsToPlot) + 2)
        start = df.first_valid_index() - 12 
        stop = df.last_valid_index()
        if stop - start > maxDuration:
            stop = start + maxDuration
        if start < 1:
            start = 1

        cap = cv2.VideoCapture(videoIn)
        fps = cap.get(cv2.CAP_PROP_FPS)
        cap.set(1, start)
        ret, frame = cap.read()
        (h, w) = frame.shape[:2]
        dpi = 100
        writer = FFMpegWriter(fps=fps, codec='h264')

        ind = df.columns.get_level_values(1).unique()[0]
        x_min = int(df[scorer, ind].loc[:,(slice(None),'x')].min().min())
        x_max = int(df[scorer, ind].loc[:,(slice(None),'x')].max().max())
        y_min = int(df[scorer, ind].loc[:,(slice(None),'y')].min().min())
        y_max = int(df[scorer, ind].loc[:,(slice(None),'y')].max().max())

        if x_min < cropExtra + 1:
            x_min = 0
        else:
            x_min = x_min - cropExtra
        if y_min < cropExtra + 1:
            y_min = 0
        else:
            y_min = y_min - cropExtra
        if x_max < w - cropExtra + 1:
            x_max = x_max + cropExtra
        else:
            x_max = w
        if y_max < h - cropExtra + 1:
            y_max = y_max + cropExtra
        else:
            y_max = h

        fig = plt.figure(frameon=False, figsize=((x_max - x_min) / dpi, (y_max - y_min) / dpi))
        plt.set_cmap('plasma')

        with writer.saving(fig, videoOut, dpi=dpi):
            for i in tqdm(range(start, stop, everyNthFrame)):
                cap.set(1, i)
                ret, frame = cap.read()
                plt.imshow(frame[y_min:y_max, x_min:x_max])
                j = 0
                try:
                    for ind in df.columns.get_level_values(1).unique():
                        for bp in bptsToPlot:
                            plt.scatter(df[scorer, ind, bp].x[i] - x_min, df[scorer, ind, bp].y[i] - y_min, color=cmap(j), alpha=0.5)
                            j += 1
                except:
                    print('no data for frame ' + str(i))
                plt.axis('off')
                plt.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
                writer.grab_frame()
                plt.clf()

        plt.close(fig)
