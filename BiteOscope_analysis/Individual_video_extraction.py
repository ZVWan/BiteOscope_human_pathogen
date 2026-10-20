import biteOscope_utils

hdf = '/individual_hdf.h5'
videoIn = '/original_video.mp4'
videoOut = '/output_individual_tracking_video.mp4'
minDuration = 25
maxDuration = 4500
everyNthFrame = 3
pCutOff = 0.8
biteOscope_utils.createVideoDLCTrackCrop(hdf, videoIn, videoOut, minDuration, maxDuration, everyNthFrame, pCutOff)