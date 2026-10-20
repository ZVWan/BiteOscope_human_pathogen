##### Package import
import os
from pathlib import Path
import glob
import tensorflow as tf
import deeplabcut
import psutil
import time
import multiprocessing

os.environ["DLClight"]="True"

##### Restricting use of GPU
gpus = tf.config.experimental.list_physical_devices('GPU')
if gpus:
  try:
    for gpu in gpus:
      tf.config.experimental.set_memory_growth(gpu, True)
  except RuntimeError as e:
    print(e)

##### Restricting use of CPU, FOR CONVERTING TRACKLETS
multiprocessing.set_start_method('spawn')

##### Setting path
newvideo1 = glob.glob('/Your/Path/To/Your/Video/*.mp4')
config_path = '/anopheles-zw-2023-07-02/config.yaml' 

##### Analyze videos
deeplabcut.analyze_videos(config_path, newvideo1, shuffle=1, batchsize=16, auto_track=False, dynamic=(True, 0.5, 100))
deeplabcut.create_video_with_all_detections(config_path, newvideo1)

##### Extract Outliers
pickle = glob.glob('/Your/Path/To/Your/Video/*.pickle')
print(pickle)
deeplabcut.find_outliers_in_raw_data(config_path, pickle[0], newvideo1[0], extraction_algo='uniform', pcutoff=0.1, with_annotations=True)

#### Generate & stitch tracklets
video_tracking = glob.glob('/Your/Path/To/Your/Video/*DLC_resnet50_AnophelesJul2023shuffle1_1030000*mp4*')
print(video_tracking, sep = "\n")
deeplabcut.convert_detections2tracklets(config_path, newvideo1, identity_only=False, ignore_bodyparts = [
       'rightForeLeg1', 'rightForeLeg2', 'rightForeLeg3', 'rightForeLeg4',
       'leftForeLeg1', 'leftForeLeg2', 'leftForeLeg3', 'leftForeLeg4',
       'rightMidleg1', 'rightMidleg2', 'rightMidleg3', 'rightMidleg4',
       'leftMidleg1', 'leftMidleg2', 'leftMidleg3', 'leftMidleg4',
       'rightHindleg1', 'rightHindleg2', 'rightHindleg3', 'rightHindleg4',
       'leftHindleg1', 'leftHindleg2', 'leftHindleg3', 'leftHindleg4',
       'leftMPalp', 'rightMPalp','antennaR','antennaL','proboscis1','proboscis2']  ,overwrite=True)

deeplabcut.stitch_tracklets(config_path, newvideo1, shuffle=1, track_method='ellipse', max_gap=None, min_length=25, n_tracks=25, prestitch_residuals=True)
deeplabcut.create_labeled_video(config_path, ['*.mp4'], track_method = "box", color_by = 'individual')

##### Extract outliers
deeplabcut.extract_outlier_frames(config_path, video_tracking, shuffle=1, outlieralgorithm='uncertain', track_method='ellipse', p_bound=0.5, comparisonbodyparts=['rightForeLeg1','leftForeLeg1'], extractionalgorithm='uniform', copy_videos=False, automatic=True)

##### Merge for new training session
# deeplabcut.merge_datasets(config_path)
# deeplabcut.convertcsv2h5(config_path, userfeedback=True)

