##### Package import
import os
from pathlib import Path
import glob
import tensorflow as tf
import deeplabcut

gpus = tf.config.experimental.list_physical_devices('GPU')
if gpus:
  try:
    for gpu in gpus:
      tf.config.experimental.set_memory_growth(gpu, True)
  except RuntimeError as e:
    print(e)

##### Creating New project
video_path = glob.glob('/Your/Path/To/Your/Video/*.mp4')
print(video_path, sep = "\n")
config_path = deeplabcut.create_new_project('anopheles','zw', video_path, copy_videos=False, multianimal=True) 
print(config_path, sep = "\n")


##### Add videos & Extract frames for labelling
deeplabcut.add_new_videos(config_path, video_path, copy_videos=False)
deeplabcut.extract_frames(config_path, mode='automatic', algo='uniform', userfeedback=False, crop=True)
deeplabcut.label_frames(config_path)

##### Create new training session
deeplabcut.merge_datasets(config_path)
deeplabcut.create_multianimaltraining_dataset(config_path, net_type='dlcrnet_ms5', num_shuffles=2, crop_size=(1000,1000))
deeplabcut.train_network(config_path, allow_growth=True, displayiters=500, saveiters=2500, shuffle=1, maxiters=80000)
deeplabcut.evaluate_network(config_path, Shuffles=[1,2], plotting=True, gputouse=0)
deeplabcut.extract_save_all_maps(config_path, shuffle=shuffle, Indices=[0, 1, 2, 3])

##### Refine training
deeplabcut.refine_labels(config_path)

##### Analyze videos
deeplabcut.analyze_videos(config_path, videos, shuffle=1, batchsize=16, auto_track=False, dynamic=(True, 0.5, 100))

##### Extract Outliers
pickle = glob.glob('Your/Path/To/Your/Video/*.pickle')
print(pickle)
deeplabcut.find_outliers_in_raw_data(config_path, pickle[0], videos[0], extraction_algo='uniform', pcutoff=0.1, with_annotations=True)
