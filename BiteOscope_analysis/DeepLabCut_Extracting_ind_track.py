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
import biteOscope_utils

pkl = glob.glob('/mnt/DATA2/zhong/BiteOScope_Data/20250401_proline/*el.pickle')

min_length = 25

for pkls in pkl:
    biteOscope_utils.save_individual_tracklets(pkls, min_length)
