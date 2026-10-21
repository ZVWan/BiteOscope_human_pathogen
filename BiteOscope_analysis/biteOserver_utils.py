'''
author: zw
date: 2024 oct
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
import os
import pandas as pd

def process_behaviour_file(behaviour_label, output_dir, confidence_threshold, min_duration):
    behaviour_label_process = pd.read_csv(behaviour_label)
    
    if behaviour_label_process.empty:
        consecutive_behaviour_df = pd.DataFrame(columns=['start_frame', 'Duration', 'Behaviour'])
        output_file_name = f"(blank)_{os.path.splitext(os.path.basename(behaviour_label))[0]}_consecutive_behaviour.csv"
    else:
        filtered_behaviour = behaviour_label_process[behaviour_label_process['Confidence'] >= confidence_threshold].reset_index(drop=True)
        
        if len(filtered_behaviour) == 0:
            consecutive_behaviour_df = pd.DataFrame(columns=['start_frame', 'Duration', 'Behaviour'])
            output_file_name = f"(blank)_{os.path.splitext(os.path.basename(behaviour_label))[0]}_consecutive_behaviour.csv"
        else:
            consecutive_behaviour = []
            count = 1
            start_frame = 0
            current_behaviour = filtered_behaviour.loc[0, 'Behaviour']

            for i in range(1, len(filtered_behaviour)):
                if filtered_behaviour.loc[i, 'Behaviour'] == filtered_behaviour.loc[i-1, 'Behaviour']:
                    count += 1
                else:
                    if count > min_duration:
                        consecutive_behaviour.append((start_frame, count, current_behaviour))
                    start_frame = i
                    count = 1
                    current_behaviour = filtered_behaviour.loc[i, 'Behaviour']
            
            if count > min_duration:
                consecutive_behaviour.append((start_frame, count, current_behaviour))
            
            consecutive_behaviour_df = pd.DataFrame(consecutive_behaviour, columns=['start_frame', 'Duration', 'Behaviour'])
            output_file_name = f"{os.path.splitext(os.path.basename(behaviour_label))[0]}_consecutive_behaviour.csv"
    
    output_csv_path = os.path.join(output_dir, output_file_name)
    os.makedirs(os.path.dirname(output_csv_path), exist_ok=True)
    consecutive_behaviour_df.to_csv(output_csv_path, index=False)

def process_behaviour_sequence(behaviour_label, output_dir):
    df = pd.read_csv(behaviour_label)

    if df.empty:
        result_df = pd.DataFrame(columns=['start_frame', 'Duration', 'Behaviour'])
        new_file_name = f"(blank)_{os.path.splitext(os.path.basename(behaviour_label))[0]}_sequence.csv"
    else:
        result = []
        if len(df) > 0:
            current_start = df.loc[0, 'start_frame']
            current_duration = df.loc[0, 'Duration']
            current_behaviour = df.loc[0, 'Behaviour']

            for i in range(1, len(df)):
                if df.loc[i, 'Behaviour'] == current_behaviour:
                    current_duration = df.loc[i, 'start_frame'] + df.loc[i, 'Duration'] - current_start
                else:
                    result.append((current_start, current_duration, current_behaviour))
                    current_start = df.loc[i, 'start_frame']
                    current_duration = df.loc[i, 'Duration']
                    current_behaviour = df.loc[i, 'Behaviour']

            result.append((current_start, current_duration, current_behaviour))

        result_df = pd.DataFrame(result, columns=['Start_frame', 'Duration', 'Behaviour'])
        new_file_name = f"{os.path.splitext(os.path.basename(behaviour_label))[0]}_sequence.csv"
    
    output_file_path = os.path.join(output_dir, new_file_name)
    os.makedirs(output_dir, exist_ok=True)
    
    result_df.to_csv(output_file_path, index=False)


def behaviour_frequency (input_dir, output_csv):
    csv_files = glob.glob(os.path.join(input_dir, '*.csv'))
    all_data = []

    for file in csv_files:
        df = pd.read_csv(file)

        behaviour_summary = df.groupby('Behaviour').agg(
            count = ('Behaviour','size'),
            duration = ('Duration', 'sum')
        ).reset_index()

        # behaviour_counts = df['Behaviour'].value_counts().reset_index()
        # behaviour_counts.columns = ['Behaviour', 'Count']
        total_counts = behaviour_summary['count'].sum()
        behaviour_summary['Frequency'] = behaviour_summary['count'] / total_counts

        # duration_sum = df.groupby('Behaviour')['Duration'].sum().reset_index()
        total_duration = behaviour_summary['duration'].sum()
        behaviour_summary['Duration_ratio'] = behaviour_summary['duration']/total_duration

        # result = pd.merge(behaviour_counts, duration_sum, on = 'Behaviour')

        behaviour_summary['File'] = os.path.basename(file)
        all_data.append(behaviour_summary)
    
    result_df = pd.concat(all_data, ignore_index=True)
    result_df.to_csv(output_csv, index=False)
    return result_df

def aggregate_behaviour_duration(input_folder):
    csv_files = glob.glob(os.path.join(input_folder, '*.csv'))
    
    all_data = []
    
    for file in csv_files:
        df = pd.read_csv(file)
        behaviour_duration_sum = df.groupby('Behaviour')['Duration'].sum().reset_index()
        behaviour_duration_sum['File'] = os.path.basename(file)
        all_data.append(behaviour_duration_sum)
    
    result_df = pd.concat(all_data, ignore_index=True)
    
    return result_df