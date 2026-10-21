# BiteOscope_analysis

DeepLabCut pipeline and downstream statistics for the BiteOscope videos of the
manuscript.

## Layout

    DeepLabCut_New_project.py           new DeepLabCut project and training
    DeepLabCut_Analyze_video.py         poses and tracklets per video
    DeepLabCut_Extracting_ind_track.py  individual h5 files
    Deeplabcut_engorge_trackstat.ipynb  per-individual track statistics
    Statistics_BiteOscope.ipynb         landing / engorgement statistics
    Statistics_BiteOserve.ipynb         behavior-classification statistics
    Individual_video_extraction.py      crop video per individual
    test_clip_2.mp4                     short clip to test the pipeline
    biteOscope_utils.py                 helper functions for the track steps
    biteOserver_utils.py                helper functions for the behavior statistics

## Video pipeline

    video
      -> DeepLabCut_New_project.py           create new DeepLabCut project for training
      -> DeepLabCut_Analyze_video.py         generating poses and tracklets per video
      -> DeepLabCut_Extracting_ind_track.py  extract individual h5 files
      -> Deeplabcut_engorge_trackstat.ipynb  extract morphological and kinetic information for each individual
      -> Statistics_BiteOscope.ipynb         landing / engorgement statistics
      -> Statistics_BiteOserve.ipynb         behavior-classification and probing statistics

The behavioral sequence (walk, probe, engorge) comes from the external
classification tool; Statistics_BiteOserve.ipynb reads its output.

## Running

Python 3 with deeplabcut (multi-animal, maDLC), tensorflow, opencv, pandas, scipy,
statsmodels and the two helper packages `biteOscope_utils` and
`biteOserver_utils`, which sit in this folder. The scripts and
notebooks carry the absolute paths of the Linux analysis machine
; set them to your own before running.

    python DeepLabCut_New_project.py          # once per project, then label
    python DeepLabCut_Analyze_video.py        # per batch of videos
    python DeepLabCut_Extracting_ind_track.py # per pickle
    python Individual_video_extraction.py     # optional, per individual
    jupyter notebook Deeplabcut_engorge_trackstat.ipynb
    jupyter notebook Statistics_BiteOscope.ipynb
    jupyter notebook Statistics_BiteOserve.ipynb

## Settings

Parameters and thresholds are set at the top of each script and in the first
cells of the notebooks.
