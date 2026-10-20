# BiteOscope_human_pathogen

Repository for the study "Deep behavioral phenotyping of pathogen infected
mosquitoes reveals species-specific behavior changes enhancing transmission"

The repository holds the DeepLabCut-assisted pose estimation, the per-mosquito
behaviour statistics and the transmission model. The behavioural classification
tool itself is a separate repository: https://github.com/RUMCvectorBiology/

## Layout

    BiteOscope_analysis/          pose estimation and per-mosquito statistics,
                                  see BiteOscope_analysis/README.md
    MBITES/                       transmission model, see MBITES/README.md

## Pipeline

    video
      -> DeepLabCut_New_project.py           project, labels, training
      -> DeepLabCut_Analyze_video.py         poses and tracklets per video
      -> DeepLabCut_Extracting_ind_track.py  one h5 per individual
      -> Deeplabcut_engorge_trackstat.ipynb  one row per individual
      -> Statistics_BiteOscope.ipynb         landing / engorgement statistics
      -> Statistics_BiteOserve.ipynb         behaviour-classification statistics

The behaviour sequence (walk, probe, engorge) comes from the external
classification tool; Statistics_BiteOserve.ipynb reads its output.

## Running

Python 3 with deeplabcut (multi-animal), tensorflow, opencv, pandas, scipy,
statsmodels and the two helper packages `biteOscope_utils` and
`biteOserver_utils`, which are not part of this repository. The scripts and
notebooks carry the absolute paths of the Linux analysis machine
(`/mnt/DATA2/zhong/...`); set them to your own before running.

## Settings

  - Multi-animal DeepLabCut, dlcrnet_ms5, two shuffles, 80,000 iterations,
    crop size 1000 x 1000 px; legs, mouthparts and antennae are labelled but
    ignored when the tracklets are stitched.
  - 25 frames per second and 0.04 mm per pixel: in the statistics notebooks
    frame counts are divided by 25 and pixel distances multiplied by 0.04.
  - Engorgement thresholds and the track statistics are listed in
    BiteOscope_analysis/README.md.
  - The transmission-model code and its settings are described in
    MBITES/README.md.
