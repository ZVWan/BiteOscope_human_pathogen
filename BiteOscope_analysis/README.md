# BiteOscope analysis: pose estimation and per-mosquito statistics

DeepLabCut pipeline and the downstream statistics for the BiteOscope videos of
the manuscript on pathogen-induced behavioural changes in mosquitoes. Every
script is one step; the steps exchange only the h5 and csv files they write.

## Layout

    DeepLabCut_New_project.py         creates the multi-animal project
                                      ('anopheles', 'zw', multianimal), extracts
                                      frames uniformly with cropping, opens the
                                      labelling GUI and trains dlcrnet_ms5
                                      (2 shuffles, crop 1000 x 1000 px,
                                      maxiters 80000)

    DeepLabCut_Analyze_video.py       analyses new videos with shuffle 1, writes
                                      a video with all detections, stitches the
                                      tracklets per individual (ellipse method,
                                      min_length 25, n_tracks 25, residual
                                      stitching on), renders a labelled video
                                      and extracts uncertain frames for
                                      relabelling (p_bound 0.5)

    DeepLabCut_Extracting_ind_track.py
                                      cuts the stitched tracks into one h5 per
                                      individual with
                                      biteOscope_utils.save_individual_tracklets

    Individual_video_extraction.py    renders a cropped video per individual from
                                      an individual h5 (biteOscope_utils.
                                      createVideoDLCTrackCrop)

    Deeplabcut_engorge_trackstat.ipynb
                                      one row per individual: distance, time,
                                      speed, moving fraction, abdominal width and
                                      engorgement score

    Statistics_BiteOscope.ipynb       statistics on that table for the landing and
                                      engorgement videos (Wilcoxon rank-sum,
                                      Fisher exact, Cliff's delta)

    Statistics_BiteOserve.ipynb       statistics on the output of the behaviour
                                      classification (walk, probe and engorge
                                      durations, probe counts) through
                                      biteOserver_utils

    test_clip_2.mp4                   short clip to run the pipeline end to end

## Running

    python DeepLabCut_New_project.py          # once per project, then label
    python DeepLabCut_Analyze_video.py        # per batch of videos
    python DeepLabCut_Extracting_ind_track.py # per pickle
    python Individual_video_extraction.py     # optional, per individual
    jupyter notebook Deeplabcut_engorge_trackstat.ipynb
    jupyter notebook Statistics_BiteOscope.ipynb
    jupyter notebook Statistics_BiteOserve.ipynb

The first two scripts need deeplabcut (multi-animal), tensorflow and opencv;
the last three need pandas, scipy, statsmodels and the helpers
`biteOscope_utils` and `biteOserver_utils`, which are not part of this
repository. The notebooks set their paths as absolute paths of the Linux
analysis machine, and the track statistics notebook ends by writing
`saveDir/'/your/name.h5'`, so both have to be edited before a run.

## Settings

  - Tracking on 25 fps video; totDist is the mean of the thorax, head and
    abdomen distances, and an individual counts as moving when the smoothed
    thorax distance exceeds 2.5 px per frame.
  - Tracks shorter than 25 frames are dropped. Smoothing windows are 3 frames
    for x and y, 5 frames for distance and 25 frames for abdominal width.
  - Abdominal width is the distance between abdomenL and abdomenR, smoothed over
    25 frames and masked outside 7-42 px. An individual is scored engorged (1)
    when the width passes 18 px, starts below 18.5 px and rises by more than
    1.2x within the first 400 frames; 2 means it was already engorged when the
    track started. timeToEngorge is in frames and is divided by 25 in the
    statistics, where pixel distances are converted with 0.04 mm per pixel.
  - Detections below a likelihood of 0.8 are dropped, x and y are interpolated
    over gaps of at most 2 frames, and outlier frames are extracted with a
    pcutoff of 0.1.
  - Statistics are Wilcoxon rank-sum tests for the group comparisons and Fisher
    exact tests for the landing and engorgement contingency tables; medians and
    means are reported per individual, with the per-individual values coming
    from one row of the track statistics table.
