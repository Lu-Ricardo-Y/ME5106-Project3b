# ME5106 Project 3b — Acoustic monitoring

MATLAB analysis of the laser additive manufacturing recordings. The repository includes every generated analysis script and result from steps 1–6, including editable MATLAB figures and intermediate data.

| Step | Included results |
| --- | --- |
| 1 · Inspect recording | [Full and local plots, envelopes](me3_step1_results/) |
| 2 · Locate acoustic activity | [Spectra, high-band envelopes, track close-ups](me3_step2_results/) |
| 3 · Segment | [245 candidate tracks and 49-group overview](me3_step3_results/) |
| 4 · Check boundaries and labels | [Boundary checks, label mapping, duration comparison](me3_step4_results/) |
| 5 · Extract features | [Track and group features, controlled comparisons](me3_step5_results/) |
| 6 · Validate and predict | [Full error and prediction tables](me3_step6_results/validation_and_predictions.txt) and [MATLAB results](me3_step6_results/validation_and_predictions.mat) |

## Reproduce the analysis

Open the repository folder in MATLAB. Step 6 can be rerun directly with the included step 5 results. To regenerate all results, obtain the course-provided recording and place it at `data/Mini Project 3 Data.mat`. The recording and course PDFs are distributed separately through the course.

```matlab
me3_step1
me3_step2
me3_step2_zoom
me3_step3_diagnostic
me3_step3_export
me3_step4_qa_labels
me3_step5_features
me3_step6_validate_predict
```

The recording is sampled at 160 kHz. Segmentation used the high-frequency RMS envelope and was checked across four thresholds and six plotted groups. The five tracks per sample were aggregated before model validation. Step 6 compares the median baseline, duration-only ridge and acoustic-feature ridge, using 46 fully labeled groups and both sample-group and parameter-combination holdouts.

| Sample | Estimated power (W) | Scan speed (mm/s) |
| --- | ---: | ---: |
| 41 | 236.73 | 708.12 (estimated) |
| 43 | 269.30 | 592.37 (estimated) |
| 49 | 270.61 | 850 (given) |

For the acoustic-feature model, leave-one-combination-out MAE is 50.19 W and 69.34 mm/s. Speed improves over the median baseline (81.52 mm/s), while power does not (30.22 W baseline). Treat the power outputs as exploratory. The assignment table marks sample 49 power unknown, while a figure caption lists 425 W. We excluded its power from training. One of its five input features falls outside the training range.
