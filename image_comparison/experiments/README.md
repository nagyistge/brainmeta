## Image Comparison

### Experiment1
This experiment pertains to testing the similarity metrics. The data has been produced for this, but the current focus is on the following experiments, for how to deal with missing data (previously known as "masking strategy")

### Experiment2
This first experiment was testing a "masking strategy," and is better termed "How to deal with missing data." We wrote a maniscript and I [gave a talk](http://www.vbmis.com/bmi/media/talks/03032015Sochat.mp4) and after substantial feedback, am working toward the same goal of assessing how to deal with missing data, but using task-based data for the base standard (see experiment 3)

### Experiment3
Is version 2 of experiment2! For this experiment I will also assess the distributions and ranking of similarity scores, however I will be using HCP task data. [CURRENTLY IN PROGRESS]

### Experiment4
We want to know if some weighting of the voxels in the map can improve the similarity calculation. We will first investigate if there is a meaningful relationship between residuals in cope images and the cope values.  If so, we can predict residuals from the randomise cope group map, and then possibly use the residuals as a weighting for the map (the other option is to edit the fsl source code to output actual residuals). We will then test different weightings of the voxels toward returning the gold standard ordering of images (eg, same contrast first). We will need to define different gold standards of ordering (tasks, contrasts groups).
