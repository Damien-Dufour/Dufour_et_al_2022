Workflow_PLA
================
Damien, Dufour

# Description

How to analyse data from PLA on tissues

# Protocol

Adrenals were fixed in PFA 4% for 6 hours and embedded in paraffin. 5µm
sections were used according to Duolink PLA protocol. Pictures were
taken with a Zeiss Imager M2 with enough stacks to cover all the visible
*foci*.

# Image analysis

## Merge stacks into one based on maximal values

run this script in FIJI. It makes you choose an image, then merge the
stacks based on max values, save in .tif and close the picture.

``` ijm
 path = File.openDialog("Select a File");
  //open(path); // open the file
  dir = File.getParent(path);
  name = File.getName(path);
 open(path); 
 run("Z Project...", "projection=[Max Intensity]");
 saveAs(".tif", path);
 close("*");
```

## Detect cells and *foci* in QuPath

The first plug-in detects the *nuclei* based on DAPI (threshold = 100)
and expend cytoplasm to 10 µm (default = 5µm).

The second one detects subcellular *foci* of high intensity (400) (Can
be specific to *nuclei* by adding a threshold on DAPI channel) between
0.1 and 3 µm<sup>2</sup>

``` qupath
setImageType('FLUORESCENCE');
selectAnnotations(); 
runPlugin('qupath.imagej.detect.cells.WatershedCellDetection', '{"detectionImage": "Channel 1", "requestedPixelSizeMicrons": 0.5, "backgroundRadiusMicrons": 8.0, "medianRadiusMicrons": 0.0, "sigmaMicrons": 1.5, "minAreaMicrons": 10.0, "maxAreaMicrons": 400.0, "threshold": 100.0, "watershedPostProcess": true, "cellExpansionMicrons": 10.0, "includeNuclei": true, "smoothBoundaries": true, "makeMeasurements": true}'); 
runPlugin('qupath.imagej.detect.cells.SubcellularDetection', '{"detection[Channel 1]": -1.0, "detection[Channel 2]": -1.0, "detection[Channel 3]": 400.0, "doSmoothing": false, "splitByIntensity": false, "splitByShape": false, "spotSizeMicrons": 1.0, "minSpotSizeMicrons": 0.1, "maxSpotSizeMicrons": 3.0, "includeClusters": true}');
```

## Clean data and plot

``` r
# open the file and only keep the interesting columns (Image refers to the full name of the picture, Parent is the name of the annotation (here refers to genotype) and the last one is the number of foci per cells)

df <- read.csv("measurements.csv")[,c(1,4,66)] #the number of the third column to keep depends on the channel where PLA foci are

# rename the columns

names(df) <- c("Image", 
               "Genotype",
               "Foci_nb")

# sometimes, QuPath yields estimations of the nb of foci that are decimal hence they have to be rounded 

df$Foci_nb <- round(df$Foci_nb)

# put in the right order if necessary

df$Foci_nb <- factor(df$Foci_nb,
                     levels = c("genotype1",
                                "genotype2"))

# plot 

ggplot(df)+
  aes(x = Genotype,
      fill = Foci_nb)+
  geom_bar(position = "fill")+
  scale_fill_viridis_d()+
  theme_classic()

# compare

chisq.test(df$Foci_nb,
           df$Genotype)
```
