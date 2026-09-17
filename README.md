# WSP
Accompanying repository of the paper "Wildfire Suppression: Complexity, Models, and Instances".

## Data

**Instances**. The instances are in `instances`. This directory contains two zipped files, one for each experimental section (literature and instance generator). For a comprehensive description of the instance format, refer to the [instance generator README](scripts/README.md).

**Experiments**. Directory `data` contains two folders, one for each experimental section (literature or instance generator). Each folder contains the experimental data and additional information about the instances, such as the best known value, the best known lower bound, and the number of vertices. The generated-instance results include one row per algorithm run and explicit experimental-factor columns. The scripts and instructions to reproduce the generated-instance Skillings--Mack analyses in Table 5 and Figure 3 are in [`analysis`](analysis/README.md).

## Instance Generator
The instance generator is implemented as a python script, which can be found in the directory `scripts`. The same directory also contains detailed instructions on how to use the generator.
