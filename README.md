# SigHotSpotter

## About

A key goal of regenerative medicine is to restore structure and function of damaged tissues and organs due to disease and ageing. Despite recent advances, the identification of key factors for efficient control of cell phenotype to enable cellular rejuvenation is still a challenge. 
Several computational methods have been useful in designing cellular conversion strategies, however, they mainly rely on gene regulatory network (GRN) models without accounting for the cellular microenvironment (niche), which is important for maintenance of in vivo cellular phenotype.
Our main assumption here is that a cellular phenotype is maintained by the sustained effect of the niche via key signaling molecules. Accordingly, we propose a probabilistic method `SigHotSpotter` to identify such signaling molecules by integrating signaling and transcriptional networks.
Application of SigHotSpotter on several niche dependent cellular systems correctly predicted known niche induced signaling molecules. Importantly, SigHotSpotter outperforms commonly employed signaling pathway enrichment/inference methods.

### Authors

This software was developed in the [Computational Biology Group](https://wwwfr.uni.lu/lcsb/research/computational_biology) by
- [Srikanth Ravichandran](https://wwwen.uni.lu/lcsb/people/srikanth_ravichandran)
- [Andras Hartmann](https://wwwfr.uni.lu/lcsb/people/andras_hartmann)
- [Antonio del Sol](https://wwwfr.uni.lu/lcsb/people/antonio_del_sol_mesa)

## System requirements 
### Hardware requirements
There are no special hardware requirements, the tool was tested on a virtual machine with 4GB of memory, a 4 core Intel Xeon 2.5 GHz, and 50 GB allocated space on a hard disk. 

### Software dependencies
The package supports all main operating systems, and requires `R` version 3.3 or later. Tested configurations:
- Windows: Windows 7 / R 3.5.1
- linux: Ubuntu 16.04 TLS / R 3.4.4
- macOS: Sierra 10.12.6 / R 3.5.3

The web interface depends on shiny webserver (deployed on version 1.5.9.923), all further software dependencies are contained in the package description. In order to use the software as a standalone application rstudio is recommended (tested on Version 1.0.143) .

## Installation guide (linux/unix)
From an `RStudio` terminal, type: 

```R
require("devtools")
install_git("https://gitlab.com/srikanth.ravichandran/sighotspotter")
```

This will install all the package dependencies as well as the SigHotSpotter package, it takes about 3 minutes.
If the install would fail because the dependencies can not be installed, you can install the dependencies before installing the package by typing:
```R
install.packages( c("plyr", "igraph", "Matrix", "reshape2", "RSpectra", "dplyr", "snow", "shiny", "DT", "shinyjs", "shinythemes", "shinyBS", "rintrojs", "openxlsx", "markdown") )
```
Takes about 3-5 minutes to install.


## Instruction to use

### Running the interface locally

From an `RStudio` terminal, type: 

```R
shiny::runApp(paste(find.package("SigHotSpotter"), "webapp", "app.R", sep = .Platform$file.sep))
```

This will open the graphical user interface locally

### Web interface
An online web interface of the SigHotSpotter application (SOFTWARE for short) is accessible at
[https://sighotspotter.lcsb.uni.lu](https://sighotspotter.lcsb.uni.lu), and is provided free of charge for academic, non-profit use.
Using the SOFTWARE means you accept the terms and conditions of the Disclaimer below.

### Application walkthrough

Click on the "Take a tour" button for walkthrough.

### Input file formats

The required input of SigHotSpotter are tab separated value files. 
- The files for condition 1 and condition 2 shall contain gene expression data (normalized bulk or single-cell measurements).
Rows represent genes labeled with gene symbols, and columns represent replicates or data for single cells.
- The differential expression file shall contain only the significantly differentially expressed genes where 1 denotes up-regulated genes in condition 1 and -1 denotes up-regulated genes in condition 2
See an example dataset
<a href="https://webdav-r3lab.uni.lu/public/cbg/SigHotSpotter/data/SigHotSpotter_datasets.zip" target="_blank">here</a>.

## Demo

Download and unzip the test dataset (e.g. in your tmp folder) from <a href="https://webdav-r3lab.uni.lu/public/cbg/SigHotSpotter/data/SigHotSpotter_datasets.zip" target="_blank">here</a>, or use the following `shell` commands to download the dataset:

```bash
cd /tmp
wget https://webdav-r3lab.uni.lu/public/cbg/SigHotSpotter/data/SigHotSpotter_datasets.zip
unzip SigHotSpotter_datasets.zip
```
Access the application as specified in the previous section, and upload the test data.
For example, to test the mESC predictions, choose the following files:
- condition 1 = 2i.txt
- condition 2 = lif.txt
- differential expression = DE_BOOLEAN_2i_lif.txt

Parameters:
- species = Mouse
- cutoff = 50
- percentile = 90

## Time estimates for running

Calculation of hotspots depends heavily on the settings and the input data type, and ranges from 30 seconds to 10 minutes. In case of cutoff set to 30% on single-cell data, the calculation takes less than two minutes on the configuration described above.


## Disclaimer

THE SOFTWARE IS NOT TO BE USED FOR TREATING OR DIAGNOSING HUMAN SUBJECTS.

THE SOFTWARE AND ALL CONTENT ON THIS SERVER ARE PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
