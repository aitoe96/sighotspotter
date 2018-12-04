# NicheSIG

## About

A key goal of regenerative medicine is to restore structure and function of damaged tissues and organs due to disease and ageing. Despite recent advances, the identification of key factors for efficient control of cell phenotype to enable cellular rejuvenation is still a challenge. 
Several computational methods have been useful in designing cellular conversion strategies, however, they mainly rely on gene regulatory network (GRN) models without accounting for the cellular microenvironment (niche), which is important for maintenance of in vivo cellular phenotype.
Our main assumption here is that a cellular phenotype is maintained by the sustained effect of the niche via key signaling molecules. Accordingly, we propose a probabilistic method `NicheSIG` to identify such signaling molecules by integrating signaling and transcriptional networks.
Application of NicheSIG on several niche dependent cellular systems correctly predicted known niche induced signaling molecules. Importantly, NicheSIG outperforms commonly employed signaling pathway enrichment/inference methods.

### Authors

This software was developed in the [Computational Biology Group] (https://wwwfr.uni.lu/lcsb/research/computational_biology) by
- [Srikanth Ravichandran](https://wwwen.uni.lu/lcsb/people/srikanth_ravichandran)
- [Andras Hartmann](https://wwwfr.uni.lu/lcsb/people/andras_hartmann)
- [Antonio del Sol](https://wwwfr.uni.lu/lcsb/people/antonio_del_sol_mesa)

## System requirements 
### Hardware requirements
The `NicheSIG` package requires only a standard PC, for optimal performance, at least 8 GB of RAM and 4 cores of 2.5 GHz each (or above) are recommended.

### Software dependencies
The package supports linux and macOS operating systems, and requires `R` version 3.2.3 or later. Tested platforms:
- linux: Ubuntu 16.04 TLS / R 3.2.3
- macOS: Sierra 10.12.6 / R 3.5.1

Running the interface locally requires `RStudio`, we suggest to use the latest build (RStudio Desktop 1.1.463)

## Installation guide (linux/unix)
From an `RStudio` terminal, type: 

```R
require("devtools")
install_git("https://gitlab.uni.lu/sravichandran/nichesig.git")
```

This will install all the package dependencies as well as the NicheSIG package, it takes about 3 minutes.



## Instruction to use

### Running the interface locally

From an `RStudio` terminal, type: 

```R
shiny::runApp(paste(find.package("NicheSIG"), "webapp", "app.R", sep = .Platform$file.sep))
```

This will open the graphical user interface locally

### Web interface
An online web interface of the application is accessible at
[https://nichesig.lcsb.uni.lu](https://nichesig.lcsb.uni.lu), and is provided free of charge for academic, non-profit use.


### Application walkthrough

Click on the "Take a tour" button for walkthrough.

### Input file formats

The required input of NicheSIG are tab separated value files. 
- The files for condition 1 and condition 2 shall contain gene expression data (normalized bulk or single-cell measurements).
Rows represent genes labeled with gene symbols, and columns represent replicates or data for single cells.
- The differential expression file shall contain only the significantly differentially expressed genes where 1 denotes up-regulated genes in condition 1 and -1 denotes up-regulated genes in condition 2
See an example dataset
<a href="https://webdav-r3lab.uni.lu/public/cbg/NicheSIG/data/NicheSIG_datasets.zip" target="_blank">here</a>.

## Demo

Download and unzip the test dataset (e.g. in your tmp folder) from <a href="https://webdav-r3lab.uni.lu/public/cbg/NicheSIG/data/NicheSIG_datasets.zip" target="_blank">here</a>, or use the following `shell` commands to download the dataset:

```bash
cd /tmp
wget https://webdav-r3lab.uni.lu/public/cbg/NicheSIG/source/NicheSIG.zip
unzip NicheSIG.zip
```
Access the application as specified in the previous section, and upload the test data.
For example, to test the mESC predictions, choose the following files:
- condition 1 = ESC_data.txt
- condition 2 = Epi_LC_D2_data.txt
- differential expression = DEG_ESC_Epi_LC_D2.txt

Parameters:
- species = Mouse
- cutoff = 100
- percentile = 90

The expected runtime is 6 minutes, see the expected output in the corresponding sheet of the Supplementary Table 1. to Computational approach to control cell phenotype for cellular rejuvenation strategies

---
The [NicheSIG](https://nichesig.lcsb.uni.lu/webapp/) web-based tool (SOFTWARE for short) is provided free of charge for academic, non-profit use.
For commercial use, please contact the authors for a license.
Using the SOFTWARE means you accept the terms and conditions of the Disclaimer below.

### Disclaimer

THE SOFTWARE IS NOT TO BE USED FOR TREATING OR DIAGNOSING HUMAN SUBJECTS.

THE SOFTWARE AND ALL CONTENT ON THIS SERVER ARE PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
