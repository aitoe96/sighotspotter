## About
---

A key goal of regenerative medicine is to restore structure and function of damaged tissues and organs due to disease and ageing. Despite recent advances, the identification of key factors for efficient control of cell phenotype to enable cellular rejuvenation is still a challenge. 
Several computational methods have been useful in designing cellular conversion strategies, however, they mainly rely on gene regulatory network (GRN) models without accounting for the cellular microenvironment (niche), which is important for maintenance of in vivo cellular phenotype.
Our main assumption here is that a cellular phenotype is maintained by the sustained effect of the niche via key signaling molecules. Accordingly, we propose a probabilistic method (NicheSIG) to identify such signaling molecules by integrating signaling and transcriptional networks.
Application of NicheSIG on several niche dependent cellular systems correctly predicted known niche induced signaling molecules. Importantly, NicheSIG outperforms commonly employed signaling pathway enrichment/inference methods

### Authors
---

This software was developed in the [Computational Biology Group] (https://wwwfr.uni.lu/lcsb/research/computational_biology) by
- [Srikanth Ravichandran](https://wwwen.uni.lu/lcsb/people/srikanth_ravichandran)
- [Andras Hartmann](https://wwwfr.uni.lu/lcsb/people/andras_hartmann)
- [Antonio del Sol](https://wwwfr.uni.lu/lcsb/people/antonio_del_sol_mesa)

## Documentation
---
### File formats


The required input of NicheSIG are tab separated value files. 
- The files for condition 1 and condition 2 shall contain gene expression data (normalized bulk or single-cell measurements).
Rows represent genes labeled with gene symbols, and columns represent replicates or data for single cells.
- The differential expression file shall contain only the significantly differentially expressed genes where 1 denotes up-regulated genes in condition 1 and -1 denotes up-regulated genes in condition 2
See an example dataset
<a href="https://webdav-r3lab.uni.lu/public/cbg/NicheSIG/data/NicheSIG_datasets.zip" target="_blank">here</a>.

### Application walkthrough


Click on the "Take a tour" button for walkthrough.

## Use
---

The NicheSIG web-based tool (SOFTWARE for short) is provided free of charge for academic, non-profit use.
For commercial use, please contact the authors for a license.
Using the SOFTWARE means you accept the terms and conditions of the Disclaimer below.


### Disclaimer

THE SOFTWARE IS NOT TO BE USED FOR TREATING OR DIAGNOSING HUMAN SUBJECTS.

THE SOFTWARE AND ALL CONTENT ON THIS SERVER ARE PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
