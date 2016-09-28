Purpose:
===

This is the distribution of the ClusPath code. It is distributed under the GNU General Public License v3.0.
This code was developed as a research prototype and it should be treated as such. It comes with no guarantees or support, however questions can be addressed to [Marian-Andrei@rizoiu.eu](mailto:Marian-Andrei@rizoiu.eu).

Contains:
===
The **ClusPath** source code and the CPDS1 dataset. For the description of the ClusPath algorithm, please check the paper:  
[M.-A. Rizoiu, J. Velcin, S. Bonnevay and S. Lallich, "ClusPath: A Temporal-driven Clustering to Infer Typical Evolution Paths," Data Mining and Knowledge Discovery, pp. 1–26, 2015.](http://arxiv.org/pdf/1512.03501.pdf)

Requirements:
===
The current code was developed and tested under Linux using [Octave](https://www.gnu.org/software/octave/) v3.8.1, with the following packages installed:
* control
* econometrics
* general
* io
* missing-functions
* nan
* optim
* struct
* statistics

Install additional tools for plotting (in Ubuntu):
``` $> sudo apt-get install graphviz epstool ```


Octave is a fast-changing environment and new functionalities appear constantly, while old functionalities might be retired without notice.
Workarounds might be required to make ClusPath function with the newest version of Octave.

Example usage:
===

Start `octave` in the main folder of the project. Upon lauching `octave`, the ClusPath code is loaded into the Octave path (check the `.octaverc` file).
```
$> octave
```

Next, execute ClusPath with parameters: `Alpha = 0.95, Beta = 0.0002, Delta = 3, lamb1 = 1, lamb2 = 10, lamb3 = 1000`
```
octave:2> [Adj, clusassign, dtimestamp, centroids, ctimestamp] = temporalKMeans('file', 'data/cpds1.csv', 'typeT', 5, 'Alpha', 0.95, 'Beta', 0.0002, 'Delta', 3, 'lamb1', 1, 'lamb2', 10, 'lamb3', 1000, 'saveResults', true );
```

Evaluate the clustering using the four measures:
```
octave:3> evaluateTemporalClustering('inputfile', 'data/cpds1.csv', 'clusassign', clusassign, 'centroids', centroids, 'ctimestamp', ctimestamp, 'Alpha', 0.95);
```
