# CNVMetrics 0.1.6

NEW FEATURES

* `calculateLog2ratioMetric()` method enables log2 ratio metric calculation using similar workflow than `calculateOverlapRegionsMetric()` method.

SIGNIFICANT USER-VISIBLE CHANGES

* `plotMetric()` replaces `plotOverlapMetric()` method. The new method can plot all metrics (state call metrics and log2 ratio metrics).
* Vignette section 'Workflow for metrics calculated using the level of amplification/deletion' is complete.
* New citing section in README and vignette refering to published F1000Research poster (http://www.doi.org/10.7490/f1000research.1118704.1). 
* Instead of calculating distance, log2 ratio metrics are calculated distance-based metrics (1/(1+distance)).

BUG FIXES

* None


# CNVMetrics 0.1.4

NEW FEATURES

* None

SIGNIFICANT USER-VISIBLE CHANGES

* New website https://krasnitzlab.github.io/CNVMetrics/index.html associated to package.
* Vignette section 'Workflow for metrics calculated using CNV status calls' is complete.

BUG FIXES

* None


# CNVMetrics 0.1.3

NEW FEATURES

* None

SIGNIFICANT USER-VISIBLE CHANGES

* `plotOneOverlapMetric()` method has a new argument `silent=TRUE` so that the plot is not drawn by default.

BUG FIXES

* `plotOneOverlapMetric()` method now uses sample distance for clustering as default clustering method.


# CNVMetrics 0.1.2

NEW FEATURES

* `plotOneOverlapMetric()` method enables plotting result of overlapping metric calculation.

SIGNIFICANT USER-VISIBLE CHANGES

* `calculateOverlapRegionsMetric()` method changed to `calculateOverlapMetric()`.

BUG FIXES

* None


# CNVMetrics 0.1.1

NEW FEATURES

* Added a `NEWS.md` file to track changes to the package.
* `calculateOverlapRegionsMetric()` method enables calculation of similarity metrics using overlapping amplified/deleted regions.

SIGNIFICANT USER-VISIBLE CHANGES

* None.

BUG FIXES

* None
