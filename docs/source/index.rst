benchmarKIN - Evaluation of methods for kinase activity inference
===================================

**benchmarKIN** is a package for evaluation of kinase activity inference tools. It
provides two complementary benchmarking approaches for evaluation: The perturbation-based
benchmark approach and the tumor-based benchmark approach.

.. figure:: graphical_abstract.png
   :height: 300px
   :alt: benchmarKIN's workflow
   :align: center
   :class: no-scaled-link

   benchmarKIN contains two complementary approaches for kinase activity inference


Check out the `perturbation-based <https://benchmarkin.readthedocs.io/en/latest/notebooks/perturbBench.html>`_ or `tumor-based <https://benchmarkin.readthedocs.io/en/latest/notebooks/tumorBench.html>`_ benchmark sections for further information,
of how test your own kinase activity inference method.

If you have any question or problem do not hesitate to open an `issue <https://github.com/saezlab/benchmarKIN/issues>`_.

Citation
--------
Mueller-Dott, Sophia, Eric J. Jaehnig, Khoi Pham Munchic, Wen Jiang, Tomer M. Yaron-Barir, Sara R. Savage, Martin Garrido-Rodriguez, et al. 2024. “Comprehensive Evaluation of Phosphoproteomic-Based Kinase Activity Inference.” bioRxiv. https://doi.org/10.1101/2024.06.27.601117.

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Main

   installation

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Vignettes

   notebooks/perturbBench
   notebooks/tumorBench
