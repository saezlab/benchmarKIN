Installation
============

To install ``benchmarKIN`` please run the following lines.


.. code-block:: console

   # Install specific versions of dependencies (currently has to be done separate)
   remotes::install_version("Matrix", version="1.6-5", repos="http://cran.us.r-project.org") 
   remotes::install_version("MASS", version="7.3-60", repos="http://cran.us.r-project.org")

   # Install the package from GitHub
   remotes::install_github("saezlab/benchmarKIN")

