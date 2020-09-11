ecmapper
========

Just a test to see if it is possible to extract map EC number from an organism to BIGG and ModelSEED reactions.

Requirements
~~~~~~~~~~~~

It requires `biopython <https://github.com/biopython/biopython>`__ and `cobrapy <https://github.com/opencobra/cobrapy>`__.

Installation
~~~~~~~~~~~~

In the ecmapper git folder:

.. code:: sh
    python3 setup.py develop

This tool takes as input a genbank file (containing the EC number). It also uses database files from BIGG and ModelSEED.

To download the database into a folder use the command:

.. code:: sh
    ecmapper database -d database_folder

To use ecmapper:

.. code:: sh
    ecmapper map -g genbank_file -d database_folder -o output_folder
