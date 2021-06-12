Force update
************

Warning: you will lose any of your local changes to the repository   

.. code-block:: bash
   
   cd star-forming-regions
   git fetch --all
   git reset --hard origin/master
   git submodule update --force --remote -- lime
   python setup.py install
