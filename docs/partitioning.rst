Partitioning
============

.. image:: ../images/partitioning.png
  :width: 200

Graph partitioners such as METIS, partitions each component mesh independently. Therefore, a partition does not contain more than one mesh. Also, METIS with default settings produces partitions that contain approximately the same number of cells. In overset mesh technique, overlapping mesh-cells need to be in the processors. In order to bring overlapping mesh-cells to the same partitions, geometric partitioning is applied. An octree is used to registed mesh-cells based on their spatial locations.

Load balancing
--------------

For balanced distrubution of work-load, the octree is refined adaptively at the most-loaded octree-bins until balanced distribution of octree-bins to processors is possible. At each refinement step, octree is converted to a graph and passed to METIS to find optimum distribution.
