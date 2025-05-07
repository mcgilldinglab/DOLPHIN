DOLPHIN Preprocess Module
=========================

This module provides functions for processing GTF files and generating non-overlapping exon annotations.

Main Function
-------------

This is the recommended entry point for users.

.. autofunction:: DOLPHIN.preprocess.generate_exon_gtf.generate_nonoverlapping_exons
.. autofunction:: DOLPHIN.preprocess.generate_adj_index.generate_adj_index_table

Helper Functions
----------------

The following functions support different steps of the GTF processing pipeline. Advanced users may call these directly.

.. autofunction:: DOLPHIN.preprocess.generate_exon_gtf.prepare_exon_gtf

.. autofunction:: DOLPHIN.preprocess.generate_exon_gtf.exon_uniq

.. autofunction:: DOLPHIN.preprocess.generate_exon_gtf.save_by_batch

.. autofunction:: DOLPHIN.preprocess.generate_exon_gtf.combine_saved_batches

.. autofunction:: DOLPHIN.preprocess.generate_exon_gtf.check_exon_overlap

.. autofunction:: DOLPHIN.preprocess.generate_exon_gtf.save_gtf_outputs