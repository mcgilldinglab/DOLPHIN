Alternative Splicing Analysis
===================================

DOLPHIN performs alternative splicing analysis using the `Outrigger module <https://github.com/YeoLab/outrigger>`_ from the `Expedition toolkit <https://github.com/YeoLab/Expedition>`_. This method quantifies alternative splicing events from aggregated BAM files and computes percent spliced-in (PSI) values for exon-exon junctions [1]_.

.. [1] Song, Y., Botvinnik, O. B., Lovci, M. T., Kakaradov, B., Liu, P., Xu, J. L., & Yeo, G. W. (2017). *Single-cell alternative splicing analysis with Expedition reveals splicing dynamics during neuron differentiation*. Molecular Cell, 67(1), 148–161. https://pubmed.ncbi.nlm.nih.gov/28673540/

.. autofunction:: DOLPHIN.as_analysis.convert_psi_to_h5ad.convert_psi_to_h5ad
.. autofunction:: DOLPHIN.as_analysis.differential_as_data.prepare_data_for_differential_AS_analysis
