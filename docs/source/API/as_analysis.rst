Alternative Splicing Analysis
===================================

DOLPHIN performs alternative splicing analysis using the `Outrigger module <https://github.com/YeoLab/outrigger>`_ from the `Expedition toolkit <https://github.com/YeoLab/Expedition>`_. This method quantifies alternative splicing events from aggregated BAM files and computes percent spliced-in (PSI) values for exon-exon junctions [1]_.

.. [1] Song, Y., Botvinnik, O. B., Lovci, M. T., Kakaradov, B., Liu, P., Xu, J. L., & Yeo, G. W. (2017). *Single-cell alternative splicing analysis with Expedition reveals splicing dynamics during neuron differentiation*. Molecular Cell, 67(1), 148–161. https://pubmed.ncbi.nlm.nih.gov/28673540/

.. autofunction:: DOLPHIN.alternative_splicing.convert_psi_to_h5ad.run_convert_psi
.. autofunction:: DOLPHIN.alternative_splicing.convert_random_psi.run_psi_random
.. autofunction:: DOLPHIN.alternative_splicing.generate_differential_as.run_differential_as

The historical ``convert_modality_ohe`` helper remains available from the
``legacy-v1`` branch. DOLPHIN v2 performs route and modality handling through
the unified alternative-splicing pipeline instead.
