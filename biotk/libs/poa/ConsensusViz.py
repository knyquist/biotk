import ConsensusTensor

class ConsensusViz(ConsensusTensor.SubreadSetCircularConsensusTensors):
    """

    """

    def __init__(self, sset_path,
                       ref=None,
                       n_zmws_sample=None,
                       min_coverage_depth=3,
                       coverage_depth=None,
                       tensors_context_width=0,
                       tensors_collection_mode='standard',
                       tensors_per_poa=None,
                       kinetics_summary=None):
        """

        """
        ConsensusTensor.SubreadSetCircularConsensusTensors.__init__(
            sset_path=sset_path,
            ref=ref,
            n_zmws_sample=n_zmws_sample,
            min_coverage_depth=min_coverage_depth,
            coverage_depth=coverage_depth,
            tensors_context_width=tensors_context_width,
            tensors_collection_mode=tensors_collection_mode,
            tensors_per_poa=tensors_per_poa,
            kinetics_summary=kinetics_summary
        )