from biotk.libs.poa.ConsensusTensor import SubreadSetCircularConsensusTensors
from biotk.libs.QuickKinetics import kinetics
from pbcore.io import ReferenceSet

class TestSubreadSetCircularConsensusTensors:
    def __init__(self, sset_path,
                       ref=None,
                       n_zmws_sample=3,
                       min_coverage_depth=3,
                       coverage_depth=None,
                       tensors_context_width=1,
                       tensors_collection_mode='equal-state',
                       tensors_per_poa=5,
                       kinetics_summary=None):

        self.sset_path = sset_path
        self.ref = ref
        self.n_zmws_sample = n_zmws_sample
        self.min_coverage_depth = min_coverage_depth
        self.coverage_depth = coverage_depth
        self.tensors_context_width = tensors_context_width
        self.tensors_collection_mode = tensors_collection_mode
        self.tensors_per_poa = tensors_per_poa
        self.kinetics_summary = kinetics_summary
        self.sset_ccs_tensors = SubreadSetCircularConsensusTensors(
            self.sset_path,
            ref=self.ref,
            n_zmws_sample=self.n_zmws_sample,
            min_coverage_depth=self.min_coverage_depth,
            coverage_depth=self.coverage_depth,
            tensors_context_width=self.tensors_context_width,
            tensors_collection_mode=self.tensors_collection_mode,
            tensors_per_poa=self.tensors_per_poa,
            kinetics_summary=self.kinetics_summary
        )

def test_subreadset_circular_consensus_tensors():
    sset_path = 'data/tiny_set_internal.subreadset.xml'
    kin_summary = kinetics(sset_path=sset_path,
                           nreads=10,
                           samples_per_read=100,
                           unique_zmws=True).summarizeKinetics()

    tscct = TestSubreadSetCircularConsensusTensors(
        sset_path=sset_path,
        ref=ReferenceSet('data/references/All4mers_InsertOnly.ReferenceSet.xml')[0],
        n_zmws_sample=3,
        min_coverage_depth=3,
        coverage_depth=5,
        tensors_context_width=1,
        tensors_collection_mode='standard',
        tensors_per_poa=None,
        kinetics_summary=kin_summary
    )
    print 'hello'
