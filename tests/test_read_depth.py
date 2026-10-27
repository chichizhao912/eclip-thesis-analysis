from __future__ import annotations

import numpy as np
import pysam

from eclip_thesis_analysis.read_depth import simulate


def test_read_thinning_reproduces_full_depth_before_sampling(tmp_path):
    profile = np.zeros((1, 300), dtype=np.float32)
    profile[0, 100:150] = 1
    dataset = tmp_path / "windows.npz"
    np.savez_compressed(
        dataset, ip1=profile, ip2=profile, sminput=profile,
        sequence=np.zeros((1, 4, 300), dtype=np.float32),
        window_id=np.array(["w1"]), chrom=np.array(["chr1"]),
        start=np.array([0]), end=np.array([300]), strand=np.array(["+"]),
        split=np.array(["test"]), window_class=np.array(["peak"]),
    )
    paths = {}
    for track in ("ip1", "ip2", "sminput"):
        path = tmp_path / f"{track}.bam"
        with pysam.AlignmentFile(str(path), "wb",
                                 header={"HD": {"VN": "1.6", "SO": "coordinate"},
                                         "SQ": [{"SN": "chr1", "LN": 1000}]}) as bam:
            read = pysam.AlignedSegment(bam.header)
            read.query_name = f"{track}_read"
            read.flag = 128
            read.reference_id = 0
            read.reference_start = 100
            read.mapping_quality = 60
            read.cigar = ((0, 50),)
            read.query_sequence = "A" * 50
            read.query_qualities = pysam.qualitystring_to_array("I" * 50)
            bam.write(read)
        pysam.index(str(path))
        paths[track] = path
    output = tmp_path / "thinning"
    simulate(dataset, paths, output, fractions=(1.0, 0.5), samples=1, seed=1)
    with np.load(output / "depth_1.00_sample1.npz", allow_pickle=False) as full:
        assert np.array_equal(full["ip1"], profile)
    with np.load(output / "depth_0.50_sample1.npz", allow_pickle=False) as thin:
        assert set(np.unique(thin["ip1"])) <= {0, 2}
    assert (output / "read_coverage_audit.json").exists()

