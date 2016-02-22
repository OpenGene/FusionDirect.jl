"""
IMPORTANT_FUSIONS is used to whitelist those important fusions to not be filtered by read support or quality value
For example, fusions with only one read support will be discarded, but if the two genes are listed below, it will be kept
This feature is useful for detecting low frequency fusions in cancer sequencing applications
"""
# "*" means any
const IMPORTANT_FUSIONS = [
    (("ALK", "intron", "19"), ("EML4", "intron", "*")),
    (("ALK", "intron", "19"), ("KIF5B", "intron", "*")),
    (("RET", "intron", "11"), ("CCDC6", "intron", "*")),
    (("RET", "intron", "*"), ("KIF5B", "intron", "*")),
    (("ROS1", "intron", "33"), ("CD74", "intron", "*")),
    (("ROS1", "intron", "33"), ("EZR", "intron", "*")),
    (("ROS1", "intron", "*"), ("GOPC", "intron", "*")),
    (("NTRK1", "intron", "*"), ("TPM3", "intron", "*")),
]