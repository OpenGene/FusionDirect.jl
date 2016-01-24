"""
IMPORTANT_FUSIONS is used to whitelist those important fusions to not be filtered by read support or quality value
For example, fusions with only one read support will be discarded, but if the two genes are listed below, it will be kept
This feature is useful for detecting low frequency fusions in cancer sequencing applications
"""
const IMPORTANT_FUSIONS = [
    ("ALK", "EML4"),
]