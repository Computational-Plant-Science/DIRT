from pathlib import Path


class DIRTOptions:
    def __init__(
            self,
            input_path,
            trait_file_path,
            excised_roots: int = 0,
            mask_threshold: float = 1.0,
            marker_diameter: float = 0.0,
            crown_root: bool = False,
            segmentation: bool = False,
            stem_reconstruction: bool = False,
            plot: bool = False,
            traits: list = None):
        self.input_path = input_path
        self.input_name = Path(input_path).name
        self.input_stem = Path(input_path).stem
        self.trait_file_path = trait_file_path
        self.mask_threshold = mask_threshold
        self.excised_roots = excised_roots
        self.marker_diameter = marker_diameter
        self.crown_root = crown_root
        self.segmentation = segmentation
        self.stem_reconstruction = stem_reconstruction
        self.plot = plot
        self.traits = traits
