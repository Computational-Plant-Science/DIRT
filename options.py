from pathlib import Path


class TraitOptions:
    def __init__(
            self,
            crown: bool = True,
            medial: bool = True,
            rtp: bool = True,
            rdist: bool = True,
            hypocotol: bool = True,
            lateral: bool = True,
            sta_quantile: bool = True,
            rta: bool = True,
            sta: bool = True,
            cp_diameter: bool = True,
            sta_dominant: bool = True,
            sta_dominant_25: bool = True,
            sta_dominant_50: bool = True,
            sta_dominant_75: bool = True,
            sta_dominant_90: bool = True,
            rta_dominant: bool = True,
            rtp_skeleton: bool = True,
            branching_frequency: bool = True,
            rtp_lateral: bool = True,
            lateral_length: bool = True,
            lateral_angles: bool = True):
        self.crown = crown
        self.medial = medial
        self.rtp = rtp
        self.rdist = rdist
        self.hypocotol = hypocotol
        self.lateral = lateral
        self.sta_quantile = sta_quantile
        self.rta = rta
        self.sta = sta
        self.cp_diameter = cp_diameter
        self.sta_dominant = sta_dominant
        self.sta_dominant_25 = sta_dominant_25
        self.sta_dominant_50 = sta_dominant_50
        self.sta_dominant_75 = sta_dominant_75
        self.sta_dominant_90 = sta_dominant_90
        self.rta_dominant = rta_dominant
        self.rtp_skeleton = rtp_skeleton
        self.branching_frequency = branching_frequency
        self.rtp_lateral = rtp_lateral
        self.lateral_length = lateral_length
        self.lateral_angles = lateral_angles

    @staticmethod
    def from_dict(**entries):
        options = TraitOptions()
        options.__dict__.update(entries)
        return options


class DIRTOptions:
    def __init__(
            self,
            input_file,
            output_directory: str = '',
            marker_diameter: float = 25.4):
        self.input_file = input_file
        self.input_name = Path(input_file).name
        self.input_stem = Path(input_file).stem
        self.output_directory = output_directory
        self.marker_diameter = marker_diameter
