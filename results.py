from typing import TypedDict


class DIRTResults(TypedDict, total=False):
    image_name: str
    seconds: int
    tag: str
    marker_ratio: float
    marker_width: float
    marker_height: float
    x_scale: float
    y_scale: float
    branching_points: int
    depth: float
    width: float
    stem_diameter: float
    tip_median_diameter: float
    tip_mean_diameter: float
    paths: int
    path_seg_1: float
    path_seg_2: float
    projected_area: float
    mean_density: float
    median_width: float
    max_width: float
    acc_width_10: float
    acc_width_20: float
    acc_width_30: float
    acc_width_40: float
    acc_width_50: float
    acc_width_60: float
    acc_width_70: float
    acc_width_80: float
    acc_width_90: float
    acc_width_slope_10: float
    acc_width_slope_20: float
    acc_width_slope_30: float
    acc_width_slope_40: float
    acc_width_slope_50: float
    acc_width_slope_60: float
    acc_width_slope_70: float
    acc_width_slope_80: float
    acc_width_slope_90: float
    max_diameter_90: float
    drop_50: float
    spatial_dist_x: float
    spatial_dist_y: float
    top_angle: float
    top_dominant_angle_1: float
    top_dominant_angle_2: float
    bottom_angle: float
    min_angle: float
    max_angle: float
    median_angle: float
    soil_tissue_angle_range: float
    soil_tissue_min_angle: float
    soil_tissue_max_angle: float
    soil_tissue_median_angle: float
    soil_tissue_dominant_angle_1: float
    soil_tissue_dominant_angle_2: float
    soil_tissue_dominant_angle_1_25: float
    soil_tissue_dominant_angle_1_50: float
    soil_tissue_dominant_angle_1_75: float
    soil_tissue_dominant_angle_1_90: float
    soil_tissue_dominant_angle_2_25: float
    soil_tissue_dominant_angle_2_50: float
    soil_tissue_dominant_angle_2_75: float
    soil_tissue_dominant_angle_2_90: float
    adventitious: int
    adventitious_angles: float
    basal: int
    basal_angles: float
    hypocotol_diameter: float
    taproot_diameter: float
    taproot_diameter_25: float
    taproot_diameter_50: float
    taproot_diameter_75: float
    taproot_diameter_90: float
    nodal_path_length: float
    nodal_mean_diameter: float
    lateral_mean_length: float
    lateral_branching_frequency: float
    lateral_mean_angle: float
    lateral_angular_range: float
    lateral_min_angle: float
    lateral_max_angle: float
    lateral_distance_to_first: float
    lateral_median_diameter: float
    lateral_mean_diameter: float







