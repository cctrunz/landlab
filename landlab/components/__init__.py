from .chi_index import ChiFinder
from .depth_dependent_diffusion import DepthDependentDiffuser
from .depth_dependent_taylor_soil_creep import DepthDependentTaylorDiffuser
from .detachment_ltd_erosion import DepthSlopeProductErosion, DetachmentLtdErosion
from .diffusion import LinearDiffuser
from .drainage_density import DrainageDensity
from .erosion_deposition import ErosionDeposition
from .fire_generator import FireGenerator
from .flexure import Flexure
from .flow_accum import FlowAccumulator, LossyFlowAccumulator
from .flow_director import (
    FlowDirectorD8,
    FlowDirectorDINF,
    FlowDirectorMFD,
    FlowDirectorSteepest,
)
from .flow_routing import DepressionFinderAndRouter, FlowRouter
from .gflex import gFlex
from .lake_fill import LakeMapperBarnes
from .landslides import LandslideProbability
from .lithology import LithoLayers, Lithology
from .nonlinear_diffusion import PerronNLDiffuse
from .normal_fault import NormalFault
from .overland_flow import (
    KinwaveImplicitOverlandFlow,
    KinwaveOverlandFlowModel,
    OverlandFlow,
    OverlandFlowBates,
)
from .pet import PotentialEvapotranspiration
from .plant_competition_ca import VegCA
from .potentiality_flowrouting import PotentialityFlowRouter
from .radiation import Radiation
from .sink_fill import SinkFiller, SinkFillerBarnes
from .soil_moisture import SoilInfiltrationGreenAmpt, SoilMoisture
from .space import Space
from .spatial_precip import SpatialPrecipitationDistribution
from .steepness_index import SteepnessFinder
from .stream_power import (
    FastscapeEroder,
    SedDepEroder,
    StreamPowerEroder,
    StreamPowerSmoothThresholdEroder,
)
from .taylor_nonlinear_hillslope_flux import TaylorNonLinearDiffuser
from .transport_length_diffusion import TransportLengthHillslopeDiffuser
from .uniform_precip import PrecipitationDistribution
from .vegetation_dynamics import Vegetation
from .weathering import ExponentialWeatherer
<<<<<<< HEAD
from .depth_dependent_diffusion import DepthDependentDiffuser
from .cubic_nonlinear_hillslope_flux import CubicNonLinearDiffuser
from .depth_dependent_cubic_soil_creep import DepthDependentCubicDiffuser
from .erosion_deposition import ErosionDeposition
from .space import Space
from .landslides import LandslideProbability
from .conduit_networks import PresFlowNetwork
from .conduit_networks import MeltCreep

COMPONENTS = [ChiFinder, LinearDiffuser,
              Flexure, FlowRouter, DepressionFinderAndRouter,
              PerronNLDiffuse, OverlandFlowBates, OverlandFlow,
              KinwaveImplicitOverlandFlow,
              PotentialEvapotranspiration, PotentialityFlowRouter,
              Radiation, SinkFiller,
              StreamPowerEroder, StreamPowerSmoothThresholdEroder,
              FastscapeEroder, SedDepEroder,
              PrecipitationDistribution,
              SteepnessFinder, DetachmentLtdErosion, gFlex,
              SoilInfiltrationGreenAmpt, FireGenerator,
              SoilMoisture, Vegetation, VegCA, DrainageDensity,
              ExponentialWeatherer, DepthDependentDiffuser,
              CubicNonLinearDiffuser, DepthSlopeProductErosion,
              FlowDirectorD8, FlowDirectorSteepest, FlowDirectorMFD,
              FlowDirectorDINF, FlowAccumulator, Space, ErosionDeposition,
              LandslideProbability, DepthDependentCubicDiffuser,
              PresFlowNetwork, MeltCreep]
=======

COMPONENTS = [
    ChiFinder,
    LinearDiffuser,
    Flexure,
    FlowRouter,
    DepressionFinderAndRouter,
    PerronNLDiffuse,
    OverlandFlowBates,
    OverlandFlow,
    KinwaveImplicitOverlandFlow,
    KinwaveOverlandFlowModel,
    PotentialEvapotranspiration,
    PotentialityFlowRouter,
    Radiation,
    SinkFiller,
    SinkFillerBarnes,
    LakeMapperBarnes,
    StreamPowerEroder,
    StreamPowerSmoothThresholdEroder,
    FastscapeEroder,
    SedDepEroder,
    PrecipitationDistribution,
    SpatialPrecipitationDistribution,
    SteepnessFinder,
    DetachmentLtdErosion,
    gFlex,
    SoilInfiltrationGreenAmpt,
    FireGenerator,
    SoilMoisture,
    Vegetation,
    VegCA,
    DrainageDensity,
    ExponentialWeatherer,
    DepthDependentDiffuser,
    TaylorNonLinearDiffuser,
    DepthSlopeProductErosion,
    FlowDirectorD8,
    FlowDirectorSteepest,
    FlowDirectorMFD,
    FlowDirectorDINF,
    FlowAccumulator,
    LossyFlowAccumulator,
    Space,
    ErosionDeposition,
    LandslideProbability,
    DepthDependentTaylorDiffuser,
    NormalFault,
    Lithology,
    LithoLayers,
    TransportLengthHillslopeDiffuser,
]
>>>>>>> origin/master

__all__ = [cls.__name__ for cls in COMPONENTS]
