from pydantic import BaseModel
from typing import List, Optional

class BHAComponent(BaseModel):
    name: str
    od: float
    id: float
    length: float

class Fluid(BaseModel):
    mw: float
    fann_600: float
    fann_300: float
    fann_200: float
    fann_100: float
    fann_6: float
    fann_3: float

class TemperatureProfile(BaseModel):
    surface_temp: float
    bhct: float
    geothermal_grad: float
    beta: float = 0.0003

class FluidSection(BaseModel):
    top: float
    bottom: float
    mw: float

class SurveyPoint(BaseModel):
    md: float
    inc: float
    azi: float

class WellSection(BaseModel):
    type: str
    top_md: float
    end_md: float
    casing_id: float
    hole_d: float

class Cuttings(BaseModel):
    rop: float
    size: float
    density: float
    rpm: float
    concentration: Optional[float] = None

class SimulationInput(BaseModel):
    flowrate: float
    sbp: float = 0.0
    target_bhp: Optional[float] = None
    mpd_mode: str

    depth: float
    hole_diameter: float

    bha: List[BHAComponent]
    fluid: Fluid
    temperature: TemperatureProfile

    gradient_mode: str
    fluid_sections: List[FluidSection] = []

    trajectory: List[SurveyPoint]
    well_sections: List[WellSection]
    cuttings: Cuttings