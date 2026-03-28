from fastapi import FastAPI
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel
from typing import List
from hydraulics import run_simulation

app = FastAPI()

# -----------------------------
# SERVE STATIC FILES
# -----------------------------
app.mount("/static", StaticFiles(directory="static"), name="static")

@app.get("/")
def home():
    return FileResponse("static/index.html")

# -----------------------------
# MODELS
# -----------------------------
class Fluid(BaseModel):
    mw: float
    fann_600: float
    fann_300: float
    fann_200: float
    fann_100: float
    fann_6: float
    fann_3: float

class Temperature(BaseModel):
    surface_temp: float
    bhct: float
    geothermal_grad: float

class Trajectory(BaseModel):
    md: float
    inc: float
    azi: float
    tvd: float

class WellSection(BaseModel):
    type: str
    top_md: float
    end_md: float
    casing_id: float
    hole_d: float

class BHAComponent(BaseModel):
    name: str
    od: float
    id: float
    length: float

class Cuttings(BaseModel):
    rop: float
    concentration: float
    size: float
    density: float
    rpm: float

class SimulationInput(BaseModel):
    flowrate: float
    depth: float
    mpd_mode: str
    gradient_mode: str
    fluid: Fluid
    temperature: Temperature
    trajectory: List[Trajectory]
    well_sections: List[WellSection]
    bha: List[BHAComponent]
    cuttings: Cuttings

# -----------------------------
# API
# -----------------------------
@app.post("/simulate")
def simulate(data: SimulationInput):
    return run_simulation(data)