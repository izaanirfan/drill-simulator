from fastapi import FastAPI
from fastapi.staticfiles import StaticFiles
from fastapi.responses import FileResponse
from pydantic import BaseModel
from typing import List
import uvicorn

from hydraulics import run_simulation  # make sure this matches your file

app = FastAPI()

# -------------------------------
# SERVE FRONTEND (STATIC FILES)
# -------------------------------
app.mount("/static", StaticFiles(directory="static"), name="static")

@app.get("/")
def serve_ui():
    return FileResponse("static/index.html")


# -------------------------------
# DATA MODELS
# -------------------------------

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

class Temperature(BaseModel):
    surface_temp: float
    bhct: float
    geothermal_grad: float

class TrajectoryPoint(BaseModel):
    md: float
    inc: float
    azi: float
    tvd: float

class WellSection(BaseModel):
    end_md: float
    hole_d: float
    pipe_od: float

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

    trajectory: List[TrajectoryPoint]
    well_sections: List[WellSection]
    bha: List[BHAComponent]
    fluid: Fluid
    temperature: Temperature
    cuttings: Cuttings


# -------------------------------
# API ENDPOINT
# -------------------------------

@app.post("/simulate")
def simulate(data: SimulationInput):

    try:
        result = run_simulation(data)
        return result

    except Exception as e:
        return {
            "error": str(e)
        }


# -------------------------------
# HEALTH CHECK
# -------------------------------
@app.get("/health")
def health():
    return {"status": "ok"}


# -------------------------------
# LOCAL RUN (OPTIONAL)
# -------------------------------
if __name__ == "__main__":
    uvicorn.run(app, host="0.0.0.0", port=10000)
