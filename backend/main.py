from fastapi import FastAPI
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles

# Import models directly from your models.py file
from models import SimulationInput
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
# API
# -----------------------------
@app.post("/simulate")
def simulate(data: SimulationInput):
    return run_simulation(data)