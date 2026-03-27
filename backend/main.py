from fastapi import FastAPI
from fastapi.staticfiles import StaticFiles
from fastapi.responses import FileResponse

from models import SimulationInput
from hydraulics import run_simulation
from export_excel import generate_excel
from export_pdf import generate_pdf

app = FastAPI()

# ---------------------------------------------------
# API ROUTES FIRST
# ---------------------------------------------------
@app.post("/simulate")
def simulate(data: SimulationInput):
    print("Simulation started")
    result = run_simulation(data)
    print("Simulation completed")
    return result


@app.post("/export/excel")
def export_excel(data: SimulationInput):
    results = run_simulation(data)
    file_path = generate_excel(results)
    return FileResponse(file_path, filename="drilling_results.xlsx")


@app.post("/export/pdf")
def export_pdf(data: SimulationInput):
    results = run_simulation(data)
    file_path = generate_pdf(results)
    return FileResponse(file_path, filename="drilling_report.pdf")

# ---------------------------------------------------
# SERVE FRONTEND (IMPORTANT: NOT "/")
# ---------------------------------------------------
app.mount("/static", StaticFiles(directory="static"), name="static")

# Root route to serve UI
@app.get("/")
def serve_ui():
    return FileResponse("static/index.html")