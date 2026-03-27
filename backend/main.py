from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse

from models import SimulationInput
from hydraulics import run_simulation
from export_excel import generate_excel
from export_pdf import generate_pdf

app = FastAPI()

from fastapi.middleware.cors import CORSMiddleware

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # allow all
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

@app.post("/simulate")
def simulate(data: SimulationInput):
    return run_simulation(data)

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