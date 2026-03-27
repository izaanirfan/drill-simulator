from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse

from models import SimulationInput
from hydraulics import run_simulation
from export_excel import generate_excel
from export_pdf import generate_pdf

app = FastAPI()

# ✅ CORS CONFIG (FIXED)
app.add_middleware(
    CORSMiddleware,
    allow_origins=[
        "https://drill-simulator.vercel.app"
    ],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# ✅ SIMULATION ENDPOINT
@app.post("/simulate")
def simulate(data: SimulationInput):
    return run_simulation(data)

# ✅ EXPORT EXCEL
@app.post("/export/excel")
def export_excel(data: SimulationInput):
    results = run_simulation(data)
    file_path = generate_excel(results)
    return FileResponse(file_path, filename="drilling_results.xlsx")

# ✅ EXPORT PDF
@app.post("/export/pdf")
def export_pdf(data: SimulationInput):
    results = run_simulation(data)
    file_path = generate_pdf(results)
    return FileResponse(file_path, filename="drilling_report.pdf")
