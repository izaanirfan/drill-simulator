from fastapi import FastAPI
from fastapi.staticfiles import StaticFiles
from fastapi.responses import FileResponse

from models import SimulationInput
from hydraulics import run_simulation
from export_excel import generate_excel
from export_pdf import generate_pdf

# ---------------------------------------------------
# INIT APP
# ---------------------------------------------------
app = FastAPI()

# ---------------------------------------------------
# SERVE FRONTEND (STATIC FILES)
# This removes ALL CORS issues
# ---------------------------------------------------
app.mount("/", StaticFiles(directory="static", html=True), name="static")

# ---------------------------------------------------
# HEALTH CHECK (optional but useful)
# ---------------------------------------------------
@app.get("/health")
def health_check():
    return {"status": "ok"}

# ---------------------------------------------------
# SIMULATION ENDPOINT
# ---------------------------------------------------
@app.post("/simulate")
def simulate(data: SimulationInput):
    print("Simulation started")
    result = run_simulation(data)
    print("Simulation completed")
    return result

# ---------------------------------------------------
# EXPORT EXCEL
# ---------------------------------------------------
@app.post("/export/excel")
def export_excel_endpoint(data: SimulationInput):
    results = run_simulation(data)
    file_path = generate_excel(results)
    return FileResponse(file_path, filename="drilling_results.xlsx")

# ---------------------------------------------------
# EXPORT PDF
# ---------------------------------------------------
@app.post("/export/pdf")
def export_pdf_endpoint(data: SimulationInput):
    results = run_simulation(data)
    file_path = generate_pdf(results)
    return FileResponse(file_path, filename="drilling_report.pdf")