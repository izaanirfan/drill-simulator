from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
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
# CORS CONFIG (VERY IMPORTANT)
# ---------------------------------------------------
origins = [
    "https://drill-simulator.vercel.app",  # your frontend
    "http://localhost:3000",               # optional (local testing)
    "*"                                   # fallback (safe for now)
]

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# ---------------------------------------------------
# PREFLIGHT HANDLER (FIXES FAILED FETCH)
# ---------------------------------------------------
@app.options("/{full_path:path}")
async def preflight_handler(full_path: str):
    return {"message": "OK"}

# ---------------------------------------------------
# ROOT HEALTH CHECK
# ---------------------------------------------------
@app.get("/")
def read_root():
    return {"message": "Drilling Simulator API is running"}

# ---------------------------------------------------
# SIMULATION ENDPOINT
# ---------------------------------------------------
@app.post("/simulate")
def simulate(data: SimulationInput):
    print("Simulation started")  # debug log
    result = run_simulation(data)
    print("Simulation completed")
    return result

# ---------------------------------------------------
# EXPORT EXCEL
# ---------------------------------------------------
@app.post("/export/excel")
def export_excel(data: SimulationInput):
    results = run_simulation(data)
    file_path = generate_excel(results)
    return FileResponse(file_path, filename="drilling_results.xlsx")

# ---------------------------------------------------
# EXPORT PDF
# ---------------------------------------------------
@app.post("/export/pdf")
def export_pdf(data: SimulationInput):
    results = run_simulation(data)
    file_path = generate_pdf(results)
    return FileResponse(file_path, filename="drilling_report.pdf")