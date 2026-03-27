from openpyxl import Workbook

def generate_excel(results):
    wb = Workbook()
    ws = wb.active

    ws.append(["Depth", "ECD"])

    for i in range(len(results["depth_profile"])):
        ws.append([
            results["depth_profile"][i],
            results["ecd_profile"][i]
        ])

    file_path = "output.xlsx"
    wb.save(file_path)

    return file_path