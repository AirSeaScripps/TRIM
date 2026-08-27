import cdsapi

def get_era5_data(params):
    dataset = "reanalysis-era5-single-levels"

    request = {
        "product_type": ["reanalysis"],
        "variable": [params["variable"]],
        "year": params["year"],
        "month": params["month"],
        "day": params["day"],
        "time": params["time"],
        "data_format": "netcdf",
        "download_format": "unarchived",
        "area": params["area"]
    }

    client = cdsapi.Client()
    print("Submitting request to CDS API...", flush=True)

    client.retrieve(dataset, request).download()
