import pandas as pd
import xarray as xr
import rioxarray as rio
import geopandas as gpd


nuts21 = gpd.read_file("/home/roncmic/MultiHazard/data/NUTS_RG_20M_2021_3035.shp")
nuts3 = nuts21[nuts21["LEVL_CODE"]==3]

ws_path = "/eos/jeodpp/data/projects/CLIMEX/ERA5-WindGust/DailyMax"


def combine_nc_files(edo_path):
    nc_files = [os.path.join(edo_path, f) for f in sorted(os.listdir(edo_path)) if f.endswith('.nc')]
    combined_ds = xr.open_mfdataset(nc_files, concat_dim='time', combine='nested')
    combined_ds.rio.write_crs("EPSG:4326", inplace=True)
    combined_ds_projected = combined_ds.rio.reproject("EPSG:3035")
    return combined_ds_projected

# Combine the dataset
projected_ds = combine_nc_files(ws_path)


def detect_wind_speed_events(nuts_ws, threshold):
    events = []
    nuts_ws = nuts_ws.fillna(0)
    max_ws = nuts_ws.max(dim=['x', 'y'], skipna=True)
    ws_series = max_ws["i10fg"].to_series()
    event_in_progress = False
    start_date = None
    end_date = None
    max_intensity = 0

    for date, value in ws_series.items():
        if value > threshold and not event_in_progress:
            # Start a new wind speed event
            event_in_progress = True
            start_date = date
            max_intensity = value
        elif value > threshold and event_in_progress:
            # Continue the current wind speed event
            max_intensity = max(max_intensity, value)
        elif value <= threshold and event_in_progress:
            # End the wind speed event
            end_date = date
            duration = (end_date - start_date).days
            intensity_percentile = calculate_full_percentile(subws, max_intensity, percentile_99)
            events.append({
                'begin': start_date,
                'end': end_date,
                'intensity': max_intensity,
                'intensity_percentile': intensity_percentile,
                'duration': duration
            })
            event_in_progress = False
            start_date = None
            end_date = None
            max_intensity = 0

    # Handle the case where the event is still in progress at the end of the series
    if event_in_progress:
        end_date = ws_series.index[-1]
        duration = (end_date - start_date).days
        intensity_percentile = calculate_full_percentile(subws, max_intensity, percentile_99)
        events.append({
            'begin': start_date,
            'end': end_date,
            'intensity': max_intensity,
            'intensity_percentile': intensity_percentile,
            'duration': duration
        })

    return events

# Iterate over each NUTS3 region and compute wind speed statistics
results_list = []

for index, row in nuts3.iterrows():
    nuts_id = row['NUTS_ID']  # Adapt the 'nuts_id' field based on your shapefile
    geometry = row['geometry']
    print(nuts_id)
    
    # Clip or mask the wind speed data to the current NUTS3 region
    try:
        nuts_ws = projected_ds.rio.clip([geometry], projected_ds.rio.crs, all_touched=True)
        #mask_values_if_exceeds_threshold(nuts_ws, 1000)
        wind_speed_events = detect_wind_speed_events(nuts_ws, percentile_99)
    except Exception as e:
        print(f"No data found for NUTS ID: {nuts_id}")
        wind_speed_events = [{'begin': 0, 'end': 0, 'intensity': 0, 'intensity_percentile':0, 'duration': 0}]
    
    # Append results for this NUTS3 region to the results list
    for event in wind_speed_events:
        results_list.append({
            'nuts_id': nuts_id,
            'begin': event['begin'],
            'end': event['end'],
            'intensity': event['intensity'],
            'intensity_percentile': event['intensity_percentile'],
            'duration': event['duration']
        })

# Convert the results list to a Pandas DataFrame
results_df = pd.DataFrame(results_list)