import numpy as np
import pandas as pd


def process_events(
    events_df,
    reference_df,
    event_idx,
    window,
    nut_id,
    haz="flood",
    endpeaks_col="enddate",
    stpeaks_col="stardate",
    start_date_col="Start date",
    end_date_col="End date",
    id_col="eid"
):
    # Convert relevant columns to datetime (safely)
    events_df[endpeaks_col] = pd.to_datetime(events_df[endpeaks_col], errors='coerce')
    events_df[stpeaks_col] = pd.to_datetime(events_df[stpeaks_col], errors='coerce')
    reference_df[start_date_col] = pd.to_datetime(reference_df[start_date_col], errors='coerce')
    reference_df[end_date_col] = pd.to_datetime(reference_df[end_date_col], errors='coerce')

    # Extract reference event
    ref_start = reference_df.iloc[event_idx][start_date_col]
    ref_end = reference_df.iloc[event_idx][end_date_col]
    ref_id = reference_df.iloc[event_idx,][id_col]

    # Calculate date differences (equivalent to R’s d1, d2)
    d1 = (events_df[endpeaks_col] - ref_start).dt.days + 1
    d2 = (events_df[stpeaks_col] - ref_end).dt.days + 1

    # Determine which events to keep, hazard-specific logic
    if haz == "flood":
        keep_mask = (d1 <= -4) & (d1 >= -window)
    elif haz == "flmatch":
        keep_mask = (d2 <= window) & (d1 >= -window)
    elif haz in ["drought", "cold", "wind"]:
        keep_mask = (d2 <= 0) & (d1 >= -window)
    elif haz == "heat":
        keep_mask = (d1 <= 0) & (d1 >= -window)
    else:
        print("Invalid hazard type.")
        return pd.DataFrame()  # Empty dataframe

    # Keep the matching rows
    events_to_keep = events_df.loc[keep_mask].copy()

    if not events_to_keep.empty:
        events_to_keep["dtime1"] = d1[keep_mask].values
        events_to_keep["dtime2"] = d1[keep_mask].values
        events_to_keep["NUTID"] = nut_id
        events_to_keep["eventStart"] = ref_start
        events_to_keep["eventEnd"] = ref_end
        events_to_keep["eventID"] = ref_id
        return events_to_keep

    # If no matches found, return empty DataFrame
    return pd.DataFrame()


def elbow_finder(x_values, y_values, n_points=16):
    """
    Finds the 'elbow point' of a curve by computing the farthest distance 
    between the curve and a straight line connecting the min and max values.

    Parameters
    ----------
    x_values : array-like
        X coordinates of the curve.
    y_values : array-like
        Y coordinates of the curve.
    n_points : int
        Number of points to interpolate between min and max (default: 16).

    Returns
    -------
    tuple
        (x_max_dist, y_max_dist, xcmax, ycmax)
        where:
        - (x_max_dist, y_max_dist) = coordinates of the elbow on the curve
        - (xcmax, ycmax) = coordinates on the baseline line directly below the elbow
    """

    x_values = np.array(x_values, dtype=float)
    y_values = np.array(y_values, dtype=float)

    # Identify min-x and max-y points
    max_x_x = np.min(x_values)
    max_x_y = y_values[np.argmin(x_values)]
    max_y_y = np.max(y_values)
    max_y_x = x_values[np.argmax(y_values)]

    # Line between the two extreme points (max_x and max_y)
    line_x = np.linspace(max_y_x, max_x_x, n_points)
    line_y = np.linspace(max_y_y, max_x_y, n_points)

    # Normalize both line and curve
    xcurve = line_y / np.max(line_y)
    ycurve = line_x / np.max(line_x)
    x_values2 = x_values / np.max(x_values)
    y_values2 = y_values / np.max(y_values)

    # Compute minimal distances from curve to line
    distances = []
    x_min_dist = []
    y_min_dist = []

    for i in range(len(x_values2)):
        d = np.sqrt((x_values2[i] - xcurve)**2 + (y_values2[i] - ycurve)**2)
        min_idx = np.argmin(d)
        distances.append(np.min(d))
        x_min_dist.append(xcurve[min_idx])
        y_min_dist.append(ycurve[min_idx])

    distances = np.array(distances)
    x_min_dist = np.array(x_min_dist)
    y_min_dist = np.array(y_min_dist)

    # Find point with maximum distance
    max_idx = np.argmax(distances)
    x_max_dist = x_values[max_idx]
    y_max_dist = y_values[max_idx]
    xcmax = x_min_dist[max_idx]
    ycmax = y_min_dist[max_idx]

    return (x_max_dist, y_max_dist, xcmax, ycmax)
