import pandas as _pd


def get_fit_results(survey, model, params):
    """Get lightcurve fits for a given survey and model

    Args:
        survey (str): The name of the survey
        model  (str): The name of the fitted model
        params (int): The number of Salt2 params that were fit (4 or 5)

    Returns:
        A DataFrame of fits in all bands
        A DataFrame of fits in blue bands
        A DataFrame of fits in red bands
    """

    index_col = 0
    path_pattern = f'./pipeline_outputs/{survey}_{model}_{params}param_{{}}.csv'

    all_data = _pd.read_csv(path_pattern.format('all'), index_col=index_col)

    blue_data = _pd.read_csv(path_pattern.format('blue'),
                             index_col=index_col)

    red_data = _pd.read_csv(path_pattern.format('red'),
                            index_col=index_col)

    return all_data, blue_data, red_data
