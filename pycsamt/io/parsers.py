import os

from .config import Config


def read_any(path: str, **kwargs):
    """
    Read a table‐like file into a pandas object using the Config parser map.

    Parameters
    ----------
    path : str
        Path to file. Extension must be one of Config().parsers.keys().
    **kwargs
        Passed directly to the pandas reader (e.g. sep, sheet_name,
        parse_dates, etc.)

    Returns
    -------
    DataFrame or similar
    """
    ext = os.path.splitext(path)[1].lower()
    reader = Config().parsers.get(ext)
    if reader is None:
        raise ValueError(f"No parser configured for extension '{ext}'")
    return reader(path, **kwargs)


def write_any(obj, path: str, **kwargs):
    """
    Write a pandas‐like object to disk using Config writers.

    Parameters
    ----------
    obj : pandas.DataFrame or similar
    path : str
        Output path with one of the Config.writers() extensions.
    **kwargs
        Passed directly to the pandas writer (e.g. index, compression, etc.)
    """
    ext = os.path.splitext(path)[1].lower()
    writer_map = Config.writers(obj)
    writer = writer_map.get(ext)
    if writer is None:
        raise ValueError(f"No writer configured for extension '{ext}'")
    # many writer methods take (path, **kwargs)
    return writer(path, **kwargs)
