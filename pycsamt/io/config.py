
import pandas as pd

class Config:
    """
    Container of un‐modifiable I/O settings and lookup tables.
    """

    @property
    def parsers(self) -> dict[str, callable]:
        """Mapping of file extensions → pandas reader functions."""
        return {
            ".csv":   pd.read_csv,
            ".xlsx":  pd.read_excel,
            ".json":  pd.read_json,
            ".html":  pd.read_html,
            ".sql":   pd.read_sql,
            ".xml":   pd.read_xml,
            ".fwf":   pd.read_fwf,
            ".pkl":   pd.read_pickle,
            ".sas":   pd.read_sas,
            ".spss":  pd.read_spss,
            # add more as needed
        }

    @staticmethod
    def writers(obj) -> dict[str, callable]:
        """
        Mapping of file extensions → DataFrame writer methods.
        Pass in the DataFrame (or similar) instance as `obj`.
        """
        return {
            ".csv":    obj.to_csv,
            ".hdf":    obj.to_hdf,
            ".sql":    obj.to_sql,
            ".dict":   obj.to_dict,
            ".xlsx":   obj.to_excel,
            ".json":   obj.to_json,
            ".html":   obj.to_html,
            ".feather":obj.to_feather,
            ".tex":    obj.to_latex,
            ".stata":  obj.to_stata,
            ".gbq":    obj.to_gbq,
            ".rec":    obj.to_records,
            ".str":    obj.to_string,
            ".clip":   obj.to_clipboard,
            ".md":     obj.to_markdown,
            ".parq":   obj.to_parquet,
            ".pkl":    obj.to_pickle,
            # add more as needed
        }

