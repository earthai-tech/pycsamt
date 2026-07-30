# License: BSD-3-Clause

"""
`util` module provides utility functions and classes to support various
data manipulation and analysis tasks
"""

import os
import re
import shutil
import warnings

import numpy as np
import pandas as pd


class TerminalSize:
    """
    A utility class to get the size of the terminal window. This class provides
    methods to obtain the terminal size across different platforms, ensuring
    compatibility and providing a fallback size when necessary.

    Attributes
    ----------
    DEFAULT_SIZE : tuple
        A default terminal size to fallback to when the actual terminal size
        cannot be determined. The default is (80, 20).

    Methods
    -------
    get_terminal_size()
        Get the size of the terminal window as (columns, rows). Attempts multiple
        methods to ensure compatibility across different platforms.
    _get_terminal_size_windows()
        Get the terminal size for Windows platforms using ctypes.
    _get_terminal_size_unix()
        Get the terminal size for Unix-like platforms using ioctl syscall.

    Notes
    -----
    The `get_terminal_size` method first tries to use `shutil.get_terminal_size`
    with a fallback to `DEFAULT_SIZE`. If that fails, it checks the operating
    system and attempts to use platform-specific methods. On Windows, it uses
    the `ctypes` library to query the terminal size. On Unix-like systems, it
    uses the `ioctl` syscall.

    Examples
    --------
    >>> from fusionlab.api.util import TerminalSize
    >>> TerminalSize.get_terminal_size()
    (80, 25)

    See Also
    --------
    shutil.get_terminal_size : Get the size of the terminal window.

    References
    ----------
    .. [1] Python Software Foundation. Python Language Reference, version 3.9.
       Available at http://www.python.org
    .. [2] Python `shutil` module documentation. Available at
       https://docs.python.org/3/library/shutil.html#get_terminal_size
    """

    DEFAULT_SIZE = (80, 20)

    @staticmethod
    def get_terminal_size():
        """
        Get the size of the terminal window as (columns, rows).
        Attempts multiple methods to ensure compatibility across different
        platforms.

        Returns
        -------
        tuple
            The size of the terminal as (columns, rows).
        """
        try:
            size = shutil.get_terminal_size(fallback=TerminalSize.DEFAULT_SIZE)
            return size.columns, size.lines
        except Exception:
            pass  # Ignore and try other methods

        if os.name == "nt":
            return TerminalSize._get_terminal_size_windows()
        elif os.name == "posix":
            return TerminalSize._get_terminal_size_unix()

        return TerminalSize.DEFAULT_SIZE

    @staticmethod
    def _get_terminal_size_windows():
        """
        Get the terminal size for Windows platforms using ctypes.

        Returns
        -------
        tuple
            The size of the terminal as (columns, rows).
        """
        import struct
        from ctypes import create_string_buffer, windll

        try:
            h = windll.kernel32.GetStdHandle(-12)
            csbi = create_string_buffer(22)
            res = windll.kernel32.GetConsoleScreenBufferInfo(h, csbi)
            if res:
                (_, _, _, _, _, left, top, right, bottom, _, _) = (
                    struct.unpack("hhhhHhhhhhh", csbi.raw)
                )
                sizex = right - left + 1
                sizey = bottom - top + 1
                return sizex, sizey
        except Exception:
            pass  # Ignore and fallback to default

        return TerminalSize.DEFAULT_SIZE

    @staticmethod
    def _get_terminal_size_unix():
        """
        Get the terminal size for Unix-like platforms using ioctl syscall.

        Returns
        -------
        tuple
            The size of the terminal as (columns, rows).
        """
        import fcntl
        import struct
        import termios

        def ioctl_GWINSZ(fd):
            try:
                return struct.unpack(
                    "hh", fcntl.ioctl(fd, termios.TIOCGWINSZ, "1234")
                )
            except Exception:
                return None

        size = ioctl_GWINSZ(0) or ioctl_GWINSZ(1) or ioctl_GWINSZ(2)
        if not size:
            try:
                with open(os.ctermid(), "rb") as fd:
                    size = ioctl_GWINSZ(fd.fileno())
            except Exception:
                pass  # Ignore and fallback to default

        if size:
            return size[1], size[0]

        return TerminalSize.DEFAULT_SIZE


def format_value(value, precision=4):
    """
    Format a numeric value to a string, rounding floats to four decimal
    places and converting integers directly to strings.

    Parameters
    ----------
    value : int, float, np.integer, np.floating
        The numeric value to be formatted.

    Returns
    -------
    str
        A formatted string representing the value.

    Examples
    --------
    >>> from fusionlab.api.util import format_value
    >>> format_value(123)
    '123'
    >>> format_value(123.45678)
    '123.4568'
    """
    value_str = str(value)
    precision = validate_precision(precision)
    if isinstance(value, (int, float, np.integer, np.floating)):
        value = apply_precision(value, precision)
        value_str = f"{value}"  # if isinstance (
        #     value, int) else  f"{float(value):.{precision}f}"
    return value_str


def apply_precision(value, precision=4):
    """
    Applies precision to a numerical value only if the decimal part is
    larger than the specified precision. Otherwise, returns the value
    as is.

    Parameters
    ----------
    value : int, float, np.integer, np.floating
        The numerical value to be formatted.
    precision : int, optional
        The number of decimal places to format floating-point numbers.
        Default is 4.

    Returns
    -------
    int, float
        The formatted value if the decimal part is larger than the
        specified precision; otherwise, the value as is.

    Examples
    --------
    >>> from fusionlab.api.util import apply_precision
    >>> apply_precision(123.456789, 2)
    123.46
    >>> apply_precision(123.4, 2)
    123.4
    >>> apply_precision(np.int32(456), 2)
    456
    >>> apply_precision(np.float64(456.78), 2)
    456.78
    >>> apply_precision(123, 2)
    123
    >>> apply_precision(123.00, 2)
    123.0
    """
    precision = validate_precision(precision)
    if isinstance(value, (int, np.integer)):
        return int(value)
    elif isinstance(value, (float, np.floating)):
        rounded_value = round(value, precision)
        if value == rounded_value:
            return value
        return rounded_value
    return value


def validate_precision(precision, /):
    """
    Validates and converts the precision parameter to ensure it is a
    non-negative integer.

    Parameters
    ----------
    precision : int, float, np.integer, np.floating
        The precision value to be validated.

    Returns
    -------
    int
        The validated and converted precision value.

    Raises
    ------
    ValueError
        If the precision is not a non-negative number or cannot be
        converted to a non-negative integer.

    Examples
    --------
    >>> validate_precision(3)
    3
    >>> validate_precision(3.0)
    3
    >>> validate_precision(np.int32(2))
    2
    >>> validate_precision(np.float64(2.0))
    2
    >>> validate_precision(-1)
    Traceback (most recent call last):
        ...
    ValueError: Precision must be a non-negative integer.
    >>> validate_precision("three")
    Traceback (most recent call last):
        ...
    ValueError: Precision must be a non-negative integer.
    """
    precision = precision or 4
    try:
        precision = int(precision)
        if precision < 0:
            raise ValueError
    except (ValueError, TypeError):
        raise ValueError("Precision must be a non-negative integer.")
    return precision


def parse_component_kind(pc_list, kind):
    """
    Extracts specific principal component's feature names and their importance
    values from a list based on a given component identifier.

    Parameters
    ----------
    pc_list : list of tuples
        A list where each tuple contains ('pc{i}', feature_names,
                                          sorted_component_values),
        corresponding to each principal component. 'pc{i}' is a string label
        like 'pc1', 'feature_names' is an array of feature names sorted by
        their importance, and 'sorted_component_values' are the corresponding
        sorted values of component loadings.
    kind : str
        A string that identifies the principal component number to extract,
        e.g., 'pc1'. The string should contain a numeric part that corresponds
        to the component index in `pc_list`.

    Returns
    -------
    tuple
        A tuple containing two elements:
        - An array of feature names for the specified principal component.
        - An array of sorted component values for the specified principal
          component.

    Raises
    ------
    ValueError
        If the `kind` parameter does not contain a valid component number or if the
        specified component number is out of the range of available components
        in `pc_list`.

    Examples
    --------
    >>> from fusionlab.api.extension import parse_component_kind
    >>> pc_list = [
    ...     ("pc1", ["feature1", "feature2", "feature3"], [0.8, 0.5, 0.3]),
    ...     ("pc2", ["feature1", "feature2", "feature3"], [0.6, 0.4, 0.2]),
    ... ]
    >>> feature_names, importances = parse_component_kind(pc_list, "pc1")
    >>> print(feature_names)
    ['feature1', 'feature2', 'feature3']
    >>> print(importances)
    [0.8, 0.5, 0.3]

    Notes
    -----
    The function requires that the `kind` parameter include a numeric value
    that accurately represents a valid index in `pc_list`. The index is derived
    from the numeric part of the `kind` string and is expected to be 1-based.
    If no valid index is found or if the index is out of range, the function
    raises a `ValueError`.
    """
    match = re.search(r"\d+", str(kind))
    if match:
        # Convert 1-based index from `kind` to 0-based index for list access
        index = int(match.group()) - 1
        if index < len(pc_list) and index >= 0:
            return pc_list[index][1], pc_list[index][2]
        else:
            raise ValueError(
                f"Component index {index + 1} is out of the"
                " range of available components."
            )
    else:
        raise ValueError(
            "The 'kind' parameter must include an integer"
            " indicating the desired principal component."
        )


def find_maximum_table_width(summary_contents, header_marker="="):
    """
    Calculates the maximum width of tables in a summary string based on header lines.

    This function parses a multi-table summary string, identifying lines that represent
    the top or bottom borders of tables (header lines). It determines the maximum width
    of these tables by measuring the length of these header lines. The function assumes
    that the header lines consist of repeated instances of a specific marker character.

    Parameters
    ----------
    summary_contents : str
        A string containing the summarized representation of one or more tables.
        This string should include header lines made up of repeated header markers
        that denote the start and end of each table's border.
    header_marker : str, optional
        The character used to construct the header lines in the summary_contents.
        Defaults to '=', the common character for denoting table borders in ASCII
        table representations.

    Returns
    -------
    int
        The maximum width of the tables found in summary_contents, measured as the
        length of the longest header line. If no header lines are found, returns 0.

    Examples
    --------
    >>> from fusionlab.api.util import find_maximum_table_width
    >>> summary = '''Model Performance
    ... ===============
    ... Estimator : SVC
    ... Accuracy  : 0.9500
    ... Precision : 0.8900
    ... Recall    : 0.9300
    ... ===============
    ... Model Performance
    ... =================
    ... Estimator : RandomForest
    ... Accuracy  : 0.9500
    ... Precision : 0.8900
    ... Recall    : 0.9300
    ... ================='''
    >>> find_maximum_table_width(summary)
    18

    This example shows how the function can be used to find the maximum table width
    in a string containing summaries of model performances, where '=' is used as
    the header marker.
    """
    # Split the input string into lines
    lines = summary_contents.split("\n")
    # Filter out lines that consist only of the header
    # marker, and measure their lengths
    header_line_lengths = [
        len(line) for line in lines if line.strip(header_marker) == ""
    ]
    # Return the maximum of these lengths, or 0 if the list is empty
    return max(header_line_lengths, default=0)


def format_text(
    text,
    key=None,
    key_length=15,
    max_char_text=50,
    add_frame_lines=False,
    border_line="=",
    buffer_space=4,
):
    """
    Formats a block of text to fit within a specified maximum character width,
    optionally prefixing it with a key. If the text exceeds the maximum width,
    it wraps to a new line, aligning with the key or the specified indentation.

    Parameters
    ----------
    text : str
        The text to be formatted.
    key : str, optional
        An optional key to prefix the text. Defaults to None.
    key_length : int, optional
        The length reserved for the key, including following spaces.
        If `key` is provided but `key_length` is None, the length of the
        `key` plus one space is used. Defaults to 15.
    max_char_text : int, optional
        The maximum number of characters for the text width, including the key
        if present. Defaults to 50.
    add_frame_lines: bool, False
       If True, frame the text with '=' line (top and bottom)
    border_line: str, optional
      The border line to frame the text.  Default is '='
    buffer_space: int, default=4
       The extra space to force breaken the text. Using a large value will
       provide nice reading and force terminating the line sentence with no
       preprosition like `a`, 'of' or `the` etc.

    Returns
    -------
    str
        The formatted text with line breaks added to ensure that no line exceeds
        `max_char_text` characters. If a `key` is provided, it is included only
        on the first line, with subsequent lines aligned accordingly.

    Examples
    --------
    >>> from fusionlab.api.util import format_text
    >>> text_example = ("This is an example text that is supposed to wrap"
                      "around after a certain number of characters.")
    >>> print(format_text(text_example, key="Note"))
    Note           : This is an example text that is supposed to
                      wrap around after a certain number of
                      characters.

    Notes
    -----
    - The function dynamically adjusts the text to fit within `max_char_text`,
      taking into account the length of `key` if provided.
    - Text that exceeds the `max_char_text` limit is wrapped to new lines, with
      proper alignment to match the initial line's formatting.
    """
    if key is not None:
        # If key_length is None, use the length of the key + 1
        # for the space after the key
        if key_length is None:
            key_length = len(key)  # + 1
        key_str = f"{key.ljust(key_length)} : "
    elif key_length is not None:
        # If key is None but key_length is specified, use spaces
        key_str = " " * key_length + " : "
    else:
        # If both key and key_length are None, there's no key part
        key_str = ""

    # Adjust max_char_text based on the length of the key part: +3 for " : ".
    effective_max_char_text = (
        max_char_text - len(key_str) + buffer_space
        if key_str
        else max_char_text
    )
    formatted_text = ""
    text = str(text)
    while text:
        # If the remaining text is shorter than the effective
        # max length, or if there's no key part, add it as is
        if len(text) <= effective_max_char_text - buffer_space or not key_str:
            formatted_text += key_str + text
            break
        else:
            # Find the space to break the line, ensuring it doesn't
            # exceed effective_max_char_text
            break_point = text.rfind(
                " ", 0, effective_max_char_text - buffer_space
            )

            if break_point == -1:  # No spaces found, force break
                break_point = effective_max_char_text - buffer_space
            # Add the line to formatted_text
            formatted_text += key_str + text[:break_point].rstrip() + "\n"
            # Remove the added part from text
            text = text[break_point:].lstrip()

            # After the first line, the key part is just spaces
            key_str = " " * len(key_str)

    if add_frame_lines:
        frame_lines = (
            border_line * max_char_text
        )  # (effective_max_char_text + 1 )
        formatted_text = (
            frame_lines + "\n" + formatted_text + "\n" + frame_lines
        )

    return formatted_text


def get_frame_chars(frame_char):
    """
    Retrieve framing characters based on the input frame indicator.

    Parameters
    ----------
    frame_char : str
        A single character that indicates the desired framing style.

    Returns
    -------
    tuple
        A tuple containing the close character and the open-close pair
        for framing index values.

    Examples
    --------
    >>> from fusionlab.api.util import get_frame_chars
    >>> get_frame_chars("[")
    (']', '[', ']')
    >>> get_frame_chars("{")
    ('}', '{', '}')
    """
    pairs = {
        "[": ("]", "[", "]"),
        "{": ("}", "{", "}"),
        "(": (")", "(", ")"),
        "<": (">", "<", ">"),
    }
    return pairs.get(frame_char, (".", ".", "."))


def format_cell(x, max_text_length, max_width=None):
    """
    Truncates a string to the maximum specified length and appends '...'
    if needed, and right-aligns it.

    Parameters:
    x (str): The string to format.
    max_width (int): The width to which the string should be aligned.
    max_text_length (int): The maximum allowed length of the string before truncation.

    Returns:
    str: The formatted and aligned string.
    """
    x = str(x)
    if len(x) > max_text_length:
        x = x[: max_text_length - 3] + "..."
    return x.rjust(max_width) if max_width else x


def get_table_width_from(
    formatted_str,
    /,
    border_char="=",
    check_first=True,
    deep_check=True,
    width_strategy="max",
    error="warn",
):
    """
    Calculate the maximum table width from the given formatted string based on
    the border style.

    Parameters
    ----------
    formatted_str : str
        The string containing the formatted table.
    border_char : str, optional
        The character used for the table border. Default is '='.
    check_first : bool, optional
        If True, check the width of the first border line found and return it.
        If False, check for additional border lines if the first one is found.
        Default is True.
    deep_check : bool, optional
        If True, evaluate all lines containing the border character to determine
        the table width. If False, follow the `check_first` strategy. Default is
        True.
    width_strategy : {'max', 'min', 'average'}, optional
        Strategy to determine the table width if `deep_check` is True. 'max'
        returns the longest border line. 'min' returns the shortest border line.
        'average' returns the average width of all border lines. Default is 'max'.
    error : {'warn', 'raise', 'ignore'}, optional
        How to handle cases where no border line is found. 'warn' logs a warning
        message. 'raise' raises an exception. 'ignore' does nothing and returns
        None. Default is 'warn'.

    Returns
    -------
    int or None
        The calculated table width, or None if no border line is found and error
        handling is set to 'ignore'.

    Notes
    -----
    This function identifies lines in the `formatted_str` containing the
    `border_char` and calculates the width of the table based on the specified
    strategies and checks.

    The calculation follows these steps:

    1. Identify lines in `formatted_str` containing `border_char`.
    2. If no such lines are found, handle the error based on `error` parameter.
    3. If `deep_check` is False:
        - If `check_first` is True, return the length of the first border line.
        - Otherwise, return the length of the last border line.
    4. If `deep_check` is True:
        - Compute the lengths of all border lines.
        - If `width_strategy` is 'average', return the average length of border
          lines:

          .. math::
              \text{average\\_width} = \frac{\\sum_{i=1}^{n} \text{length}_i}{n}

        - If `width_strategy` is 'min', return the minimum length of border lines.
        - If `width_strategy` is 'max', return the maximum length of border lines.

    Examples
    --------
    >>> from fusionlab.api.util import get_table_width_from
    >>> formatted_str = "=======\n| Col |\n=======\n"
    >>> get_table_width_from(formatted_str)
    7

    >>> formatted_str = "=====\n| C |\n=====\n"
    >>> get_table_width_from(formatted_str, border_char="=")
    5

    See Also
    --------
    pandas.DataFrame : DataFrame structure used in pandas for data manipulation
        and analysis.

    References
    ----------
    .. [1] McKinney, Wes. "Data Structures for Statistical Computing in Python."
           Proceedings of the 9th Python in Science Conference. 2010.
    """
    border_lines = [
        line for line in formatted_str.splitlines() if border_char in line
    ]

    if not border_lines:
        if error == "warn":
            warnings.warn(
                " No border line found in the formatted string.", stacklevel=2
            )
        elif error == "raise":
            raise ValueError(
                "Error: No border line found in the formatted string."
            )
        return None

    if not deep_check:
        return len(border_lines[0]) if check_first else len(border_lines[-1])

    line_lengths = [len(line) for line in border_lines]

    if width_strategy == "average":
        return sum(line_lengths) // len(line_lengths)
    elif width_strategy == "min":
        return min(line_lengths)
    else:  # width_strategy == 'max'
        return max(line_lengths)


def generate_legend(
    custom_markers=None,
    no_corr_placeholder="...",
    hide_diag=True,
    max_width=50,
    add_frame_lines=True,
    border_line=".",
):
    """
    Generates a legend for a table (dataframe) matrix visualization, formatted
    according to specified parameters.

    This function supports both numeric and symbolic representations of table
    values. Symbolic representations, which are used primarily for visual clarity,
    include the following symbols:

    - ``'++'``: Represents a strong positive relationship.
    - ``'--'``: Represents a strong negative relationship.
    - ``'-+'``: Represents a moderate relationship.
    - ``'o'``: Used exclusively for diagonal elements, typically representing
      a perfect relationship in correlation matrices (value of 1.0).

    Parameters
    ----------
    custom_markers : dict, optional
        A dictionary mapping table symbols to their descriptions. If provided,
        it overrides the default markers. Default is None.
    no_corr_placeholder : str, optional
        Placeholder text for table values that do not meet the minimum threshold.
        Default is '...'.
    hide_diag : bool, optional
        If True, omits the diagonal entries from the legend. These are typically
        frame of a variable with itself (1.0). Default is True.
    max_width : int, optional
        The maximum width of the formatted legend text, influences the centering
        of the title. Default is 50.
    add_frame_lines : bool, optional
        If True, adds a frame around the legend using the specified `border_line`.
        Default is True.
    border_line : str, optional
        The character used to create the border of the frame if `add_frame_lines`
        is True. Default is '.'.

    Returns
    -------
    str
        The formatted legend text, potentially framed, centered according to the
        specified width, and including custom or default descriptions of correlation
        values.

    Examples
    --------
    >>> from fusionlab.api.util import generate_legend
    >>> custom_markers = {"++": "High Positive", "--": "High Negative"}
    >>> print(generate_legend(custom_markers=custom_markers, max_width=60))
    ............................................................
    Legend : ...: Non-correlated, ++: High Positive, --: High
             Negative, -+: Moderate
    ............................................................

    >>> print(generate_legend(hide_diag=False, max_width=70))
    ......................................................................
    Legend : ...: Non-correlated, ++: Strong positive, --: Strong negative,
             -+: Moderate, o: Diagonal
    ......................................................................
    >>> custom_markers = {"++": "Highly positive", "--": "Highly negative"}
    >>> legend = generate_legend(
    ...     custom_markers=custom_markers,
    ...     no_corr_placeholder="N/A",
    ...     hide_diag=False,
    ...     border_line="=",
    ... )

    >>> print(legend)

    ==================================================
    Legend : N/A: Non-correlated, ++: Highly positive,
             --: Highly negative, -+: Moderate, o:
             Diagonal
    ==================================================
    """
    # Default markers and their descriptions
    default_markers = {
        no_corr_placeholder: "Non-correlated",
        "++": "Strong positive",
        "--": "Strong negative",
        "-+": "Moderate",
        "o": "Diagonal",  # only used if hide_diag is False
    }
    if custom_markers is not None and not isinstance(custom_markers, dict):
        raise TypeError(
            "The 'custom_markers' parameter must be a dictionary."
            f" Received type: {type(custom_markers).__name__}. Please provide a dictionary"
            " where keys are the legend symbols and values"
            " are their descriptions."
        )

    # Update default markers with any custom markers provided
    markers = {**default_markers, **(custom_markers or {})}
    # If no correlation placeholder, then remove it from the markers.
    if not no_corr_placeholder:
        markers.pop(no_corr_placeholder)
    # Create legend entries
    legend_entries = [
        f"{key}: {value}"
        for key, value in markers.items()
        if not (key == "o" and hide_diag)
    ]

    # Join entries with commas and format the legend text
    legend_text = ", ".join(legend_entries)
    legend = "\n\n" + format_text(
        legend_text,
        key="Legend",
        key_length=len("Legend"),
        max_char_text=max_width,  # + len('Legend'),
        add_frame_lines=add_frame_lines,
        border_line=border_line,
    )
    return legend


def to_snake_case(name, mode="standard"):
    """
    Converts a string to snake_case.

    Parameters
    ----------
    name : str
        The string to convert to snake_case.
    mode : str, optional
        If 'soft', extra whitespace and case inconsistencies are
        handled to produce clean snake_case.

    Returns
    -------
    str
        The snake_case version of the input string.
    """
    name = str(name)

    if mode == "soft":
        # Convert to lowercase and replace multiple spaces
        # or non-word characters with a single underscore
        name = re.sub(
            r"\W+", " ", name
        )  # Replace non-word characters with spaces
        name = re.sub(r"\s+", " ", name).strip()  # Normalize whitespace
        name = name.lower().replace(" ", "_")  # Convert spaces to underscores

    else:
        # Standard snake_case conversion without additional processing
        name = re.sub(
            r"(?<!^)(?=[A-Z])", "_", name
        ).lower()  # Convert CamelCase to snake_case
        name = re.sub(
            r"\W+", "_", name
        )  # Replace non-word characters with '_'
        name = re.sub(
            r"_+", "_", name
        )  # Replace multiple underscores with a single '_'

    return name.strip("_")  # Remove any leading or trailing underscores


def generate_column_name_mapping(columns):
    """
    Generates a mapping from snake_case column names to their original names.

    Parameters
    ----------
    columns : List[str]
        The list of column names to convert and map.

    Returns
    -------
    dict
        A dictionary mapping snake_case column names to their original names.
    """
    return {to_snake_case(col): col for col in columns}


def series_to_dataframe(series):
    """
    Transforms a pandas Series into a DataFrame where the columns are the index
    of the Series. If the Series' index is numeric, the index values are converted
    to strings and used as column names.

    Parameters
    ----------
    series : pandas.Series
        The Series to be transformed into a DataFrame.

    Returns
    -------
    pandas.DataFrame
        A DataFrame where each column represents a value from the Series,
        with column names corresponding to the Series' index values.

    Raises
    ------
    TypeError
        If the input is not a pandas Series.

    Examples
    --------
    >>> import pandas as pd
    >>> from fusionlab.api.util import series_to_dataframe
    >>> series = pd.Series(data=[1, 2, 3], index=["a", "b", "c"])
    >>> df = series_to_dataframe(series)
    >>> print(df)
       a  b  c
    0  1  2  3

    >>> series_numeric_index = pd.Series(data=[4, 5, 6], index=[10, 20, 30])
    >>> df_numeric = series_to_dataframe(series_numeric_index)
    >>> print(df_numeric)
      10 20 30
    0  4  5  6
    """
    if not isinstance(series, pd.Series):
        raise TypeError("Input must be a pandas Series.")
    # Convert index to string if it's numeric
    if (
        series.index.dtype.kind in "iufc"
    ):  # Checks for int, unsigned int, float, complex
        index_as_str = series.index.astype(str)
    else:
        index_as_str = series.index

    # Create a DataFrame with a single row populated with the Series' values
    # and columns named after the Series' index.
    df = pd.DataFrame([series.values], columns=index_as_str)

    return df


def get_table_size(width="auto", error="warn", return_height=False):
    """
    Determines the appropriate width (and optionally height) for table display
    based on terminal size, with options for manual width adjustment.

    Parameters
    ----------
    width : int or str, optional
        The desired width for the table. If set to 'auto', the terminal width
        is used. If an integer is provided, it will be used as the width,
        default is 'auto'.
    error : str, optional
        Error handling strategy when specified width exceeds terminal
        width: 'warn' or 'ignore'.
        Default is 'warn'.
    return_height : bool, optional
        If True, the function also returns the height of the table.
        Default is False.

    Returns
    -------
    int or tuple
        The width of the table as an integer, or a tuple of (width, height)
        if return_height is True.

    Examples
    --------
    >>> table_width = get_table_size()
    >>> print("Table width:", table_width)
    >>> table_width, table_height = get_table_size(return_height=True)
    >>> print("Table width:", table_width, "Table height:", table_height)
    """
    auto_width, auto_height = get_terminal_size()
    if width == "auto":
        width = auto_width
    else:
        try:
            width = int(width)
            if width > auto_width:
                if error == "warn":
                    warnings.warn(
                        f"Specified width {width} exceeds terminal width {auto_width}. "
                        "This may cause display issues.",
                        stacklevel=2,
                    )
        except ValueError:
            raise ValueError(
                "Width must be 'auto' or an integer; got {type(width).__name__!r}"
            )

    if return_height:
        return (width, auto_height)
    return width


def get_terminal_size():
    """
    Retrieves the current terminal size (width and height) to help dynamically
    set the maximum width for displaying data columns.

    Returns
    -------
    tuple
        A tuple containing two integers:
        - The width of the terminal in characters.
        - The height of the terminal in lines.

    Examples
    --------
    >>> from fusionlab.api.util import get_terminal_size
    >>> terminal_width, terminal_height = get_terminal_size()
    >>> print("Terminal Width:", terminal_width)
    >>> print("Terminal Height:", terminal_height)
    """
    # Use shutil.get_terminal_size if available (Python 3.3+)
    # This provides a fallback of (80, 24) which is a common default size
    if hasattr(shutil, "get_terminal_size"):
        size = shutil.get_terminal_size(fallback=(80, 24))
    else:
        # Fallback for Python versions before 3.3
        try:
            # UNIX-based systems
            size = os.popen("stty size", "r").read().split()
            return int(size[1]), int(size[0])
        except Exception:
            # Default fallback size
            size = (80, 24)
    return size.columns, size.lines


def to_camel_case(text, delimiter=None, use_regex=False):
    """
    Converts a given string to CamelCase. The function handles strings with or
    without delimiters, and can optionally use regex for splitting based on
    non-alphanumeric characters.

    Parameters
    ----------
    text : str
        The string to convert to CamelCase.
    delimiter : str, optional
        A character or string used as a delimiter to split the input string.
        Common delimiters include underscores ('_') or spaces (' '). If None
        and use_regex is ``False``, the function tries to automatically detect
        common delimiters like spaces or underscores.
        If `use_regex` is ``True``, it splits the string at any non-alphabetic
        character.
    use_regex : bool, optional
        Specifies whether to use regex for splitting the string on non-alphabetic
        characters.
        Defaults to ``False``.

    Returns
    -------
    str
        The CamelCase version of the input string.

    Examples
    --------
    >>> from fusionlab.api.util import to_camel_case
    >>> to_camel_case("outlier_results", "_")
    'OutlierResults'

    >>> to_camel_case("outlier results", " ")
    'OutlierResults'

    >>> to_camel_case("outlierresults")
    'Outlierresults'

    >>> to_camel_case("data science rocks")
    'DataScienceRocks'

    >>> to_camel_case("data_science_rocks")
    'DataScienceRocks'

    >>> to_camel_case("multi@var_analysis", use_regex=True)
    'MultiVarAnalysis'

    >>> to_camel_case("OutlierResults")
    'OutlierResults'

    >>> to_camel_case("BoxFormatter")
    'BoxFormatter'

    >>> to_camel_case("MultiFrameFormatter")
    'MultiFrameFormatter'
    """
    # Remove any leading/trailing whitespace
    text = str(text).strip()

    # Check if text is already in CamelCase and return it as is
    if text and text[0].isupper() and not text[1:].islower():
        return text

    if use_regex:
        # Split text using any non-alphabetic character as a delimiter
        words = re.split("[^a-zA-Z]", text)
    elif delimiter is None:
        if " " in text and "_" in text:
            # Both space and underscore are present, replace '_' with ' ' then split
            text = text.replace("_", " ")
            words = text.split()
        elif " " in text:
            words = text.split(" ")
        elif "_" in text:
            words = text.split("_")
        else:
            # No common delimiter found, handle as a single word
            words = [text]
    else:
        # Use the specified delimiter
        words = text.split(delimiter)

    # Capitalize the first letter of each word and join them without spaces
    # Ensure empty strings from split are ignored
    return "".join(word.capitalize() for word in words if word)


def beautify_dict(d, space=4, key=None, max_char=None):
    """
    Format a dictionary with custom indentation and aligned keys and values.

    Parameters:
    ----------
    d : dict
        The dictionary to format.
    space : int, optional
        The number of spaces to indent the dictionary entries.
    key : str, optional
        An optional key to nest the dictionary under.
    max_char : int, optional
        Maximum characters a value can have before being truncated.

    Returns:
    -------
    str
        A string representation of the dictionary with custom formatting.

    Examples:
    --------
    >>> from fusionlab.api.util import beautify_dict
    >>> dictionary = {
    ...     3: "Home & Garden",
    ...     2: "Health & Beauty",
    ...     4: "Sports",
    ...     0: "Electronics",
    ...     1: "Fashion",
    ... }
    >>> print(beautify_dict(dictionary, space=4))
    """

    if not isinstance(d, dict):
        raise TypeError(
            "Expected input to be a 'dict',"
            f" received '{type(d).__name__}' instead."
        )

    if max_char is None:
        # get it automatically
        max_char, _ = get_terminal_size()
    # Determine the longest key for alignment
    if len(d) == 0:
        max_key_length = 0
    else:
        max_key_length = max(len(str(k)) for k in d.keys())

    # Create a list of formatted rows
    formatted_rows = []
    for dkey, value in sorted(d.items()):
        value_str = str(value)
        if max_char is not None and len(value_str) > max_char:
            value_str = value_str[:max_char] + "..."
        # Ensure all keys are right-aligned to the longest key length
        formatted_row = f"{str(dkey):>{max_key_length}}: '{value_str}'"
        formatted_rows.append(formatted_row)

    # Join all rows into a single string with custom indentation
    indent = " " * space
    inner_join = ",\n" + indent
    formatted_dict = "{\n" + indent + inner_join.join(formatted_rows) + "\n}"

    if key:
        # Prepare outer key indentation and format
        # Slightly less than the main indent
        outer_indent = " " * (space - 2 + len(key) + max_key_length)
        # Construct a new header with the key
        formatted_dict = f"{key} : {formatted_dict}"
        # Split lines and indent properly to align with the key
        lines = formatted_dict.split("\n")
        for i in range(1, len(lines)):
            lines[i] = outer_indent + lines[i]
            if max_char is not None and len(lines[i]) > max_char:
                lines[i] = lines[i][:max_char] + "..."
        # format lins -1
        # lines [-1] = outer_indent + lines [-1]
        formatted_dict = "\n".join(lines)

    return formatted_dict


def remove_extra_spaces(text):
    """
    Removes extra spaces from the input text, leaving only one
    space between words.

    Parameters
    ----------
    text : str
        The string from which extra spaces need to be removed.

    Returns
    -------
    str
        The string with only one space between each word.

    Example
    -------
    >>> from fusionlab.api.util import remove_extra_spaces
    >>> text = "this is      text that    have   extra          space."
    >>> remove_extra_spaces(text)
    'this is text that have extra space.'
    """
    # Use regular expression to replace multiple spaces with a single space
    cleaned_text = re.sub(r"\s+", " ", text).strip()
    return cleaned_text


def format_iterable(attr):
    """
    Formats an iterable with a string representation that includes
    statistical or structural information depending on the iterable's type.
    """

    def _numeric_stats(iterable):
        return {
            "min": round(np.min(iterable), 4),
            "max": round(np.max(iterable), 4),
            "mean": round(np.mean(iterable), 4),
            "len": len(iterable),
        }

    def _format_numeric_iterable(iterable):
        stats = _numeric_stats(iterable)
        return (
            f"{type(iterable).__name__} (min={stats['min']},"
            f" max={stats['max']}, mean={stats['mean']}, len={stats['len']})"
        )

    def _format_ndarray(array):
        stats = (
            _numeric_stats(array.flat)
            if np.issubdtype(array.dtype, np.number)
            else {}
        )
        details = ", ".join([f"{key}={value}" for key, value in stats.items()])
        return f"ndarray ({details}, shape={array.shape}, dtype={array.dtype})"

    def _format_pandas_object(obj):
        if isinstance(obj, pd.Series):
            stats = _numeric_stats(obj) if obj.dtype != "object" else {}
            details = ", ".join(
                [f"{key}={value}" for key, value in stats.items()]
            )
            if details:
                details += ", "
            return f"Series ({details}len={obj.size}, dtype={obj.dtype})"
        elif isinstance(obj, pd.DataFrame):
            numeric_cols = obj.select_dtypes(include=np.number).columns
            stats = (
                _numeric_stats(obj[numeric_cols].values.flat)
                if not numeric_cols.empty
                else {}
            )
            details = ", ".join(
                [f"{key}={value}" for key, value in stats.items()]
            )
            if details:
                details += ", "
            return (
                f"DataFrame ({details}n_rows={obj.shape[0]},"
                f" n_cols={obj.shape[1]}, dtypes={obj.dtypes.unique()})"
            )

    if isinstance(attr, (list, tuple, set)) and all(
        isinstance(item, (int, float)) for item in attr
    ):
        return _format_numeric_iterable(attr)
    elif isinstance(attr, np.ndarray):
        return _format_ndarray(attr)
    elif isinstance(attr, (pd.Series, pd.DataFrame)):
        return _format_pandas_object(attr)

    return str(attr)


def format_dict_result(
    dictionary, dict_name="Container", max_char=50, include_message=False
):
    """
    Formats a dictionary into a string with specified formatting rules.

    Parameters
    ----------
    dictionary : dict
        The dictionary to format.
    dict_name : str, optional
        The name of the dictionary, by default 'Container'.
    max_char : int, optional
        The maximum number of characters for each value before truncating,
        by default 50.
    include_message : bool, optional
        Whether to include a remainder message at the end, by default False.

    Returns
    -------
    str
        The formatted string representation of the dictionary.

    Examples
    --------
    >>> example_dict = {
    ...     "key1": "short value",
    ...     "key2": "a much longer value that should be truncated for readability purposes",
    ...     "key3": "another short value",
    ...     "key4": "value",
    ... }
    >>> print(
    ...     format_dict_result(
    ...         example_dict, dict_name="ExampleDict", max_char=30
    ...     )
    ... )
    ExampleDict({
        key1: short value,
        key2: a much longer value that s...,
        key3: another short value,
        key4: value,
    })

    Notes
    -----
    The function calculates the required indentation based on the length of the
    dictionary name and the maximum key length. If a value exceeds the specified
    maximum length, it truncates the value and appends an ellipsis ("...").
    """
    max_key_length = max(len(str(key)) for key in dictionary.keys())
    formatted_lines = [f"{dict_name}({{"]

    for key, value in dictionary.items():
        if (
            isinstance(value, value.__class__)
            and not hasattr(value, "__array__")
            and not isinstance(value, (str, list, tuple))
        ):
            try:
                formatted_value = value.__class__.__name__
            except:
                formatted_value = value.__name__
        else:
            formatted_value = format_iterable(value)

        if len(formatted_value) > max_char:
            formatted_value = formatted_value[: max_char - 3] + "..."
        formatted_lines.append(
            f"{' ' * (len(dict_name) + 2)}{key:{max_key_length}}: {formatted_value},"
        )

    formatted_lines.append(" " * (len(dict_name) + 1) + "})")

    remainder = f"[Use <{dict_name}.key> to get the full value ...]"

    return (
        "\n".join(formatted_lines) + f"\n\n{remainder}"
        if include_message
        else "\n".join(formatted_lines)
    )


def count_functions(
    module_name,
    include_class=False,
    return_counts=True,
    include_private=False,
    include_local=False,
):
    """
    Count and list the number of functions and classes in a specified module.

    Parameters
    ----------
    module_name : str
        The name of the module to inspect, in the format `package.module`.
    include_class : bool, optional
        Whether to include classes in the count and listing. Default is
        `False`.
    return_counts : bool, optional
        Whether to return only the count of functions and classes (if
        ``include_class`` is `True`). If `False`, returns a list of functions
        and classes in alphabetical order. Default is `True`.
    include_private : bool, optional
        Whether to include private functions and classes (those starting with
        `_`). Default is `False`.
    include_local : bool, optional
        Whether to include local (nested) functions in the count and listing.
        Default is `False`.

    Returns
    -------
    int or list
        If ``return_counts`` is `True`, returns the count of functions and
        classes (if ``include_class`` is `True`). If ``return_counts`` is
        `False`, returns a list of function and class names (if
        ``include_class`` is `True`) in alphabetical order.

    Notes
    -----
    This function dynamically imports the specified module and analyzes its
    Abstract Syntax Tree (AST) to count and list functions and classes. It
    provides flexibility to include or exclude private and local functions
    based on the parameters provided.

    The process can be summarized as:

    .. math::
        \text{total\\_count} =
        \text{len(functions)} + \text{len(classes)}

    where:

    - :math:`\text{functions}` is the list of functions found in the module.
    - :math:`\text{classes}` is the list of classes found in the module
      (if ``include_class`` is `True`).

    Examples
    --------
    >>> from fusionlab.api.util import count_functions_classes
    >>> count_functions_classes('fusionlab.api.util', include_class=True,
                                return_counts=True)
    10

    >>> count_functions('fusionlab.api.util', include_class=True,
                                return_counts=False)
    ['ClassA', 'ClassB', 'func1', 'func2', 'func3']

    >>> count_functions('fusionlab.api.util', include_class=False,
                                return_counts=True, include_private=True)
    15

    >>> count_functions('fusionlab.api.util', include_class=False,
                                return_counts=False, include_private=True)
    ['_private_func1', '_private_func2', 'func1', 'func2']

    See Also
    --------
    ast : Abstract Syntax Tree (AST) module for parsing Python source code.

    References
    ----------
    .. [1] Python Software Foundation. Python Language Reference, version 3.9.
       Available at http://www.python.org
    .. [2] Python `ast` module documentation. Available at
       https://docs.python.org/3/library/ast.html
    """

    try:
        import ast
    except ImportError as e:  # Catch the specific ImportError exception
        raise ImportError(
            "The 'ast' module could not be imported. This module is essential"
            " for analyzing Python source code to count functions and classes."
            " Ensure that you are using a standard Python distribution, which"
            " includes the 'ast' module by default."
        ) from e

    import importlib
    import inspect

    # Import the module dynamically
    module = importlib.import_module(module_name)

    # Get the source code of the module
    source = inspect.getsource(module)

    # Parse the source code into an AST
    tree = ast.parse(source)

    # Initialize lists to store function and class names
    functions = []
    classes = []

    def is_local_function(node):
        """Determine if the function is local (nested)."""
        while node:
            if isinstance(node, ast.FunctionDef):
                return True
            node = getattr(node, "parent", None)
        return False

    # Add parent references to each node
    for node in ast.walk(tree):
        for child in ast.iter_child_nodes(node):
            child.parent = node

    # Traverse the AST to find function and class definitions
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef):
            if (include_private or not node.name.startswith("_")) and (
                include_local or not is_local_function(node.parent)
            ):
                functions.append(node.name)
        elif isinstance(node, ast.ClassDef) and include_class:
            if include_private or not node.name.startswith("_"):
                classes.append(node.name)

    # Combine and sort the lists if needed
    if include_class:
        result = sorted(functions + classes)
    else:
        result = sorted(functions)

    if return_counts:
        return len(result)
    else:
        return result


def round_numeric_values(df, precision=4):
    """
    Rounds numeric floating values in a DataFrame to the specified
    precision. Integers and non-numeric values are not modified.

    Parameters
    ----------
    df : pandas.DataFrame
        The DataFrame whose numeric floating values are to be rounded.
        The DataFrame can contain mixed data types, including integers,
        floating-point numbers, and non-numeric values.
    precision : int, optional
        The number of decimal places to round to. Defaults to 4.
        - If `precision` is set to a positive integer, it specifies
          the number of decimal places.
        - If `precision` is set to a negative integer, it specifies
          the number of positions to the left of the decimal point.

    Returns
    -------
    pandas.DataFrame
        A DataFrame with numeric floating values rounded to the
        specified precision. Integers and non-numeric values remain
        unchanged.

    Raises
    ------
    ValueError
        If `precision` is not an integer.

    Notes
    -----
    This function applies rounding only to floating-point numbers.
    Integers and non-numeric values remain unchanged. The function
    uses pandas' `applymap` to apply the rounding operation element-wise
    across the DataFrame.

    The mathematical formulation for the rounding operation is given by:

    .. math::
        y =
        \begin{cases}
        \text{round}(x, \text{precision}) & \text{if } x \\in \\mathbb{R} \\
        x & \text{otherwise} \\
        \\end{cases}

    Where:
    - :math:`x` is the input value
    - :math:`y` is the output value after rounding
    - :math:`\text{precision}` is the specified number of decimal places

    Examples
    --------
    >>> from fusionlab.api.util import round_numeric_values
    >>> import pandas as pd
    >>> df = pd.DataFrame({
    ...     'A': [1.12345, 2.6789, 3],
    ...     'B': [4, 5.98765, 'text'],
    ...     'C': [7.123456, 8.654321, 9.0]
    ... })
    >>> round_numeric_values(df, precision=2)
           A      B     C
    0  1.12      4  7.12
    1  2.68  5.99  8.65
    2     3  text     9

    See Also
    --------
    pandas.DataFrame.applymap : Element-wise operation on DataFrame.

    References
    ----------
    .. [1] "NumPy Documentation", https://numpy.org/doc/stable/

    """

    def round_if_float(x):
        if isinstance(x, float):
            return round(x, precision)
        return x

    return df.applymap(round_if_float)
