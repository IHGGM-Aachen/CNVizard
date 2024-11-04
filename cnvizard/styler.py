"""
File which contains the pandas stylers used for the CNVizard
@author: Jeremias Krause, Carlos Classen, Matthias Begemann, Florian Kraft
@company: UKA Aachen (RWTH)
@mail: jerkrause@ukaachen.de
"""

from pandas.io.formats.style import Styler


def make_pretty(styler: Styler) -> Styler:
    """
    Styler function used to highlight and format the filtered .cnr/bintest DataFrame
    """
    # Set a precision of 2 for float values
    styler.format(precision=2)
    # Set a precision of 0 for the OMIMG column
    styler.format(subset="OMIMG", precision=0)
    # Set the thousand separator as ',' and apply it on the columns start and
    # end
    styler.format(subset=["start", "end"], thousands=",", decimal=".")
    # Highlight values in the weight column in green which are greater or
    # equal than 0.9
    styler.highlight_between(subset="weight", color="green", axis=0, left=0.9, right=1)
    # Highlight values in the weight column in yellow which are greater or
    # equal to 0.8 but smaller than 1
    styler.highlight_between(
        subset="weight", color="yellow", axis=0, left=0.8, right=0.9
    )
    # Highlight values in the weight column in red which are smaller than 0.8
    styler.highlight_between(subset="weight", color="red", axis=0, left=0, right=0.8)
    # Highlight values in the depth column in red which are 0
    styler.highlight_between(subset="depth", color="red", axis=0, right=0)
    # Highlight values in the log2 column in yellow which are smaller or equal
    # to -0.65
    styler.highlight_between(subset="log2", axis=0, right=-0.65)
    return styler


def mark_log2(value):
    if value <= float(0.65):
        return "background-color: pink;"
