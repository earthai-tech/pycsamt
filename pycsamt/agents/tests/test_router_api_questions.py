# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
import pytest

from pycsamt.agents.router import classify_intent_offline


@pytest.mark.parametrize("query,expected", [
    ("estimate_ss_ama parameters and return columns", "question"),
    ("ensure_sites arguments", "question"),
    ("correct_ss_ama signature", "question"),
    ("write a script showing estimate_ss_ama parameters", "code"),
    ("run static shift on L22PLT with these parameters", "workflow"),
    ("what workflows are supported?", "meta"),
    ("correct galvanic static shift with the AMA method", "workflow"),
])
def test_api_question_and_action_boundaries(query, expected):
    assert classify_intent_offline(query)[0] == expected
