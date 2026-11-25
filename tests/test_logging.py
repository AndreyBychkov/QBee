import pandas as pd
from qbee import functions, quadratize, QBEE_CONFIG


def test_logging_writes_the_search_as_strings(tmp_path):
    x = functions("x")
    log = tmp_path / "log.feather"
    enabled, file = QBEE_CONFIG.logging_enabled, QBEE_CONFIG.logging_file
    QBEE_CONFIG.logging_enabled, QBEE_CONFIG.logging_file = True, str(log)
    try:
        assert quadratize([(x, x ** 3)]) is not None
    finally:
        QBEE_CONFIG.logging_enabled, QBEE_CONFIG.logging_file = enabled, file
    df = pd.read_feather(log)
    assert list(df.columns) == ["from", "to"]
    assert len(df) > 0
    assert all(isinstance(v, str) for col in ("from", "to") for v in df[col])
