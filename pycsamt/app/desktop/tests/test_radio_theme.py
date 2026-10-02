"""Selected radio indicators remain visible under both custom themes."""
from pathlib import Path

import pytest

pytest.importorskip("PySide6")
from PySide6.QtWidgets import QRadioButton, QStyle, QStyleOptionButton


@pytest.mark.parametrize("theme,accent", [
    ("light", "#1e66f5"), ("dark", "#89b4fa"),
])
def test_radio_indicator_paints_checked_dot(qapp, theme, accent):
    from PySide6.QtGui import QColor

    path = (Path(__file__).resolve().parents[1] / "resources"
            / f"{theme}_theme.qss")
    radio = QRadioButton("Period")
    radio.setStyleSheet(path.read_text(encoding="utf-8"))
    radio.resize(180, 40)
    radio.show()
    radio.setChecked(True)
    qapp.processEvents()
    option = QStyleOptionButton()
    radio.initStyleOption(option)
    center = radio.style().subElementRect(
        QStyle.SubElement.SE_RadioButtonIndicator, option, radio
    ).center()
    pixmap = radio.grab()
    scale = pixmap.devicePixelRatio()
    selected = pixmap.toImage().pixelColor(
        int(center.x()*scale), int(center.y()*scale)
    )
    assert selected == QColor(accent)
    radio.setAutoExclusive(False)
    radio.setChecked(False)
    assert radio.grab().toImage().pixelColor(
        int(center.x()*scale), int(center.y()*scale)
    ) != selected
    radio.close()
