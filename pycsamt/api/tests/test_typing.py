from __future__ import annotations

from pycsamt.api.typing import (
    ArrayLike,
    DType,
    NDArray,
    PathLike,
    Shape,
    is_path_like,
)


def test_compat_alias_subclasses_remain_subscriptable_at_runtime():
    # Legacy annotations like ArrayLike[DType[float]] must still evaluate.
    assert ArrayLike[float] is ArrayLike
    assert ArrayLike[float, DType[float]] is ArrayLike
    assert Shape[int] is Shape
    assert NDArray[float] is NDArray


def test_is_path_like_true_for_str_bytes_and_path():
    from pathlib import Path

    assert is_path_like("some/path")
    assert is_path_like(b"some/path")
    assert is_path_like(Path("some/path"))


def test_is_path_like_false_for_non_path_values():
    assert not is_path_like(42)
    assert not is_path_like(None)
    assert not is_path_like(["not", "a", "path"])


def test_path_like_alias_is_importable():
    assert PathLike is not None
