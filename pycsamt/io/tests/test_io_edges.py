from __future__ import annotations

import gzip
import io

import pytest

from pycsamt.io import parsers
from pycsamt.io.formats import (
    TransferFunctionFormatError,
    _read_edi,
    _source_probe,
    detect_tf_format,
    get_tf_format,
    get_tf_format_for_target,
    list_tf_formats,
    register_tf_format,
)
from pycsamt.io.model_formats import (
    ModelFormatError,
    _is_pcsm,
    _is_pcsf,
    _pcsm_head,
    detect_model_format,
    get_model_format,
    get_model_format_for_target,
    list_model_formats,
    register_model_format,
)
from pycsamt.io.transfer import (
    _buffer_nonseekable,
    read_transfer_function,
    write_transfer_function,
)


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"reader": None}, "reader must be callable"),
        ({"reader": lambda source: source, "writer": 1}, "writer must be callable"),
        ({"reader": lambda source: source, "detector": 1}, "detector must be callable"),
    ],
)
def test_tf_registration_validates_callbacks(kwargs, message):
    with pytest.raises(TypeError, match=message):
        register_tf_format("bad_callback", **kwargs)


def test_tf_registration_normalizes_public_metadata_and_alias():
    spec = register_tf_format(
        " Edge Format ",
        reader=lambda source: source,
        extensions=("EDGE", ".edge", ""),
        aliases=(" Edge Alias ",),
        description=None,
        replace=True,
    )
    assert spec.name == "edge_format"
    assert spec.extensions == (".edge", ".edge")
    assert get_tf_format("edge alias") is spec
    assert list_tf_formats()["edge_format"]["writable"] is False
    with pytest.raises(TransferFunctionFormatError, match="unknown"):
        get_tf_format("missing-format")


def test_source_probe_handles_bytes_paths_inline_and_seekable_streams(tmp_path):
    assert _source_probe(b"abc") == "abc"
    assert _source_probe("<EM_TF/>") == "<EM_TF/>"
    path = tmp_path / "probe.txt"
    path.write_bytes(b"payload")
    assert _source_probe(path) == "payload"
    assert _source_probe(str(path)) == "payload"

    stream = io.StringIO("abcdef")
    stream.seek(2)
    assert _source_probe(stream, limit=2) == "cd"
    assert stream.tell() == 2
    with pytest.raises(TypeError, match="source must be"):
        _source_probe(object())


def test_source_probe_stream_without_tell_seek_returns_decoded_bytes():
    class _ReadOnlyBytesStream:
        def read(self, limit):
            return b"raw-bytes"

    assert _source_probe(_ReadOnlyBytesStream()) == "raw-bytes"


def test_source_probe_stream_with_failing_tell_and_seek_is_tolerated():
    class _FlakyStream:
        def __init__(self):
            self.seek_called_with = None

        def tell(self):
            raise OSError("tell not supported")

        def seek(self, position):
            self.seek_called_with = position
            raise OSError("seek not supported")

        def read(self, limit):
            return "text"

    stream = _FlakyStream()
    assert _source_probe(stream) == "text"
    # position was None (tell() failed), so seek() must never be called
    assert stream.seek_called_with is None


def test_source_probe_restores_position_even_when_seek_back_fails():
    class _NoRewindStream:
        def __init__(self):
            self.seek_calls = []

        def tell(self):
            return 5

        def seek(self, position):
            self.seek_calls.append(position)
            raise OSError("cannot rewind")

        def read(self, limit):
            return "abc"

    stream = _NoRewindStream()
    assert _source_probe(stream) == "abc"
    assert stream.seek_calls == [5]


def test_tf_detection_and_target_errors_include_useful_hints(tmp_path):
    assert detect_tf_format("<?xml version='1.0'?><x:EM_TF xmlns:x='u'/>") == "emtf_xml"
    assert detect_tf_format("# comment\n>HEAD\n>=DEFINEMEAS\n") == "edi"
    with pytest.raises(TransferFunctionFormatError, match="requires an explicit"):
        get_tf_format_for_target(io.BytesIO())
    with pytest.raises(TransferFunctionFormatError, match="no extension"):
        get_tf_format_for_target(tmp_path / "output")
    with pytest.raises(TransferFunctionFormatError, match="no writable"):
        get_tf_format_for_target(tmp_path / "output.unknown")
    with pytest.raises(TransferFunctionFormatError, match="suggests emtf_xml"):
        detect_tf_format(tmp_path / "invalid.xml")


def test_edi_reader_rejects_in_memory_sources():
    for source in (b">HEAD", io.StringIO(">HEAD"), ">HEAD\n>FREQ"):
        with pytest.raises(TypeError, match="filesystem path"):
            _read_edi(source)


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"reader": None}, "reader must be callable"),
        ({"reader": lambda source: source, "writer": 1}, "writer must be callable"),
        ({"reader": lambda source: source, "detector": 1}, "detector must be callable"),
    ],
)
def test_model_registration_validates_callbacks(kwargs, message):
    with pytest.raises(TypeError, match=message):
        register_model_format("bad_model_callback", **kwargs)


def test_model_registration_listing_alias_and_target_errors(tmp_path):
    spec = register_model_format(
        "Edge Model",
        reader=lambda source: source,
        aliases=("EMODEL",),
        extensions=(".edge-model",),
        replace=True,
    )
    assert get_model_format("emodel") is spec
    assert list_model_formats()["edge_model"]["writable"] is False
    with pytest.raises(ModelFormatError, match="unknown"):
        get_model_format("missing-model")
    with pytest.raises(ModelFormatError, match="requires an explicit"):
        get_model_format_for_target(io.BytesIO())
    with pytest.raises(ModelFormatError, match="no extension"):
        get_model_format_for_target(tmp_path / "model")
    with pytest.raises(ModelFormatError, match="no writable"):
        get_model_format_for_target(tmp_path / "model.unknown")
    with pytest.raises(ModelFormatError, match="suggests pcsm"):
        detect_model_format(tmp_path / "invalid.pcsm")


def test_pcsm_content_detection_plain_gzip_and_invalid_sources(tmp_path):
    plain = tmp_path / "model.txt"
    plain.write_text("comment\nPCSM_VERSION 1.0\n", encoding="utf-8")
    assert _pcsm_head(plain).startswith("comment")
    assert _is_pcsm(plain)
    assert detect_model_format(plain) == "pcsm"

    compressed = tmp_path / "model.bin"
    with gzip.open(compressed, "wt", encoding="utf-8") as stream:
        stream.write("PCSM_VERSION 1.0\n")
    assert _is_pcsm(compressed)
    assert not _is_pcsm(io.BytesIO())
    assert not _is_pcsm(tmp_path / "missing.pcsm")


def test_pcsf_detector_rejects_missing_and_non_hdf_files(tmp_path):
    assert not _is_pcsf(tmp_path / "missing.pcsf")
    path = tmp_path / "fake.pcsf"
    path.write_bytes(b"not hdf5")
    assert not _is_pcsf(path)


def test_pcsf_detector_rejects_hdf5_signature_with_corrupt_body(tmp_path):
    # Real HDF5 magic bytes, but not a well-formed HDF5 file: h5py.File
    # raises OSError while opening it, hit by _is_pcsf's except clause.
    path = tmp_path / "corrupt.pcsf"
    path.write_bytes(b"\x89HDF\r\n\x1a\n" + b"\x00" * 32)
    assert not _is_pcsf(path)


def test_is_pcsm_returns_false_when_no_version_header_present(tmp_path):
    path = tmp_path / "no_header.txt"
    path.write_text("just some unrelated text\nno marker here\n", encoding="utf-8")
    assert not _is_pcsm(path)


def test_generic_parser_dispatch_and_errors(monkeypatch, tmp_path):
    calls = []

    class FakeConfig:
        parsers = {".csv": lambda path, **kw: (path, kw)}

        @staticmethod
        def writers(obj):
            return {".csv": lambda path, **kw: calls.append((obj, path, kw))}

    monkeypatch.setattr(parsers, "Config", FakeConfig)
    assert parsers.read_any("data.CSV", sep=";") == ("data.CSV", {"sep": ";"})
    obj = object()
    assert parsers.write_any(obj, str(tmp_path / "out.csv"), index=False) is None
    assert calls[0][0] is obj and calls[0][2] == {"index": False}
    with pytest.raises(ValueError, match="No parser configured"):
        parsers.read_any("data.unknown")
    with pytest.raises(ValueError, match="No writer configured"):
        parsers.write_any(obj, "data.unknown")


class _NonSeekable:
    def __init__(self, value):
        self.value = value

    def read(self, *args):
        return self.value

    def seekable(self):
        return False


class _SeekableProbe:
    def __init__(self, *, seekable_result=None, seekable_error=False, io_error=False):
        self.seekable_result = seekable_result
        self.seekable_error = seekable_error
        self.io_error = io_error
        self.read_called = False

    def read(self, *args):
        self.read_called = True
        return b"buffered"

    def seekable(self):
        if self.seekable_error:
            raise OSError("unknown")
        return self.seekable_result

    def tell(self):
        if self.io_error:
            raise OSError("tell failed")
        return 3

    def seek(self, position):
        if self.io_error:
            raise OSError("seek failed")


def test_transfer_buffer_and_explicit_read_dispatch():
    source = _NonSeekable("payload")
    assert _buffer_nonseekable(source) == "payload"
    assert _buffer_nonseekable("plain") == "plain"

    register_tf_format(
        "edge_reader",
        reader=lambda value, **kw: (value, kw),
        aliases=("er",),
        replace=True,
    )
    assert read_transfer_function("payload", format="er", flag=True) == (
        "payload",
        {"flag": True},
    )


def test_transfer_buffer_preserves_seekable_streams_and_handles_probe_errors():
    directly_seekable = _SeekableProbe(seekable_result=True)
    assert _buffer_nonseekable(directly_seekable) is directly_seekable
    assert not directly_seekable.read_called

    seek_and_tell_work = _SeekableProbe(seekable_result=False)
    assert _buffer_nonseekable(seek_and_tell_work) is seek_and_tell_work
    assert not seek_and_tell_work.read_called

    fallback = _SeekableProbe(seekable_error=True, io_error=True)
    assert _buffer_nonseekable(fallback) == b"buffered"
    assert fallback.read_called


def test_transfer_writer_rejects_read_only_format(tmp_path):
    register_tf_format(
        "read_only_edge",
        reader=lambda value: value,
        replace=True,
    )
    with pytest.raises(TransferFunctionFormatError, match="not writable"):
        write_transfer_function(
            object(), tmp_path / "out.edge", format="read_only_edge"
        )
