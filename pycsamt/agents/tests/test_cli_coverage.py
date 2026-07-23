"""Command-dispatch coverage for ``python -m pycsamt.agents``."""

from __future__ import annotations

import builtins

import pytest

from pycsamt.agents import __main__ as cli


def test_catalogue_pricing_help_and_unknown(capsys):
    cli.main([])
    cli.main(["list"])
    cli.main(["pricing"])
    text = capsys.readouterr().out
    assert "Usage" in text
    assert "full catalogue" in text
    assert "LLM pricing" in text
    with pytest.raises(SystemExit):
        cli.main(["unknown"])
    with pytest.raises(SystemExit):
        cli.main(["preview"])


def test_web_dispatch(monkeypatch, capsys):
    import pycsamt.agents.web as web

    called = {}
    monkeypatch.setattr(web, "launch", lambda **kw: called.update(kw))
    cli.main(["web", "--share", "--port=8123"])
    assert called == {"share": True, "server_port": 8123}
    assert "8123" in capsys.readouterr().out


def test_zoo_list_download_and_errors(monkeypatch, capsys):
    import pycsamt.ai._zoo as zoo

    monkeypatch.setattr(zoo, "list_pretrained", lambda: {"tiny": "Tiny model"})
    monkeypatch.setattr(
        zoo, "download_checkpoint", lambda *a, **k: "/tmp/tiny.ckpt"
    )
    cli.main(["zoo"])
    cli.main(["zoo", "tiny", "--force"])
    assert "Tiny model" in capsys.readouterr().out

    monkeypatch.setattr(
        zoo, "download_checkpoint",
        lambda *a, **k: (_ for _ in ()).throw(KeyError("missing")),
    )
    with pytest.raises(SystemExit):
        cli.main(["zoo", "missing"])
    monkeypatch.setattr(
        zoo, "download_checkpoint",
        lambda *a, **k: (_ for _ in ()).throw(RuntimeError("offline")),
    )
    with pytest.raises(SystemExit):
        cli.main(["zoo", "tiny"])


def test_console_print_unicode_fallback(monkeypatch):
    calls = []

    def flaky(*args, **kwargs):
        if not calls:
            calls.append(args)
            raise UnicodeEncodeError("ascii", "é", 0, 1, "bad")
        calls.append(args)

    monkeypatch.setattr(builtins, "print", flaky)
    cli.print("é")
    assert len(calls) == 2
