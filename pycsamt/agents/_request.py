# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Cooperative cancellation shared by all Agent Master providers and stages."""
from __future__ import annotations

import time
from contextlib import contextmanager
from contextvars import ContextVar
from functools import wraps


class RequestCancelled(Exception):
    """The owning conversation has stopped this request."""


_CANCELLED = ContextVar('agent_request_cancelled', default=None)


@contextmanager
def request_scope(cancelled):
    token = _CANCELLED.set(cancelled)
    try:
        checkpoint()
        yield
    finally:
        _CANCELLED.reset(token)


def is_cancelled():
    predicate = _CANCELLED.get()
    return bool(predicate is not None and predicate())


def checkpoint():
    if is_cancelled():
        raise RequestCancelled('Task stopped by user; no further stages will run.')


def request_stage(function):
    """Check before and after a synchronous operation; never kill its thread."""
    @wraps(function)
    def checked(*args, **kwargs):
        checkpoint()
        result = function(*args, **kwargs)
        checkpoint()
        return result
    return checked


def cancellable_sleep(seconds):
    """Interrupt retry backoff without forcibly killing worker threads."""
    if _CANCELLED.get() is None:
        time.sleep(seconds)
        return
    end = time.monotonic() + seconds
    while time.monotonic() < end:
        checkpoint()
        time.sleep(min(0.1, max(0, end - time.monotonic())))
    checkpoint()
