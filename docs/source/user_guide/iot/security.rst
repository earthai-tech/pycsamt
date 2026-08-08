.. _user_guide_iot_security:

Security
========

The IoT security helpers define :term:`transport security` and
:term:`authentication` configuration for telemetry clients. They do not
implement cryptography directly. :term:`TLS`, certificates, and
:term:`authentication header`\ s are passed to the concrete transport
layer, while the configuration objects keep logs, manifests, and
notebooks from leaking :term:`secret` values.

Use this page for security configuration, :term:`redaction`,
environment loading, and validation. The examples use placeholder
certificate paths and tokens; replace them with deployment-specific
secrets outside version control. A useful way to read the helpers is as
a reproducible mapping from configuration state to two products: a
redacted record safe for a :term:`provenance manifest`, and a live
option dictionary for the transport.

TLS And Bearer Authentication
-----------------------------

Use :class:`~pycsamt.iot.security.TLSConfig` for :term:`TLS material` and
:class:`~pycsamt.iot.security.Credential` for authentication. A
:class:`~pycsamt.iot.security.SecurityConfig` combines both and can
restrict allowed :term:`telemetry protocol`\ s. In notation, let
:math:`s=(t,c,r,P)` be the security state: TLS configuration :math:`t`,
credential :math:`c`, the boolean :math:`r` for ``require_tls``, and the
optional allowed-protocol set :math:`P`. The serialisable summary
applies a redaction operator :math:`R(x)=\mathrm{None}` when :math:`x`
is missing or empty, otherwise :math:`R(x)=\mathtt{***REDACTED***}`.
Paths to certificates are retained because they identify files, while
token, password, and API-key values are replaced by :math:`R`.

.. code-block:: pycon

   >>> import json
   >>> from pycsamt.iot import AuthScheme, Credential, SecurityConfig, TLSConfig
   >>> security = SecurityConfig(
   ...     tls=TLSConfig(
   ...         enabled=True, ca_cert="certs/field-ca.pem",
   ...         certfile="certs/node-client.pem", keyfile="certs/node-client.key",
   ...         verify=True, min_version="TLSv1.2",
   ...     ),
   ...     credential=Credential(scheme=AuthScheme.BEARER, token="field-secret-token"),
   ...     require_tls=True, allowed_protocols=["https", "mqtt", "websocket"],
   ... )
   >>> print(json.dumps(security.as_dict(), indent=2))
   {
     "tls": {
       "enabled": true,
       "ca_cert": "certs/field-ca.pem",
       "certfile": "certs/node-client.pem",
       "keyfile": "certs/node-client.key",
       "verify": true,
       "min_version": "TLSv1.2"
     },
     "credential": {
       "scheme": "bearer",
       "token": "***REDACTED***",
       "username": null,
       "password": null,
       "api_key": null,
       "api_key_header": "X-API-Key"
     },
     "require_tls": true,
     "allowed_protocols": [
       "https",
       "mqtt",
       "websocket"
     ]
   }

Every certificate path in ``tls`` survives redaction intact — a
reviewer needs to know *which* CA and client certificate a node was
configured to trust — while ``token`` is the only credential field that
was actually secret here, and it is the only one replaced by
:math:`R`. ``username``, ``api_key``, and ``password`` all print as
``null`` simply because this credential is a bearer scheme and never
set them, not because they were redacted.

Build Client Options
--------------------

:meth:`~pycsamt.iot.security.SecurityConfig.client_options` returns the
live options that can be passed to telemetry clients. These options
include raw secrets because the transport needs them. Redact them before
printing. The mapping is deliberately asymmetric: ``as_dict()`` computes
a public record :math:`M_{\mathrm{public}}(s)` with :math:`R(c)`, but
``client_options()`` computes a transport record
:math:`M_{\mathrm{live}}(s)` that contains the real credential material.
If TLS is enabled, :math:`M_{\mathrm{live}}` contains ``tls=True`` and
any available certificate paths. The credential then contributes either
an ``Authorization`` header, an API-key header, or username/password
fields according to the selected scheme.

.. code-block:: pycon

   >>> from pycsamt.iot import redact_secret
   >>> options = security.client_options()
   >>> print(json.dumps(
   ...     {
   ...         "tls": options.get("tls"), "ca_cert": options.get("ca_cert"),
   ...         "certfile": options.get("certfile"), "keyfile": options.get("keyfile"),
   ...         "headers": {
   ...             key: redact_secret(value)
   ...             for key, value in options.get("headers", {}).items()
   ...         },
   ...         "token": redact_secret(options.get("token")),
   ...     },
   ...     indent=2,
   ... ))
   {
     "tls": true,
     "ca_cert": "certs/field-ca.pem",
     "certfile": "certs/node-client.pem",
     "keyfile": "certs/node-client.key",
     "headers": {
       "Authorization": "***REDACTED***"
     },
     "token": "***REDACTED***"
   }

Redacting this print call was a choice made at the call site, not by
``client_options()`` itself — the real dictionary it returned a moment
earlier has ``headers["Authorization"] = "Bearer field-secret-token"``
and ``token = "field-secret-token"`` in plain text, exactly what a
transport client needs to authenticate. Never log ``options`` directly;
always pass it through :func:`~pycsamt.iot.security.redact_secret` (or
skip printing it at all) the way this example does.

Other Credential Schemes
------------------------

:term:`Bearer token`\ s, :term:`API key`\ s, and :term:`basic
authentication` are supported. Use ``headers()`` only when feeding a
live transport. Use ``as_dict()`` for logs, manifests, and notebooks.
Reproducibly, the header function is a small piecewise mapping. For a
bearer credential with token :math:`b`, it returns ``Authorization:
Bearer b``. For basic authentication with username :math:`u` and
password :math:`p`, it first encodes the byte string :math:`u:p` with
Base64 and then returns ``Authorization: Basic <encoded>``. For an API
key :math:`k`, it returns ``<api_key_header>: k``. Redaction is not part
of ``headers()`` because the result is meant for the live transport
boundary.

.. code-block:: pycon

   >>> api_key = Credential(scheme="api_key", api_key="api-key-123", api_key_header="X-Field-Key")
   >>> basic = Credential(scheme="basic", username="operator", password="station-password")
   >>> print(json.dumps(
   ...     {
   ...         "api_key_redacted": {
   ...             key: redact_secret(value) for key, value in api_key.headers().items()
   ...         },
   ...         "basic_header_prefix": basic.headers()["Authorization"][:16],
   ...         "basic_summary": basic.as_dict(),
   ...     },
   ...     indent=2,
   ... ))
   {
     "api_key_redacted": {
       "X-Field-Key": "***REDACTED***"
     },
     "basic_header_prefix": "Basic b3BlcmF0b3",
     "basic_summary": {
       "scheme": "basic",
       "token": null,
       "username": "operator",
       "password": "***REDACTED***",
       "api_key": null,
       "api_key_header": "X-API-Key"
     }
   }

``basic_summary`` keeps ``username`` in the clear — a username identifies
who authenticated, not a secret in itself — while ``password`` is
redacted next to it in the same record. The custom
``api_key_header="X-Field-Key"`` set on the credential is what
``api_key.headers()`` actually keys its dictionary by, which is why the
redacted output shows ``X-Field-Key`` rather than the class default of
``X-API-Key`` still visible on ``basic_summary`` (which never set a
custom header name).

Load From Environment
---------------------

:meth:`~pycsamt.iot.security.SecurityConfig.from_env` reads
:term:`environment variable`\ s with the ``PYCSAMT_IOT_`` prefix.
Recognised variables include ``TOKEN``, ``API_KEY``, ``USERNAME``,
``PASSWORD``, ``TLS``, ``CA_CERT``, ``CERTFILE``, and ``KEYFILE``.
Credential discovery follows a deterministic precedence rule:

.. math::

   \mathrm{scheme} =
   \begin{cases}
   \mathrm{bearer}, & \text{if TOKEN is set},\\
   \mathrm{api\_key}, & \text{else if API\_KEY is set},\\
   \mathrm{basic}, & \text{else if USERNAME and PASSWORD are set},\\
   \mathrm{none}, & \text{otherwise.}
   \end{cases}

The TLS flag is parsed from ``PYCSAMT_IOT_TLS`` when present;
certificate paths are copied into the TLS configuration and remain
visible in the redacted summary.

.. code-block:: pycon

   >>> import os
   >>> os.environ["PYCSAMT_IOT_TOKEN"] = "env-token-456"
   >>> os.environ["PYCSAMT_IOT_TLS"] = "true"
   >>> os.environ["PYCSAMT_IOT_CA_CERT"] = "certs/env-ca.pem"
   >>> from_env = SecurityConfig.from_env()
   >>> for key in ["PYCSAMT_IOT_TOKEN", "PYCSAMT_IOT_TLS", "PYCSAMT_IOT_CA_CERT"]:
   ...     _ = os.environ.pop(key, None)
   >>> print(json.dumps(from_env.as_dict(), indent=2))
   {
     "tls": {
       "enabled": true,
       "ca_cert": "certs/env-ca.pem",
       "certfile": null,
       "keyfile": null,
       "verify": true,
       "min_version": null
     },
     "credential": {
       "scheme": "bearer",
       "token": "***REDACTED***",
       "username": null,
       "password": null,
       "api_key": null,
       "api_key_header": "X-API-Key"
     },
     "require_tls": false,
     "allowed_protocols": null
   }

``TOKEN`` being set is enough on its own to select the bearer scheme per
the precedence rule above, even though ``PYCSAMT_IOT_USERNAME`` and
``PYCSAMT_IOT_PASSWORD`` were never touched in this example — bearer
always wins first when present. Note also that ``require_tls`` comes
back ``False``: ``from_env`` reads a TLS *enabled* flag from the
environment, but it never sets ``require_tls`` itself, so a deployment
that wants to enforce TLS still has to set that field explicitly after
loading from the environment.

Protocol Policy
---------------

Use ``allowed_protocols`` as a small :term:`protocol policy` gate before
building telemetry clients. This is especially useful when a deployment
should only use encrypted or centrally approved transports. The check is
case-insensitive after normalisation: if :math:`P` is missing, every
protocol is allowed; otherwise a requested protocol :math:`p` is allowed
only when :math:`\mathrm{lower}(p) \in P`.

.. code-block:: pycon

   >>> print(f"allows mqtt: {security.allows('mqtt')}")
   allows mqtt: True
   >>> print(f"allows serial: {security.allows('serial')}")
   allows serial: False

``security`` still carries the ``allowed_protocols=["https", "mqtt",
"websocket"]`` list set all the way back in the first example on this
page — ``mqtt`` is a literal member of that list, while ``serial`` was
never one of the three protocols this deployment approved, so ``allows``
gates it out before a client is ever built.

Validation Checks
-----------------

Invalid security configurations fail early. This prevents field scripts
from silently running without the required credential or TLS material.
The validation predicates are intentionally simple: bearer requires
:math:`b \ne \varnothing`, basic requires both :math:`u \ne
\varnothing` and :math:`p \ne \varnothing`, API-key authentication
requires :math:`k \ne \varnothing`, client certificate paths must be
supplied as a pair, and ``require_tls=True`` is valid only when TLS is
enabled.

.. code-block:: pycon

   >>> for label, factory in [
   ...     ("missing bearer token", lambda: Credential(scheme="bearer")),
   ...     (
   ...         "require_tls without TLS",
   ...         lambda: SecurityConfig(tls=TLSConfig(enabled=False), require_tls=True),
   ...     ),
   ...     (
   ...         "certfile without keyfile",
   ...         lambda: TLSConfig(enabled=True, certfile="client.pem"),
   ...     ),
   ... ]:
   ...     try:
   ...         factory()
   ...     except Exception as exc:
   ...         print(f"{label}: {type(exc).__name__}: {exc}")
   missing bearer token: ValueError: bearer scheme requires a token.
   require_tls without TLS: ValueError: require_tls is set but TLS is not enabled.
   certfile without keyfile: ValueError: certfile and keyfile must be provided together.

All three failures raise inside ``__post_init__`` — the moment the
object is constructed, not the moment it is first used — so a
misconfigured script fails at the exact line that built the bad
``Credential``/``TLSConfig``/``SecurityConfig``, rather than surfacing a
confusing error deep inside a transport call much later.

Field Guidance
--------------

Use ``as_dict()`` for anything that may be logged, stored in a
manifest, or shown in a notebook. Use ``client_options()`` only at the
boundary where a live transport client is created. Keep certificate
files and environment variables outside the documentation tree and out
of version control.

For production deployments, prefer TLS-enabled transports and
:term:`short-lived credential`\ s. Record the redacted security summary
in the acquisition manifest so reviewers can see which security policy
was active without exposing the actual secrets. This gives later
reviewers the same policy state :math:`s` needed to reproduce client
construction, while withholding the secret components that should never
be committed.
