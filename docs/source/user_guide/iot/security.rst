.. _user_guide_iot_security:

Security
========

The IoT security helpers define transport security and authentication
configuration for telemetry clients. They do not implement cryptography
directly. TLS, certificates, and authentication headers are passed to the
concrete transport layer, while the configuration objects keep logs,
manifests, and notebooks from leaking secrets.

Use this page for security configuration, redaction, environment loading,
and validation. The examples use placeholder certificate paths and tokens;
replace them with deployment-specific secrets outside version control.

Configure TLS And Bearer Authentication
---------------------------------------

Use :class:`pycsamt.iot.TLSConfig` for TLS material and
:class:`pycsamt.iot.Credential` for authentication. A
:class:`pycsamt.iot.SecurityConfig` combines both and can restrict allowed
telemetry protocols.

.. code-block:: python
   :linenos:

   from pycsamt.iot import AuthScheme, Credential, SecurityConfig, TLSConfig

   security = SecurityConfig(
       tls=TLSConfig(
           enabled=True,
           ca_cert="certs/field-ca.pem",
           certfile="certs/node-client.pem",
           keyfile="certs/node-client.key",
           verify=True,
           min_version="TLSv1.2",
       ),
       credential=Credential(
           scheme=AuthScheme.BEARER,
           token="field-secret-token",
       ),
       require_tls=True,
       allowed_protocols=["https", "mqtt", "websocket"],
   )

   print(security.as_dict())

Output:

.. code-block:: text

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

Build Client Options
--------------------

:meth:`pycsamt.iot.SecurityConfig.client_options` returns the live options
that can be passed to telemetry clients. These options include raw secrets
because the transport needs them. Redact them before printing.

.. code-block:: python
   :linenos:

   from pycsamt.iot import redact_secret

   options = security.client_options()
   print(
       {
           "tls": options.get("tls"),
           "ca_cert": options.get("ca_cert"),
           "certfile": options.get("certfile"),
           "keyfile": options.get("keyfile"),
           "headers": {
               key: redact_secret(value)
               for key, value in options.get("headers", {}).items()
           },
           "token": redact_secret(options.get("token")),
       }
   )

Output:

.. code-block:: text

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

Other Credential Schemes
------------------------

Bearer tokens, API keys, and basic authentication are supported. Use
``headers()`` only when feeding a live transport. Use ``as_dict()`` for
logs, manifests, and notebooks.

.. code-block:: python
   :linenos:

   api_key = Credential(
       scheme="api_key",
       api_key="api-key-123",
       api_key_header="X-Field-Key",
   )
   basic = Credential(
       scheme="basic",
       username="operator",
       password="station-password",
   )

   print(
       {
           "api_key_redacted": {
               key: redact_secret(value)
               for key, value in api_key.headers().items()
           },
           "basic_header_prefix": basic.headers()["Authorization"][:16],
           "basic_summary": basic.as_dict(),
       }
   )

Output:

.. code-block:: text

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

Load From Environment
---------------------

:meth:`pycsamt.iot.SecurityConfig.from_env` reads environment variables
with the ``PYCSAMT_IOT_`` prefix. Recognised variables include ``TOKEN``,
``API_KEY``, ``USERNAME``, ``PASSWORD``, ``TLS``, ``CA_CERT``,
``CERTFILE``, and ``KEYFILE``.

.. code-block:: python
   :linenos:

   import os

   os.environ["PYCSAMT_IOT_TOKEN"] = "env-token-456"
   os.environ["PYCSAMT_IOT_TLS"] = "true"
   os.environ["PYCSAMT_IOT_CA_CERT"] = "certs/env-ca.pem"

   from_env = SecurityConfig.from_env()

   for key in ["PYCSAMT_IOT_TOKEN", "PYCSAMT_IOT_TLS", "PYCSAMT_IOT_CA_CERT"]:
       os.environ.pop(key, None)

   print(from_env.as_dict())

Output:

.. code-block:: text

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

Protocol Policy
---------------

Use ``allowed_protocols`` as a small policy gate before building telemetry
clients. This is especially useful when a deployment should only use
encrypted or centrally approved transports.

.. code-block:: python
   :linenos:

   print(f"allows mqtt: {security.allows('mqtt')}")
   print(f"allows serial: {security.allows('serial')}")

Output:

.. code-block:: text

   allows mqtt: True
   allows serial: False

Validation Checks
-----------------

Invalid security configurations fail early. This prevents field scripts
from silently running without the required credential or TLS material.

.. code-block:: python
   :linenos:

   for label, factory in [
       ("missing bearer token", lambda: Credential(scheme="bearer")),
       (
           "require_tls without TLS",
           lambda: SecurityConfig(tls=TLSConfig(enabled=False),
                                  require_tls=True),
       ),
       (
           "certfile without keyfile",
           lambda: TLSConfig(enabled=True, certfile="client.pem"),
       ),
   ]:
       try:
           factory()
       except Exception as exc:
           print(f"{label}: {type(exc).__name__}: {exc}")

Output:

.. code-block:: text

   missing bearer token: ValueError: bearer scheme requires a token.
   require_tls without TLS: ValueError: require_tls is set but TLS is not enabled.
   certfile without keyfile: ValueError: certfile and keyfile must be provided together.

Field Guidance
--------------

Use ``as_dict()`` for anything that may be logged, stored in a manifest, or
shown in a notebook. Use ``client_options()`` only at the boundary where a
live transport client is created. Keep certificate files and environment
variables outside the documentation tree and out of version control.

For production deployments, prefer TLS-enabled transports and short-lived
tokens. Record the redacted security summary in the acquisition manifest
so reviewers can see which security policy was active without exposing the
actual secrets.
