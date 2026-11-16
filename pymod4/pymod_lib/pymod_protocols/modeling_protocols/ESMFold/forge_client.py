# Copyright 2026 by Serena Rosignoli. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Minimal REST client for the EvolutionaryScale / Forge ESM API.

Adapted from the standalone ``biohub_client.py`` prototype. It deliberately does
NOT use the ESM SDK: it only relies on ``requests`` (already a PyMod
dependency), so integrating it adds nothing new to install.

The client is intentionally thin: it POSTs a JSON payload to
``<base_url>/api/v1/<endpoint>`` with a Bearer token and returns the decoded
JSON (or raw bytes for non-JSON responses). Both the API key and the base URL
are supplied by the caller, which in PyMod reads them from the preferences file
(see the ESM tool parameters), so nothing here is hard-coded.
"""

from __future__ import annotations

import json

import requests


# Default service endpoint. Editable from the PyMod preferences so users can
# point the client at whichever Forge-compatible host their key belongs to.
# Note: the older host 'forge.evolutionaryscale.ai' is deprecated (scheduled
# shut-down 2026-11-28) and the service itself directs clients to 'biohub.ai'.
DEFAULT_FORGE_BASE_URL = "https://biohub.ai"


class ForgeAPIError(RuntimeError):
    """Raised when the Forge REST API returns a non-OK response."""
    pass


# Cap the length of any server-derived error text, so an error dialog can never
# be blown up by a response that echoes the (potentially huge) request payload.
_MAX_ERROR_CHARS = 400


def _concise_error(resp):
    """
    Extract a short, human-readable error message from a failed response.

    FastAPI validation errors (422) echo the full request 'input' back, which
    for a structure request is the entire coordinate array. We pull out only the
    human-readable parts and never include that echoed input.
    """
    try:
        data = resp.json()
    except Exception:
        return (resp.text or "")[:_MAX_ERROR_CHARS]

    if not isinstance(data, dict):
        return (resp.text or "")[:_MAX_ERROR_CHARS]

    message = data.get("message")
    if isinstance(message, str):
        return message[:_MAX_ERROR_CHARS]

    # FastAPI validation errors carry a 'detail' list (either at the top level
    # or nested under 'message'); each item has 'loc'/'msg' plus an echoed
    # 'input' that we deliberately ignore (it can be the whole coordinate array).
    detail = None
    if isinstance(message, dict):
        detail = message.get("detail")
    if detail is None:
        detail = data.get("detail")

    if isinstance(detail, list):
        parts = []
        for item in detail:
            if isinstance(item, dict):
                loc = item.get("loc")
                msg = item.get("msg", "")
                where = ".".join(str(p) for p in loc) if isinstance(loc, list) else ""
                parts.append(("%s: %s" % (where, msg)).strip(": "))
        if parts:
            return "; ".join(parts)[:_MAX_ERROR_CHARS]
    if isinstance(detail, str):
        return detail[:_MAX_ERROR_CHARS]

    return (resp.text or "")[:_MAX_ERROR_CHARS]


class ForgeClient:
    """
    Thin REST client for the Forge ESM API.

    Parameters
    ----------
    api_key:
        Bearer token for the service. Required.
    base_url:
        Service root, without the ``/api/v1`` suffix.
    timeout:
        Per-request timeout, in seconds. Folding/design can be slow, so the
        default is generous.
    """

    def __init__(self, api_key, base_url=DEFAULT_FORGE_BASE_URL, timeout=300):
        if not api_key:
            raise ValueError(
                "Missing ESM/Forge API key. Set it in the PyMod preferences "
                "(Options) before using this feature."
            )
        self.api_key = api_key
        self.base_url = (base_url or DEFAULT_FORGE_BASE_URL).rstrip("/")
        self.timeout = timeout

    def post_json(self, endpoint, payload, return_bytes=False, accept=None):
        """
        POST a JSON payload to ``/api/v1/<endpoint>`` and return the response.

        Parameters
        ----------
        endpoint:
            For example ``"fold"`` or ``"inverse_fold"``.
        payload:
            JSON-serializable request body.
        return_bytes:
            If True, asks the server for a binary payload (``return-bytes``
            header). Most calls keep this False.
        accept:
            Optional ``Accept`` header value.

        Returns
        -------
        The decoded JSON object, or raw ``bytes`` if the response is not JSON.
        """
        url = "%s/api/v1/%s" % (self.base_url, endpoint.lstrip("/"))
        headers = {
            "Authorization": "Bearer %s" % self.api_key,
            "Content-Type": "application/json",
        }
        if return_bytes:
            headers["return-bytes"] = "true"
        if accept:
            headers["Accept"] = accept

        resp = requests.post(url, headers=headers, json=payload, timeout=self.timeout)

        if not resp.ok:
            raise ForgeAPIError(
                "Forge API request failed (%s %s): %s"
                % (resp.status_code, resp.reason, _concise_error(resp))
            )

        ctype = resp.headers.get("content-type", "")
        if "application/json" in ctype:
            return resp.json()

        # Fall back to JSON parsing, then to raw bytes.
        try:
            return resp.json()
        except Exception:
            return resp.content
