# Copyright 2026 by Serena Rosignoli. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Minimal REST client for the Boltz API (https://api.boltz.bio).

Like the Forge client of the ESMFold protocols, this deliberately does NOT use
the vendor SDK: the Boltz API exposes a plain REST interface, so 'requests'
(already a PyMod dependency) is enough and nothing new has to be installed.

The API is asynchronous. A job is started with a POST, which answers with an
identifier, and the job is then polled until it reaches a terminal state. The
results are not returned in the response: they are presigned URLs, valid for a
limited time, which have to be downloaded separately.

    client = BoltzClient(api_key)
    job = client.start(STRUCTURE_RESOURCE, payload)
    job = client.wait(STRUCTURE_RESOURCE, job["id"])
    data = client.download(url_taken_from(job))

Authentication is a single 'x-api-key' header: the API uses neither a Bearer
token nor request signing.

This module knows nothing about the shape of the request payloads. The caller
builds them, which keeps the transport insulated from the parts of the request
schema that may still change.
"""

import json
import time

import requests


# Service root. Editable from the PyMod preferences, so a user can point the
# client at another compatible host.
DEFAULT_BOLTZ_BASE_URL = "https://api.boltz.bio"

# All compute resources live under this prefix.
API_PREFIX = "compute/v1"

# Resources used by PyMod. Passed to the methods below, so that one client can
# serve both the structure prediction and the protein design protocols.
STRUCTURE_RESOURCE = "predictions/structure-and-binding"
PROTEIN_DESIGN_RESOURCE = "protein/design"
PROTEIN_REDESIGN_RESOURCE = "protein/sequence-redesign"

# Model served by the structure and binding resource.
DEFAULT_STRUCTURE_MODEL = "boltz-2.1"

# Job states. 'pending' and 'running' are transient, the other two terminal.
STATUS_PENDING = "pending"
STATUS_RUNNING = "running"
STATUS_SUCCEEDED = "succeeded"
STATUS_FAILED = "failed"
TERMINAL_STATUSES = (STATUS_SUCCEEDED, STATUS_FAILED)

# Downloading a result file needs far more patience than calling the API: the files
# come from object storage, they can be large, and a slow link is not an error. An
# API call should fail quickly instead, so the two have separate timeouts.
DOWNLOAD_TIMEOUT = 600

# A transient read timeout is exactly the failure worth retrying, so a download is
# attempted more than once before giving up, waiting a little in between.
DOWNLOAD_ATTEMPTS = 3
DOWNLOAD_RETRY_WAIT = 3

# A key issued for the test environment carries this marker. Jobs run with such
# a key return synthetic results and consume no GPU time, which makes them
# useful to check that the whole path works without being charged for it.
TEST_KEY_MARKER = "_test_"


class BoltzAPIError(RuntimeError):
    """Raised when the Boltz REST API returns a non-OK response."""
    pass


class BoltzJobFailed(BoltzAPIError):
    """Raised when a job reaches the 'failed' state."""
    pass


class BoltzJobCancelled(BoltzAPIError):
    """Raised when the caller asks to stop waiting for a job."""
    pass


class BoltzJobTimeout(BoltzAPIError):
    """Raised when a job does not reach a terminal state in time."""
    pass


# Cap the length of any server-derived error text, so that an error dialog can
# not be blown up by a response echoing the whole request payload back.
_MAX_ERROR_CHARS = 400


def _concise_error(resp):
    """
    Extracts a short, readable error message from a failed response, without
    including any echoed copy of the request (which for a large complex, or for
    a base64 encoded MSA, would be enormous).
    """

    try:
        data = resp.json()
    except Exception:
        return (resp.text or "")[:_MAX_ERROR_CHARS]

    if not isinstance(data, dict):
        return (resp.text or "")[:_MAX_ERROR_CHARS]

    # A single message is the common case.
    for key in ("message", "error", "detail"):
        value = data.get(key)
        if isinstance(value, str) and value:
            return value[:_MAX_ERROR_CHARS]

    # Validation errors carry a list of items, each locating the offending
    # field. The echoed 'input' of each item is deliberately ignored.
    detail = data.get("detail")
    if isinstance(data.get("error"), dict):
        detail = data["error"].get("detail", detail)
    if isinstance(detail, list):
        parts = []
        for item in detail:
            if isinstance(item, dict):
                loc = item.get("loc")
                msg = item.get("msg", "")
                where = ".".join(str(p) for p in loc) if isinstance(loc, list) else ""
                parts.append(("%s: %s" % (where, msg)).strip(": "))
            elif isinstance(item, str):
                parts.append(item)
        if parts:
            return "; ".join(parts)[:_MAX_ERROR_CHARS]

    return (resp.text or "")[:_MAX_ERROR_CHARS]


class BoltzClient:
    """
    Thin REST client for the Boltz API.

    Parameters
    ----------
    api_key:
        Key for the service, sent in the 'x-api-key' header. Required.
    base_url:
        Service root, without the '/compute/v1' suffix.
    timeout:
        Per-request timeout in seconds. This is the timeout of a single HTTP
        call, not of a prediction: a prediction is waited for by polling, see
        the 'wait' method.
    session:
        Optional object exposing 'get' and 'post' like the 'requests' module.
        Used by the tests to run the client without contacting the network.
    """

    def __init__(self, api_key, base_url=DEFAULT_BOLTZ_BASE_URL, timeout=120,
                 session=None):
        if not api_key:
            raise ValueError(
                "Missing Boltz API key. Set it in the PyMod preferences"
                " (Options) before using this feature."
            )
        self.api_key = api_key
        self.base_url = (base_url or DEFAULT_BOLTZ_BASE_URL).rstrip("/")
        self.timeout = timeout
        self._session = session or requests

    ###############################################################################################
    # Helpers.                                                                                    #
    ###############################################################################################

    def uses_test_key(self):
        """
        True when the key belongs to the test environment, whose jobs return
        synthetic results at no cost. Worth reporting in the GUI, so that a user
        does not mistake a synthetic model for a real prediction.
        """
        return TEST_KEY_MARKER in (self.api_key or "")

    def _url(self, resource, *parts):
        segments = [self.base_url, API_PREFIX, resource.strip("/")]
        segments.extend(str(p).strip("/") for p in parts if p is not None)
        return "/".join(segments)

    def _headers(self):
        return {"x-api-key": self.api_key,
                "Content-Type": "application/json",
                "Accept": "application/json"}

    def _call(self, method, url, payload=None):
        """
        Performs a request and returns the decoded JSON, raising BoltzAPIError
        with a readable message on any failure, including connection errors.
        """

        try:
            if method == "GET":
                resp = self._session.get(url, headers=self._headers(),
                                         timeout=self.timeout)
            else:
                resp = self._session.post(url, headers=self._headers(),
                                          json=payload, timeout=self.timeout)
        except requests.exceptions.RequestException as e:
            raise BoltzAPIError("Could not reach the Boltz API at %s: %s"
                                % (self.base_url, e))

        if not resp.ok:
            raise BoltzAPIError("Boltz API request failed (%s %s): %s"
                                % (resp.status_code,
                                   getattr(resp, "reason", ""),
                                   _concise_error(resp)))

        try:
            return resp.json()
        except Exception:
            raise BoltzAPIError("The Boltz API returned a response which is not"
                                " JSON (content type: %s)."
                                % resp.headers.get("content-type", "unknown"))

    ###############################################################################################
    # Jobs.                                                                                       #
    ###############################################################################################

    def start(self, resource, payload):
        """
        Starts a job and returns the decoded response, which carries at least
        the identifier of the job and its initial status.
        """
        return self._call("POST", self._url(resource), payload)

    def retrieve(self, resource, job_id):
        """Returns the current state of a job."""
        return self._call("GET", self._url(resource, job_id))

    def estimate_cost(self, resource, payload):
        """
        Returns the cost the given request would incur, without running it.
        Used to warn the user before a large complex is submitted.
        """
        return self._call("POST", self._url(resource, "estimate-cost"), payload)

    def list_jobs(self, resource):
        """Returns the jobs already submitted with this key."""
        return self._call("GET", self._url(resource))

    def list_results(self, resource, job_id):
        """
        Returns the results of a pipeline run, which are NOT part of the job
        itself: a design run reports its designs here, one entry each.
        """
        return self._call("GET", self._url(resource, job_id, "results"))

    def delete_data(self, resource, job_id):
        """Asks the service to delete the data of a job."""
        return self._call("POST", self._url(resource, job_id, "delete-data"))

    def wait(self, resource, job_id, poll_interval=5, max_wait=3600,
             should_cancel=None, on_status=None, sleep=time.sleep):
        """
        Polls a job until it succeeds, and returns its final state.

        Parameters
        ----------
        poll_interval:
            Seconds between two polls.
        max_wait:
            Give up after this many seconds and raise BoltzJobTimeout. The job
            itself keeps running on the service and can be retrieved later.
        should_cancel:
            Optional callable returning True when the wait should be abandoned,
            so that a Cancel button can interrupt it. Raises BoltzJobCancelled.
        on_status:
            Optional callable receiving the job state at every poll, to report
            progress in the GUI.
        sleep:
            Injected for the tests, so they do not have to wait.

        Raises BoltzJobFailed when the service reports a failure.
        """

        waited = 0
        while True:
            if should_cancel is not None and should_cancel():
                raise BoltzJobCancelled("Waiting for the Boltz prediction was"
                                        " cancelled. The job may still be"
                                        " running on the server (id: %s)." % job_id)

            job = self.retrieve(resource, job_id)
            status = (job or {}).get("status")

            if on_status is not None:
                on_status(job)

            if status == STATUS_SUCCEEDED:
                return job
            if status == STATUS_FAILED:
                raise BoltzJobFailed("The Boltz prediction failed: %s"
                                     % self.job_error_message(job))
            if status not in (STATUS_PENDING, STATUS_RUNNING):
                # An unknown status is not treated as terminal, to stay usable
                # if the service gains new intermediate states, but it is worth
                # reporting once.
                if on_status is None:
                    print("- Unexpected status from the Boltz API: %r." % status)

            if waited >= max_wait:
                raise BoltzJobTimeout(
                    "The Boltz prediction did not complete within %s seconds."
                    " It may still be running on the server (id: %s)."
                    % (max_wait, job_id))

            sleep(poll_interval)
            waited += poll_interval

    @staticmethod
    def job_error_message(job):
        """Best effort extraction of the reason a job failed."""
        if not isinstance(job, dict):
            return "no details available"
        for key in ("error", "error_message", "message", "failure_reason"):
            value = job.get(key)
            if isinstance(value, str) and value:
                return value[:_MAX_ERROR_CHARS]
            if isinstance(value, dict):
                nested = value.get("message")
                if isinstance(nested, str) and nested:
                    return nested[:_MAX_ERROR_CHARS]
        return "no details available"

    ###############################################################################################
    # Results.                                                                                    #
    ###############################################################################################

    def download(self, url, timeout=None, attempts=DOWNLOAD_ATTEMPTS,
                 sleep=time.sleep):
        """
        Downloads a result file and returns its bytes.

        The API hands out presigned URLs, which already carry their own
        credentials in the query string. The API key is deliberately NOT sent
        here: it is not needed, and adding headers to a presigned request can
        invalidate its signature.

        The timeout is DOWNLOAD_TIMEOUT rather than the one of the API calls: a
        file comes from object storage and may simply be slow, while an API call
        should fail quickly. A network failure is retried, since a read timing out
        once is usually transient; a refusal by the server is not retried, because
        it would fail the same way again.
        """

        timeout = timeout or DOWNLOAD_TIMEOUT
        last_error = None

        for attempt in range(1, max(1, attempts) + 1):
            try:
                resp = self._session.get(url, timeout=timeout)
            except requests.exceptions.RequestException as e:
                last_error = e
                if attempt < attempts:
                    print("- A Boltz result file could not be downloaded (%s)."
                          " Trying again (%s of %s)." % (e, attempt + 1, attempts))
                    sleep(DOWNLOAD_RETRY_WAIT)
                    continue
                break

            if resp.ok:
                return resp.content

            # A refusal is not retried: an expired link is the likely cause, and
            # trying again would not change it.
            raise BoltzAPIError(
                "Could not download a Boltz result file (%s %s). The link may"
                " have expired: the results of a job are kept for a week, so"
                " retrieve the job again to obtain new links."
                % (resp.status_code, getattr(resp, "reason", "")))

        raise BoltzAPIError(
            "Could not download a Boltz result file after %s attempts: %s.\n\n"
            "The prediction itself is safe: its results are kept for a week on the"
            " service, so nothing has to be run again." % (attempts, last_error))

    def download_to_file(self, url, filepath, timeout=None):
        """Downloads a result file and writes it, returning its path."""
        content = self.download(url, timeout)
        with open(filepath, "wb") as fh:
            fh.write(content)
        return filepath
