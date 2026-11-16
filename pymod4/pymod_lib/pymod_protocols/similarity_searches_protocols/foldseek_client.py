# Copyright 2026 by Serena Rosignoli. All rights reserved.
# This code is part of the PyMod package and governed by its license. Please
# see the LICENSE file that should have been included as part of this package
# or the main __init__.py file in the pymod3 folder.

"""
Talking to a Foldseek search server.

Foldseek finds structures which look like a structure, where BLAST and the other
searches of PyMod find sequences which look like a sequence. The server is the
one behind search.foldseek.com, which is MMseqs2-App, so the same client works
against a server a laboratory runs itself.

The conversation is three calls and no key at all:

    POST /api/ticket              the structure, the mode, the databases
         -> {"id": ..., "status": "PENDING" | "COMPLETE" | ...}
    GET  /api/ticket/<id>         until the status is COMPLETE or ERROR
    GET  /api/result/<id>/0       the hits, as JSON

and one more which says what can be searched:

    GET  /api/databases           name, the value to send, whether it has taxonomy

The public server is free and shared. PyMod therefore asks for one search at a
time, waits between polls rather than spinning, and sends only the databases
which were asked for.

**One chain per search.** A file with two chains is accepted and only its first
chain is searched -- measured against the live server: the same 56 hits came back
for a two chain file as for its chain A alone, and the query was reported as
'job_A'. The protocol therefore sends one chain and says which.
"""

import json
import time

try:
    import requests
    has_requests = True
except ImportError:
    has_requests = False


# The public server. A laboratory running MMseqs2-App of its own can be reached
# by changing this, which is why every function takes it as an argument.
DEFAULT_BASE_URL = "https://search.foldseek.com"

# What the server calls the state of a job.
STATUS_PENDING = "PENDING"
STATUS_RUNNING = "RUNNING"
STATUS_COMPLETE = "COMPLETE"
STATUS_ERROR = "ERROR"
STATUS_UNKNOWN = "UNKNOWN"
FINISHED = (STATUS_COMPLETE, STATUS_ERROR)

# How an alignment is computed. '3diaa' is the default of the server: the local
# alignment over the structural alphabet and the amino acids together.
MODE_3DI_AA = "3diaa"
MODE_TMALIGN = "tmalign"
MODE_3DI = "3di"
MODES = ((MODE_3DI_AA, "3Di and amino acids (default)"),
         (MODE_TMALIGN, "TM-align (global, sorted by TM-score)"),
         (MODE_3DI, "3Di only"))

# The databases to tick when nothing has been chosen: the PDB, because a hit
# there is a structure which exists, and the AlphaFold models of Swiss-Prot,
# because a hit there is a protein with a name.
DEFAULT_DATABASES = ("pdb100", "afdb-swissprot")


class FoldseekError(Exception):
    """
    Raised when a search can not be run or its answer can not be read. The
    message is meant to be shown to the user as it is.
    """
    pass


def _require_requests():
    if not has_requests:
        raise FoldseekError("The 'requests' module is needed to reach the Foldseek"
                            " server, and it is not installed in this PyMOL.")


def available_databases(base_url=DEFAULT_BASE_URL, timeout=20,
                        searchable_only=True):
    """
    What this server can be searched against, as [{'name', 'value', 'taxonomy'}].

    Asked of the server rather than written down here: the list has grown from
    nine to eighteen while this was being written.

    Not everything the server lists can answer a structure search. The list also
    carries the indexes of two other tools, and the server marks them:

      * 'motif' -- a Folddisco index, for searching a structural motif rather
        than a whole chain. Sending a chain to one of these ends the search as
        ERROR, and one of them among several good databases ends the whole
        search as ERROR, not only its own part (both measured against the live
        server with a 264 residue chain of 3UP7);
      * 'rna' -- an RNA database, which refuses a protein query at once with
        "selected databases are not valid".

    They are dropped unless 'searchable_only' is set to False. It matters more
    than it looks: the motif indexes are named exactly like the databases they
    are built from ("PDB", "AlphaFold/UniProt50"), so a list which keeps them
    offers the same name twice and the wrong one of the two fails.
    """

    _require_requests()
    try:
        answer = requests.get("%s/api/databases" % base_url.rstrip("/"),
                              timeout=timeout)
    except Exception as e:
        raise FoldseekError("The Foldseek server could not be reached: %s" % e)

    if not answer.ok:
        raise FoldseekError("The Foldseek server answered %s when asked which"
                            " databases it has." % answer.status_code)

    try:
        entries = answer.json().get("databases") or []
    except ValueError:
        raise FoldseekError("The list of databases could not be read.")

    databases = []
    for entry in entries:
        value = (entry or {}).get("path")
        if not value:
            continue
        database = {"name": entry.get("name") or value,
                    "value": value,
                    "taxonomy": bool(entry.get("taxonomy")),
                    "version": entry.get("version") or "",
                    "status": entry.get("status") or "",
                    "motif": bool(entry.get("motif")),
                    "rna": bool(entry.get("rna")),
                    "interface": bool(entry.get("interface"))}
        if searchable_only and (database["motif"] or database["rna"]):
            continue
        databases.append(database)
    if not databases:
        raise FoldseekError("The Foldseek server offered no database to search.")
    return databases


def submit(structure_path, databases, mode=MODE_3DI_AA,
           base_url=DEFAULT_BASE_URL, timeout=120):
    """
    Sends one structure and returns the identifier of the search.

    'structure_path' is a file PyMOL wrote: both the PDB and the mmCIF it writes
    are accepted, which was checked against the live server.
    """

    _require_requests()

    if not databases:
        raise FoldseekError("Choose at least one database to search.")

    data = [("mode", mode)] + [("database[]", name) for name in databases]
    try:
        with open(structure_path, "rb") as handle:
            answer = requests.post("%s/api/ticket" % base_url.rstrip("/"),
                                   files={"q": handle}, data=data, timeout=timeout)
    except OSError as e:
        raise FoldseekError("The structure could not be read: %s" % e)
    except Exception as e:
        raise FoldseekError("The structure could not be sent to the Foldseek"
                            " server: %s" % e)

    if answer.status_code == 429:
        raise FoldseekError("The Foldseek server is limiting requests: it asks"
                            " for fewer searches. Please try again in a while.")
    if not answer.ok:
        raise FoldseekError("The Foldseek server refused the search (%s):\n\n%s"
                            % (answer.status_code, (answer.text or "")[:300]))

    try:
        ticket = answer.json()
    except ValueError:
        raise FoldseekError("The Foldseek server answered something which is not"
                            " a search ticket.")

    if ticket.get("status") == STATUS_ERROR or not ticket.get("id"):
        raise FoldseekError("The Foldseek server did not accept the structure:"
                            "\n\n%s" % json.dumps(ticket)[:300])
    return ticket["id"]


def status(ticket, base_url=DEFAULT_BASE_URL, timeout=30):
    """The state of a search: PENDING, RUNNING, COMPLETE, ERROR or UNKNOWN."""

    _require_requests()
    try:
        answer = requests.get("%s/api/ticket/%s" % (base_url.rstrip("/"), ticket),
                              timeout=timeout)
    except Exception as e:
        raise FoldseekError("The Foldseek server could not be reached while"
                            " waiting for the search: %s" % e)
    if not answer.ok:
        raise FoldseekError("The Foldseek server answered %s while the search was"
                            " waited for." % answer.status_code)
    try:
        return answer.json().get("status") or STATUS_UNKNOWN
    except ValueError:
        return STATUS_UNKNOWN


def wait(ticket, base_url=DEFAULT_BASE_URL, poll_interval=5, max_wait=1800,
         report=None, should_stop=None):
    """
    Waits for a search to finish and returns its final status.

    'report' is called with the status each time it is asked for, so a window can
    say what is happening; 'should_stop' is asked whether to give up, so a dialog
    can be cancelled. The waiting is deliberately unhurried: this is a free
    service shared by everyone who uses it.
    """

    waited = 0
    while True:
        state = status(ticket, base_url=base_url)
        if report is not None:
            report(state)
        if state in FINISHED:
            return state
        if should_stop is not None and should_stop():
            return state
        if waited >= max_wait:
            raise FoldseekError("The search is still running after %s minutes."
                                " Its results can be fetched later from the"
                                " Foldseek website with the ticket %s."
                                % (int(max_wait / 60), ticket))
        time.sleep(poll_interval)
        waited += poll_interval


def results(ticket, entry=0, base_url=DEFAULT_BASE_URL, timeout=180,
            report=None, max_wait=900):
    """
    The hits of a finished search, as the server returns them.

    The answer is read in pieces rather than in one call. A search of many
    databases answers with tens of megabytes, and a server under load sends them
    slowly: 'timeout' is the time allowed between two pieces, not for the whole
    of it, so a trickle which never stops would otherwise be waited on for ever.
    That is what it did -- the search said COMPLETE in the shell and the window
    which says "please wait" stayed open with nothing else happening.

    'report' is called with the megabytes so far, so a window can say what is
    being waited for; 'max_wait' is the longest the whole download may take.
    """

    _require_requests()

    url = "%s/api/result/%s/%s" % (base_url.rstrip("/"), ticket, entry)
    started = time.time()
    try:
        answer = requests.get(url, timeout=(20, timeout), stream=True)
    except Exception as e:
        raise FoldseekError("The results could not be fetched: %s" % e)

    if answer.status_code == 404:
        raise FoldseekError("The Foldseek server has no results for this search"
                            " any more: they are kept only for a while.")
    if not answer.ok:
        raise FoldseekError("The Foldseek server answered %s when asked for the"
                            " results." % answer.status_code)

    pieces, so_far, said = [], 0, 0
    try:
        for piece in answer.iter_content(chunk_size=1024 * 256):
            if not piece:
                continue
            pieces.append(piece)
            so_far += len(piece)
            if report is not None and so_far - said > 2 * 1024 * 1024:
                said = so_far
                report(so_far / (1024.0 * 1024.0))
            if time.time() - started > max_wait:
                raise FoldseekError(
                    "The results of the search are still arriving after %s"
                    " minutes (%.1f MB so far). The Foldseek server is slow"
                    " today.\n\nThe search itself finished: its results can be"
                    " fetched later from %s with the ticket %s, or the search"
                    " can be run again against fewer databases."
                    % (int(max_wait / 60), so_far / (1024.0 * 1024.0),
                       base_url, ticket))
    except FoldseekError:
        raise
    except Exception as e:
        raise FoldseekError(
            "The results stopped arriving after %.1f MB: %s\n\nThe search"
            " itself finished, and its results can be fetched later from %s"
            " with the ticket %s."
            % (so_far / (1024.0 * 1024.0), e, base_url, ticket))

    if report is not None:
        report(so_far / (1024.0 * 1024.0))

    try:
        return json.loads(b"".join(pieces).decode("utf-8"))
    except ValueError:
        raise FoldseekError("The results of the search could not be read.")
