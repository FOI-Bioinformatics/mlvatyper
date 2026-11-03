"""Logging setup shared by the command line interface and library code.

Besides the human-readable console output, ``--log-json FILE`` attaches a
second handler that writes one JSON object per line (JSON Lines), so a
pipeline monitor can read the events without parsing the text format.
"""

import datetime
import json
import logging

LOGGER_NAME = "mlvatyper"


def get_logger() -> logging.Logger:
    return logging.getLogger(LOGGER_NAME)


class JsonLineFormatter(logging.Formatter):
    """One JSON object per record: time (UTC, ISO 8601), level, logger and message."""

    def format(self, record: logging.LogRecord) -> str:
        event = {
            "time": datetime.datetime.fromtimestamp(record.created, datetime.UTC).isoformat(timespec="milliseconds"),
            "level": record.levelname,
            "logger": record.name,
            "message": record.getMessage(),
        }
        if record.exc_info:
            event["exception"] = self.formatException(record.exc_info)
        return json.dumps(event)


def configure_logging(verbose: bool = False, json_path: str | None = None) -> None:
    """Configure console logging, and JSON Lines to json_path when given.

    The JSON handler is attached to the root logger, so records from merPCR
    and other libraries are captured as well. Calling this again replaces a
    JSON handler added by an earlier call, so repeated calls in one process
    (for example across tests) do not stack handlers or leak file handles.
    """
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%H:%M:%S",
    )
    root = logging.getLogger()
    root.setLevel(level)
    for handler in [h for h in root.handlers if getattr(h, "_mlvatyper_json", False)]:
        root.removeHandler(handler)
        handler.close()
    if json_path:
        handler = logging.FileHandler(json_path, mode="w", encoding="utf-8")
        handler.setFormatter(JsonLineFormatter())
        handler.setLevel(level)
        handler._mlvatyper_json = True  # type: ignore[attr-defined]  # tag so a later call replaces this one
        root.addHandler(handler)
