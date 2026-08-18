from __future__ import annotations

import argparse
from pathlib import Path

from mutflow import __version__
from mutflow.preflight import PreflightError, run_preflight
from mutflow.runner import RunExecutionError, run_workflow
from mutflow.wizard import InitError, run_init


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="mutflow",
        description="Deterministic batch protein mutation workflow toolkit.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    subparsers = parser.add_subparsers(dest="command", required=True)

    init = subparsers.add_parser(
        "init",
        help="Interactively create one deterministic workflow YAML.",
    )
    init.add_argument(
        "destination",
        nargs="?",
        type=Path,
        default=Path("workflow.yaml"),
    )
    preflight = subparsers.add_parser(
        "preflight",
        help="Validate and expand a workflow without running scientific backends.",
    )
    preflight.add_argument("workflow", type=Path)
    run = subparsers.add_parser(
        "run",
        help="Execute an implemented workflow after mandatory preflight.",
    )
    run.add_argument(
        "--resume",
        action="store_true",
        help="Resume only compatible NOT_RUN variants in the configured output directory.",
    )
    run.add_argument("workflow", type=Path)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "init":
        try:
            run_init(args.destination)
        except (InitError, KeyboardInterrupt, EOFError) as exc:
            print("MUTFLOW INIT")
            print("status: REFUSED" if isinstance(exc, InitError) else "state: INTERRUPTED")
            if isinstance(exc, InitError):
                print(f"error: {exc}")
                return 2
            return 130
        return 0

    if args.command == "preflight":
        try:
            report = run_preflight(args.workflow)
        except PreflightError as exc:
            print("MUTFLOW PREFLIGHT")
            print("status: REFUSED")
            print(f"error: {exc}")
            return 2
        print(report.render())
        return {"READY": 0, "CHECK": 1, "REFUSED": 2}[report.status]

    if args.command == "run":
        command = (
            f"mutflow run --resume {args.workflow}"
            if args.resume
            else f"mutflow run {args.workflow}"
        )
        try:
            store = run_workflow(
                args.workflow,
                command=command,
                resume=args.resume,
            )
        except KeyboardInterrupt:
            print("MUTFLOW RUN")
            print("state: INTERRUPTED")
            return 130
        except RunExecutionError as exc:
            print("MUTFLOW RUN")
            print("status: REFUSED")
            print(f"error: {exc}")
            return 2
        print("MUTFLOW RUN")
        print(f"output: {store.output}")
        print(f"state: {store.run_document['state']}")
        print(f"requested: {store.run_document['counts']['requested']}")
        print(f"ok: {store.run_document['counts']['ok']}")
        print(f"failed: {store.run_document['counts']['failed']}")
        return 0 if store.run_document["state"] == "COMPLETED" else 2

    raise AssertionError(f"Unhandled command: {args.command}")
