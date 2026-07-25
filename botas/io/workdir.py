from __future__ import annotations
from pathlib import Path
from datetime import datetime
import json
import os

class WorkDir:
    def __init__(self, root: Path):
        self.root = root
        self.logs = root / "logs"
        self.tmp = root / "tmp"
        self.results = root / "results"
        self.manifest = root / "run_manifest.json"
        self.readme = root / "README.txt"
        self.logfile = self.logs / "botas.log"


def _auto_name() -> str:
    return f"botas_run_{datetime.now().strftime('%Y%m%d_%H%M%S')}"


def init_workdir(workdir: str | None, out: str | None = None) -> WorkDir:
    if workdir:
        root = Path(workdir)

    elif out:
        p = Path(out)

        if p.name.endswith(".botas.idx"):
            stem = p.name.removesuffix(".botas.idx")
        else:
            stem = p.stem

        root = p.parent / f"{stem}.botas"

    else:
        root = Path.cwd() / _auto_name()

    root.mkdir(parents=True, exist_ok=True)

    wd = WorkDir(root)

    for d in (wd.logs, wd.tmp, wd.results):
        d.mkdir(exist_ok=True)

    if not wd.readme.exists():
        wd.readme.write_text(
            "BOTAS working directory\n"
            "logs/     runtime logs\n"
            "tmp/      temporary files\n"
            "results/  reports and tables\n"
            "run_manifest.json  command + environment\n"
        )

    if not wd.manifest.exists():
        wd.manifest.write_text("{}")

    return wd


def update_manifest(wd: WorkDir, payload: dict) -> None:
    try:
        data = json.loads(wd.manifest.read_text())
    except Exception:
        data = {}

    data.update(payload)
    wd.manifest.write_text(json.dumps(data, indent=2))