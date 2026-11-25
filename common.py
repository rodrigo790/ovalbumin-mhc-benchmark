"""Shared conventions for this standalone analysis."""
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
RAW = ROOT / "raw"
RESULTS = ROOT / "results"


def download(filename, url, *, directory=RAW):
    """Cache a public input and verify its frozen SHA-256 when recorded."""
    import requests

    directory.mkdir(parents=True, exist_ok=True)
    target = directory / filename
    if not target.is_file():
        temporary = target.with_name(target.name + ".part")
        with requests.get(url, stream=True, timeout=(20, 150)) as response:
            response.raise_for_status()
            if response.status_code != 200:
                raise ValueError(f"Unexpected status {response.status_code} for {url}")
            with temporary.open("wb") as stream:
                for chunk in response.iter_content(1024 * 1024):
                    stream.write(chunk)
        check_hash(temporary, filename)
        temporary.replace(target)
    else:
        check_hash(target, filename)
    return target


def check_hash(path, filename):
    manifest = DATA / "sources.json"
    if manifest.is_file():
        expected = json.loads(manifest.read_text(encoding="utf-8")).get(filename, {}).get("sha256")
        if expected:
            with path.open("rb") as stream:
                actual = hashlib.file_digest(stream, "sha256").hexdigest()
            if actual != expected:
                raise ValueError(f"Unexpected SHA-256 for {filename}: {actual}")


def save_summary(values):
    RESULTS.mkdir(exist_ok=True)
    (RESULTS / "summary.json").write_text(
        json.dumps(values, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )


def plot_style():
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 10,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.titleweight": "bold", "figure.dpi": 160,
        "axes.labelcolor": "#1b2c43", "text.color": "#1b2c43",
        "axes.grid": True, "grid.alpha": 0.18,
        "savefig.bbox": "tight",
    })
    return plt
