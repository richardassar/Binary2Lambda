#!/usr/bin/env python3
"""Binary2Lambda performance plots: C++ and Rust benchmarks, rendered with Python.

Runs the `--bench` mode of the C++ and Rust binaries over two sweeps —
bit-string length (for several de Bruijn caps) and de Bruijn cap (at fixed
length) — then renders the scaling curves for every operation: table build,
decode (number -> lambda), encode (lambda -> number), table size on disk,
and peak process memory. Each benchmark block runs as its own process, so
the peak-RSS reading (VmHWM) is a true per-block high-water mark; it
includes the language runtime baseline of a few MB.

Outputs: plots/*.png and plots/bench_results.csv.

Build the binaries first:
  g++ -std=c++17 -O2 -o cpp/lambda_bijection_cpp cpp/lambda_bijection.cpp
  rustc -O -o rust/lambda_bijection_rs rust/lambda_bijection.rs
"""

import csv
import subprocess
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
PLOTS = ROOT / "plots"
BINARIES = {
    "cpp": ROOT / "cpp" / "lambda_bijection_cpp",
    "rust": ROOT / "rust" / "lambda_bijection_rs",
}
COLUMNS = ("lang,cap,length,n_max,build_s,dec_mean_us,dec_std_us,dec_min_us,"
           "dec_max_us,enc_mean_us,enc_std_us,enc_min_us,enc_max_us,entries,"
           "table_bits,disk_bytes,peak_rss_kb").split(",")

LENGTHS = [16, 24, 32, 48, 64, 96, 128, 192, 256, 384, 512, 768, 1024]
LENGTH_CAPS = ["inf", "8", "32", "128"]
CAP_SWEEP = ["4", "8", "16", "32", "64", "128", "256", "inf"]
CAP_LENGTH = 192

CAP_COLORS = {"inf": "tab:blue", "128": "tab:orange", "32": "tab:green",
              "8": "tab:red"}
LANG_STYLE = {"cpp": "-", "rust": "--"}


def run_block(lang: str, cap: str, length: int) -> dict:
    out = subprocess.run([str(BINARIES[lang]), "--bench", cap, str(length)],
                         capture_output=True, text=True, check=True,
                         cwd=ROOT).stdout.strip()
    row = dict(zip(COLUMNS, out.split(",")))
    for key in COLUMNS[2:]:
        row[key] = float(row[key])
    return row


def cap_label(cap: str) -> str:
    return "unbounded" if cap == "inf" else f"K = {cap}"


def collect() -> list[dict]:
    rows = []
    for lang in BINARIES:
        for cap in LENGTH_CAPS:
            for length in LENGTHS:
                rows.append(run_block(lang, cap, length))
                print(f"  {lang} cap={cap} L={length}: "
                      f"build {rows[-1]['build_s']:.3f}s")
        for cap in CAP_SWEEP:
            rows.append(run_block(lang, cap, CAP_LENGTH))
    return rows


def select(rows, lang=None, cap=None, length=None):
    picked = [r for r in rows
              if (lang is None or r["lang"] == lang)
              and (cap is None or r["cap"] == cap)
              and (length is None or r["length"] == length)]
    return sorted(picked, key=lambda r: r["length"])


def plot_time_vs_length(rows):
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
    panels = [("build_s", "table build time (s)"),
              ("dec_mean_us", "decode, mean over 64 samples (µs)"),
              ("enc_mean_us", "encode, mean over 64 samples (µs)")]
    for axis, (key, title) in zip(axes, panels):
        for cap in LENGTH_CAPS:
            for lang in BINARIES:
                data = select(rows, lang, cap)
                axis.plot([r["length"] for r in data], [r[key] for r in data],
                          LANG_STYLE[lang], color=CAP_COLORS[cap],
                          label=f"{cap_label(cap)}, {lang}")
        axis.set_xscale("log", base=2)
        axis.set_yscale("log")
        axis.set_xlabel("bit-string length L")
        axis.set_ylabel("seconds" if key == "build_s" else "microseconds")
        axis.set_title(title)
        axis.grid(True, which="both", alpha=0.3)
    axes[1].legend(fontsize=7, ncols=2)
    fig.suptitle("Binary2Lambda time scaling vs bit-string length "
                 "(C++ solid, Rust dashed; one-time build, per-call codec)")
    fig.tight_layout()
    fig.savefig(PLOTS / "time_vs_length.png", dpi=150)
    plt.close(fig)


def plot_space_vs_length(rows):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
    for cap in LENGTH_CAPS:
        data = select(rows, "cpp", cap)
        axes[0].plot([r["length"] for r in data],
                     [r["disk_bytes"] / 1024 for r in data],
                     "-o", color=CAP_COLORS[cap], markersize=3,
                     label=cap_label(cap))
        for lang in BINARIES:
            data = select(rows, lang, cap)
            axes[1].plot([r["length"] for r in data],
                         [r["peak_rss_kb"] / 1024 for r in data],
                         LANG_STYLE[lang], color=CAP_COLORS[cap],
                         label=f"{cap_label(cap)}, {lang}")
    axes[0].set_title("table size on disk (KiB) — identical in all languages")
    axes[0].set_ylabel("KiB")
    axes[1].set_title("peak process memory (MB, incl. runtime baseline)")
    axes[1].set_ylabel("MB")
    for axis in axes:
        axis.set_xscale("log", base=2)
        axis.set_yscale("log")
        axis.set_xlabel("bit-string length L")
        axis.grid(True, which="both", alpha=0.3)
    axes[0].legend(fontsize=8)
    axes[1].legend(fontsize=7, ncols=2)
    fig.suptitle("Binary2Lambda space scaling vs bit-string length")
    fig.tight_layout()
    fig.savefig(PLOTS / "space_vs_length.png", dpi=150)
    plt.close(fig)


def plot_vs_cap(rows):
    caps_numeric = [int(c) for c in CAP_SWEEP if c != "inf"]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
    panels = [("build_s", "table build time (s)"),
              ("dec_mean_us", "decode mean (µs)"),
              ("enc_mean_us", "encode mean (µs)")]
    for axis, (key, title) in zip(axes, panels):
        for lang in BINARIES:
            swept = [select(rows, lang, str(cap), CAP_LENGTH)[0]
                     for cap in caps_numeric]
            axis.plot(caps_numeric, [r[key] for r in swept],
                      LANG_STYLE[lang], marker="o", markersize=3,
                      label=f"capped, {lang}")
            unbounded = select(rows, lang, "inf", CAP_LENGTH)[0][key]
            axis.axhline(unbounded, linestyle=":", alpha=0.7,
                         color="tab:gray" if lang == "cpp" else "tab:brown",
                         label=f"unbounded, {lang}")
        axis.set_xscale("log", base=2)
        axis.set_xlabel(f"de Bruijn cap K   (L = {CAP_LENGTH} bits)")
        axis.set_ylabel("seconds" if key == "build_s" else "microseconds")
        axis.set_title(title)
        axis.grid(True, which="both", alpha=0.3)
    axes[0].legend(fontsize=8)
    fig.suptitle("Binary2Lambda time vs de Bruijn cap at fixed length "
                 "(capped tables trade term-language coverage for speed)")
    fig.tight_layout()
    fig.savefig(PLOTS / "time_vs_cap.png", dpi=150)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
    swept = [select(rows, "cpp", str(cap), CAP_LENGTH)[0]
             for cap in caps_numeric]
    unbounded = select(rows, "cpp", "inf", CAP_LENGTH)[0]
    axes[0].plot(caps_numeric, [r["disk_bytes"] / 1024 for r in swept],
                 "-o", markersize=3)
    axes[0].axhline(unbounded["disk_bytes"] / 1024, linestyle=":",
                    color="tab:gray", label="unbounded")
    axes[0].set_title("table size on disk (KiB)")
    axes[0].set_ylabel("KiB")
    axes[1].plot(caps_numeric, [r["entries"] for r in swept],
                 "-o", markersize=3)
    axes[1].axhline(unbounded["entries"], linestyle=":", color="tab:gray",
                    label="unbounded")
    axes[1].set_title("table entries  (Θ(K·n) capped vs Θ(n²) unbounded)")
    axes[1].set_ylabel("entries")
    for axis in axes:
        axis.set_xscale("log", base=2)
        axis.set_xlabel(f"de Bruijn cap K   (L = {CAP_LENGTH} bits)")
        axis.grid(True, which="both", alpha=0.3)
        axis.legend(fontsize=8)
    fig.suptitle("Binary2Lambda table size vs de Bruijn cap at fixed length")
    fig.tight_layout()
    fig.savefig(PLOTS / "space_vs_cap.png", dpi=150)
    plt.close(fig)


def load_csv() -> list[dict]:
    with open(PLOTS / "bench_results.csv", newline="") as handle:
        rows = list(csv.DictReader(handle))
    for row in rows:
        for key in COLUMNS[2:]:
            row[key] = float(row[key])
    return rows


def main() -> None:
    import sys

    if "--from-csv" in sys.argv:  # re-render plots from existing data only
        rows = load_csv()
    else:
        for lang, binary in BINARIES.items():
            if not binary.exists():
                raise SystemExit(f"{binary} missing - build the {lang} binary "
                                 f"first (see module docstring)")
        PLOTS.mkdir(exist_ok=True)
        rows = collect()
        with open(PLOTS / "bench_results.csv", "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=COLUMNS)
            writer.writeheader()
            writer.writerows(rows)
    plot_time_vs_length(rows)
    plot_space_vs_length(rows)
    plot_vs_cap(rows)
    print(f"wrote 4 plots to {PLOTS}")


if __name__ == "__main__":
    main()
