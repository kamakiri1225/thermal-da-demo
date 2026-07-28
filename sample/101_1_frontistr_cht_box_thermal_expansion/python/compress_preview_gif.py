#!/usr/bin/env python3
"""ParaViewで作った大きなGIFからGitHub README用の軽量版を作る。"""

from __future__ import annotations

import argparse
from pathlib import Path

from PIL import Image


HERE = Path(__file__).resolve().parent
DEFAULT_INPUT = HERE.parent / "ani" / "ani.gif"
DEFAULT_OUTPUT = HERE.parent / "data" / "thermal_expansion_temperature_preview.gif"


def main() -> int:
    parser = argparse.ArgumentParser(description="README掲載用にGIFを縮小・減色する")
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--width", type=int, default=560)
    parser.add_argument("--frame-step", type=int, default=4)
    parser.add_argument("--colors", type=int, default=64)
    parser.add_argument("--duration-ms", type=int, default=200)
    args = parser.parse_args()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    frames: list[Image.Image] = []
    with Image.open(args.input) as source:
        height = round(source.height * args.width / source.width)
        for frame_index in range(0, source.n_frames, args.frame_step):
            source.seek(frame_index)
            resized = source.convert("RGB").resize((args.width, height), Image.Resampling.LANCZOS)
            frames.append(
                resized.quantize(colors=args.colors, method=Image.Quantize.MEDIANCUT)
            )

    if not frames:
        raise ValueError(f"GIFにフレームがありません: {args.input}")
    frames[0].save(
        args.output,
        save_all=True,
        append_images=frames[1:],
        duration=args.duration_ms,
        loop=0,
        optimize=True,
        disposal=2,
    )
    print(
        f"saved {args.output}: {len(frames)} frames, "
        f"{args.width}x{height}, {args.output.stat().st_size / 1024 / 1024:.2f} MiB"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
