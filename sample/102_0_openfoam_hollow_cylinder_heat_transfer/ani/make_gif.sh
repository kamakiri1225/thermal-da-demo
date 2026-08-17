#!/bin/sh
set -eu

cd "${0%/*}"

command -v convert >/dev/null 2>&1 || {
    echo "ImageMagick convert was not found." >&2
    exit 1
}

convert \
    -delay 10 \
    -loop 0 \
    ani.*.png \
    -resize 900x \
    -colors 96 \
    -layers Optimize \
    ani.gif

identify ani.gif | sed -n '1p'
du -h ani.gif
