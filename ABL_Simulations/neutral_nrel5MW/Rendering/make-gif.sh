#!/usr/bin/env bash
set -euo pipefail

# This is an input option for an image sequence. It tells FFmpeg how long each input image should remain on screen and assigns timestamps to the images.
framerate=12
output="vorticity-mag-preview.gif"
prefix="vorticity-mag"

tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

shopt -s nullglob
frames=( "${prefix}".[0-9][0-9][0-9][0-9].png )

if (( ${#frames[@]} < 2 )); then
    echo "Not enough PNG frames found." >&2
    exit 1
fi

# Sort the frame names numerically.
mapfile -t frames < <(
    printf '%s\n' "${frames[@]}" | sort -V
)

# Exclude the newest PNG because the renderer may still be writing it.
frames=( "${frames[@]:0:${#frames[@]}-1}" )

# Create a fixed snapshot of the current sequence using symlinks.
for i in "${!frames[@]}"; do
    printf -v number '%06d' "$i"
    ln -s "$PWD/${frames[$i]}" "$tmp/frame-${number}.png"
done

# Read the image dimensions from the first frame.
IFS=x read -r width height < <(
    ffprobe \
        -v error \
        -select_streams v:0 \
        -show_entries stream=width,height \
        -of csv=s=x:p=0 \
        "$tmp/frame-000000.png"
)

echo "Using ${#frames[@]} frames"
echo "Dimensions: ${width}x${height}"
echo "Frame rate: ${framerate} "
echo "Generating palette..."

# Generate the palette after compositing the transparent PNGs over black.
ffmpeg \
    -hide_banner \
    -loglevel warning \
    -y \
    -framerate "$framerate" \
    -i "$tmp/frame-%06d.png" \
    -filter_complex \
    "color=c=black:s=${width}x${height}:r=${framerate}[bg];
     [bg][0:v]overlay=shortest=1:eof_action=endall:format=auto,
     format=rgb24,
     palettegen=max_colors=256:stats_mode=full:reserve_transparent=0" \
    -frames:v 1 \
    -update 1 \
    "$tmp/palette.png"

echo "Generating GIF..."

# Apply the same black background before using the generated palette.
ffmpeg \
    -hide_banner \
    -loglevel warning \
    -y \
    -framerate "$framerate" \
    -i "$tmp/frame-%06d.png" \
    -i "$tmp/palette.png" \
    -filter_complex \
    "color=c=black:s=${width}x${height}:r=${framerate}[bg];
     [bg][0:v]overlay=shortest=1:eof_action=endall:format=auto,
     format=rgb24[opaque];
     [opaque][1:v]paletteuse=dither=sierra2_4a:diff_mode=rectangle[out]" \
    -map "[out]" \
    -loop 0 \
    "$output"

echo "Created: $output"
echo "Frames: ${#frames[@]}"
echo "Duration: ${#frames[@]} seconds"
