"""ani.*.png（5 秒ごと 121 コマ）から GIF を作る。
  10 秒ごと（1 コマおき、61 コマ）にし、全コマ共通の 255 色の表＋ディザで色の段を目立たなくする（約 5 MB）。
  再現：python3 ani_heater_front/make_gif.py upper（時刻が左上）／lower（左下）
"""
from pathlib import Path
from PIL import Image
import sys
here = Path(__file__).resolve().parent
pos = sys.argv[1] if len(sys.argv) > 1 else 'upper'      # upper（時刻が左上）／lower（左下）
files = sorted((here / f'frames_{pos}').glob('ani.*.png'))
frames = [Image.open(f).convert('RGB').resize((900, 460), Image.LANCZOS) for f in files]
pick = frames[::10]
big = Image.new('RGB', (900, 460 * len(pick)))
for i, f in enumerate(pick):
    big.paste(f, (0, 460 * i))
pal = big.quantize(colors=255, method=Image.Quantize.MEDIANCUT)
gif = [f.quantize(palette=pal, dither=Image.Dither.FLOYDSTEINBERG) for f in frames[::2]]
out = here / f'ani_heater_front_time_{pos}.gif'
gif[0].save(out, save_all=True, append_images=gif[1:], duration=200, loop=0, optimize=True)
print(out, len(gif), out.stat().st_size // 1024, 'KB')
