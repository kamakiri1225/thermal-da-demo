"""blog_013の説明確認。熱FEM・ROMを比較し、既存の校正結果は上書きしない。

実行: OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MPLCONFIGDIR=/tmp/mpl-fisher \
      python3 run/check_blog013_explanation.py
"""
from pathlib import Path
import json
import time
import numpy as np
import krylov_rom as kr

ROOT = Path(__file__).resolve().parents[1]

def main():
    fem = kr.build_fem()
    times = np.arange(0.0, 601.0, 30.0)
    cases = ["learned", "intermittent"]
    begin = time.perf_counter()
    full = {case: kr.run_full(fem, kr.QTRUE[case], times) for case in cases}
    out = {
        "note": "Explanation audit: maximum absolute ROM/full-FEM difference over all nodes and 21 output times. No CHT rerun or data assimilation.",
        "n_nodes": len(fem["f"]),
        "integration_dt_s": kr.DT,
        "output_times_s": times.tolist(),
        "h_scale": 1.0,
        "initial_temperature_rise_K": 0.0,
        "full_two_cases_time_s": time.perf_counter()-begin,
        "orders": {},
    }
    for r in [2, 3, 5, 10, 20]:
        V = kr.krylov_basis(fem, r)
        orth_error = float(np.linalg.norm(V.T @ (fem["C"] @ V)-np.eye(r), 2))
        item = {"C_orthogonality_error_2norm": orth_error, "cases": {}}
        assert orth_error < 1e-10, (r, orth_error)
        for case in cases:
            begin = time.perf_counter()
            prediction = kr.run_rom(fem, V, kr.QTRUE[case], times)
            item["cases"][case] = {
                "max_abs_ROM_minus_full_K": float(np.abs(prediction-full[case]).max()),
                "time_rom_s": time.perf_counter()-begin,
            }
        out["orders"][str(r)] = item
    path = ROOT / "results/blog013_explanation_check.json"
    path.write_text(json.dumps(out, ensure_ascii=False, indent=2)+"\n")
    print(path)
    print("C-orthogonality checked for r=2,3,5,10,20.")
    for r,item in out["orders"].items():
        print(r, {case:v["max_abs_ROM_minus_full_K"] for case,v in item["cases"].items()})

if __name__=="__main__":
    main()
