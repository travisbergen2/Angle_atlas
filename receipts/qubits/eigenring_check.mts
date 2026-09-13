import { runComparison } from "./src/lib/sim/experiment.ts";
import { initSwarm, buildLinearU, DEFAULT_PARAMS } from "./src/lib/sim/swarm.ts";
import { eigenvalues } from "./src/lib/sim/spectrum.ts";
import { writeFileSync } from "node:fs";
for (const evolution of ["additive", "multiplicative"] as const) {
  const series = runComparison({ steps: 160, evolution });
  for (const s of series) {
    const pts = s.points;
    const pick = (t: number) => pts.find((p) => p.t >= t) ?? pts[pts.length - 1];
    console.log(evolution.padEnd(15), s.case.label.padEnd(18), "C0", pts[0].coherence.toFixed(3), "C20", pick(20).coherence.toFixed(3), "C80", pick(80).coherence.toFixed(3), "C160", pick(160).coherence.toFixed(3), "loss160", pick(160).loss.toFixed(3));
  }
}
// spectrum cross-check: dump the additive linear operator for the split, inert and ramified cases
const out: Record<string, unknown> = {};
for (const [label, prime, env] of [["split13", 13, "gaussian"], ["inert11", 11, "gaussian"], ["ram2", 2, "gaussian"]] as const) {
  const st = initSwarm({ ...DEFAULT_PARAMS, prime, env, evolution: "additive" });
  const U = buildLinearU(st);
  const ev = eigenvalues(U);
  out[label] = { n: U.n, re: Array.from(U.re), im: Array.from(U.im), qr_mags: ev.slice(0, 5).map((e) => +e.mag.toFixed(6)), qr_min: +ev[ev.length - 1].mag.toFixed(6) };
  // multiplicative H sign check: eigen-moduli of exp(eta*H)
  const stm = initSwarm({ ...DEFAULT_PARAMS, prime, env, evolution: "multiplicative" });
  const Um = buildLinearU(stm);
  (out[label] as any).mult = { re: Array.from(Um.re), im: Array.from(Um.im) };
}
writeFileSync("/tmp/eigenring_ops.json", JSON.stringify(out));
