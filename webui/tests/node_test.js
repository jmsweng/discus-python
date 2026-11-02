const DB = require("../src/core.js"), EX = require("../src/examples.js"), fs = require("fs");
const tpl = fs.readFileSync(__dirname + "/../engine.py", "utf8");
for (const [k, e] of Object.entries(EX)) {
  const s = e.make();
  const A = DB.latticeMatrix(s.cell);
  const ot = DB.orientationTable(s, A);
  console.log(k, Object.values(ot).map(t => `${t.group.name}: ${t.orientations.length} distinct of ${t.nOps} ops, errors ${t.errors}`));
  const v = DB.validate(s); console.log("  errs", v.errs, "warns", v.warns.length);
  const t0 = Date.now(); const r = DB.realize(s, { n: [6, 6, 6] }); console.log("  atoms", r.atoms.length, "ms", Date.now() - t0, r.report.map(x => [x.before.toFixed(2), x.after.toFixed(2)]));
  const t1 = Date.now(); const pl = DB.diffusePlaneRows(r, s.calc.preview); const row = pl.row(Math.floor(pl.P / 2)); console.log("  plane P", pl.P, "row ms", Date.now() - t1, "max", Math.max(...row).toFixed(1));
  const outDir = process.env.OUT || require("os").tmpdir();
  fs.writeFileSync(`${outDir}/gen_${k}.py`, DB.generateScript(s, tpl));
}
// parser checks
const A = DB.latticeMatrix({ a: 4, b: 4, c: 4, alpha: 90, beta: 90, gamma: 90 });
for (const op of ["-y,x,z", "C3 [1 1 1]", "m [1 -1 0]", "S4 <0 0 1>", "-1", "R 30 [0 0 1]", "x-y,x,z"]) {
  try { const o = DB.parseOp(op, A); console.log(op, "ok improper", o.improper); } catch (e) { console.log(op, "ERR", e.message); }
}
const hexA = DB.latticeMatrix({ a: 4, b: 4, c: 6, alpha: 90, beta: 90, gamma: 120 });
console.log("hex x-y,x,z", DB.parseOp("x-y,x,z", hexA).R.map(r=>r.map(x=>x.toFixed(3))).join(" | "));
console.log("Oh closure", DB.closure(["C4 [0 0 1]","C3 [1 1 1]","-1"].map(o=>DB.parseOp(o,A).R)).length);
