import BasinMark
import Lean

/- Inspect compiled declarations, so `lemma` and generated declarations cannot
silently escape a handwritten theorem inventory. Audit definitions as well. -/
run_cmd do
  let env ← Lean.getEnv
  let mut count : Nat := 0
  for (name, info) in env.constants.toList do
    if (`BasinMark).isPrefixOf name then
      let kind := match info with
        | .thmInfo _ => "theorem"
        | .defnInfo _ => "definition"
        | .axiomInfo _ => "axiom"
        | _ => "other"
      let axioms ← Lean.collectAxioms name
      let names := String.intercalate "," (axioms.toList.map toString)
      Lean.logInfo m!"BASINMARK_AUDIT {name} {kind} [{names}]"
      count := count + 1
  if count == 0 then
    throwError "No BasinMark declarations found"
  Lean.logInfo m!"BASINMARK_AUDIT_COUNT {count}"
