(* ::Package:: *)

(* ::Title:: *)
(*Binary2Lambda: Lambda Binarization*)


(* ::Text:: *)
(*A bijection between binary strings and closed untyped lambda terms: every binary string (including the empty string) denotes exactly one closed lambda term, and every closed lambda term has exactly one binary string. Encoding and decoding share one incrementally-built counting table, which can be saved to and loaded from disk.*)


(* ::Text:: *)
(*Specification (identical to the companion Python, C++ and Rust implementations, cross-validated against them):*)
(*  terms - LambdaVar[i] with i >= 1 (de Bruijn index), LambdaAbs[body], LambdaApp[fun, arg];*)
(*  size - |Var i| = i+1, |Abs b| = |b|+2, |App f a| = |f|+|a|+2 (bits of Tromp's binary lambda calculus);*)
(*  order - ascending size; within a size class Var < Abs < App, abstractions ordered by body rank, applications by (left size, left rank, right rank);*)
(*  numeration - string s corresponds to the integer N = FromDigits["1" <> s, 2] - 1, so the leading 1 is implicit and the empty string is N = 0.*)


(* ::Text:: *)
(*A finite de Bruijn index cap K selects the sublanguage of closed terms whose indices never exceed K. Each cap value (Infinity included) defines a DIFFERENT bijection, so encode and decode must use the same cap; capped and unbounded bijections agree on all terms of size at most K+1. A finite cap shrinks the counting table from quadratically many to linearly many entries per size, which is what makes very long strings affordable.*)


(* ::Text:: *)
(*Incrementality is native here: counts are memoised on demand, so both a string-length bump and a cap bump only ever compute the missing entries. Memo keys are normalised (effective context Min[m, n-1, cap], effective cap Min[cap, n-1]), so count entries for small sizes are shared between caps automatically (the cumulative-count entries are keyed by raw cap and are not shared).*)


(* ::Section:: *)
(*Interface*)


(* ::Code:: *)
BeginPackage["LambdaBinarization`"];

LambdaVar::usage =
  "LambdaVar[i] is the de Bruijn variable with index i (an integer >= 1).";
LambdaAbs::usage =
  "LambdaAbs[body] is a lambda abstraction with the given body term.";
LambdaApp::usage =
  "LambdaApp[fun, arg] is the application of fun to arg.";
LambdaTermQ::usage =
  "LambdaTermQ[expr] gives True if expr is a well-formed lambda term.";
LambdaTermSize::usage =
  "LambdaTermSize[term] gives the binary-lambda-calculus bit size of term.";
MaxDeBruijnIndex::usage =
  "MaxDeBruijnIndex[term] gives the largest de Bruijn index occurring in \
term (0 if no variable occurs). Useful for choosing a sufficient index cap.";
TermCount::usage =
  "TermCount[n, m] or TermCount[n, m, cap] gives the number of terms of \
size n whose free indices are all <= m and whose indices nowhere exceed \
cap (default Infinity).";
ClosedTermCumulative::usage =
  "ClosedTermCumulative[n] or ClosedTermCumulative[n, cap] gives the \
number of closed terms of size at most n.";
BuildLambdaTable::usage =
  "BuildLambdaTable[sizeLimit] or BuildLambdaTable[sizeLimit, cap] forces \
the counting table up to the given size and returns a summary association. \
Optional: every function builds the entries it needs on demand.";
ClearLambdaTable::usage =
  "ClearLambdaTable[] discards all memoised counting-table entries.";
SaveLambdaTable::usage =
  "SaveLambdaTable[path] writes every memoised counting-table entry to the \
given file as plain Wolfram Language data; returns a summary association.";
LoadLambdaTable::usage =
  "LoadLambdaTable[path] merges the counting-table entries in the given file \
(written by SaveLambdaTable) into the current table, keeping existing \
entries; returns a summary association.";
EncodeLambdaTerm::usage =
  "EncodeLambdaTerm[term] or EncodeLambdaTerm[term, cap] gives the binary \
string (possibly empty) denoting the given closed lambda term.";
DecodeBitString::usage =
  "DecodeBitString[bits] or DecodeBitString[bits, cap] gives the closed \
lambda term denoted by the given binary string. Total: every string of 0s \
and 1s, including \"\", denotes a term.";
LambdaTermForm::usage =
  "LambdaTermForm[term] gives a readable one-line rendering of term, for \
example \"\[Lambda]1 (\[Lambda]\[Lambda]1)\".";
BitStringForm::usage =
  "BitStringForm[bits] gives bits itself, or \[CurlyEpsilon] for the empty string.";
LambdaTermTree::usage =
  "LambdaTermTree[term] gives the expression tree of term as a Tree object: \
\[Lambda] for abstractions, @ for applications, indices at the leaves.";
LambdaBijectionSelfTest::usage =
  "LambdaBijectionSelfTest[] runs consistency checks against constants \
cross-validated with the Python, C++ and Rust implementations; returns an \
association of booleans.";

EncodeLambdaTerm::notterm = "`1` is not a well-formed lambda term.";
EncodeLambdaTerm::notclosed = "`1` is not a closed term (it has a free variable).";
EncodeLambdaTerm::overcap =
  "Term uses de Bruijn indices above the table cap `1`.";
DecodeBitString::badbits = "`1` is not a string of 0s and 1s.";
LoadLambdaTable::badfile = "`1` does not contain a saved lambda table.";


(* ::Section:: *)
(*Terms*)


(* ::Text:: *)
(*Terms are ordinary inert expressions built from LambdaVar, LambdaAbs and LambdaApp. The size measure is the bit length of the term's code in Tromp's binary lambda calculus; it is the yardstick by which the bijection orders terms. MaxDeBruijnIndex supports the capped workflow: measure a term's deepest index, then use any table whose cap is at least that value.*)


(* ::Code:: *)
Begin["`Private`"];

LambdaTermQ[LambdaVar[i_Integer]] := i >= 1;
LambdaTermQ[LambdaAbs[body_]] := LambdaTermQ[body];
LambdaTermQ[LambdaApp[fun_, arg_]] := LambdaTermQ[fun] && LambdaTermQ[arg];
LambdaTermQ[_] := False;

LambdaTermSize[LambdaVar[i_]] := i + 1;
LambdaTermSize[LambdaAbs[body_]] := LambdaTermSize[body] + 2;
LambdaTermSize[LambdaApp[fun_, arg_]] :=
  LambdaTermSize[fun] + LambdaTermSize[arg] + 2;

MaxDeBruijnIndex[LambdaVar[i_]] := i;
MaxDeBruijnIndex[LambdaAbs[body_]] := MaxDeBruijnIndex[body];
MaxDeBruijnIndex[LambdaApp[fun_, arg_]] :=
  Max[MaxDeBruijnIndex[fun], MaxDeBruijnIndex[arg]];

(* True if term is closed: every de Bruijn index is bound by an enclosing
   abstraction (1 <= index <= depth).  EncodeLambdaTerm requires this: the
   bijection ranges over closed terms, so only a closed term has a string. *)
closedTermQ[LambdaVar[i_], depth_] := 1 <= i <= depth;
closedTermQ[LambdaAbs[body_], depth_] := closedTermQ[body, depth + 1];
closedTermQ[LambdaApp[fun_, arg_], depth_] :=
  closedTermQ[fun, depth] && closedTermQ[arg, depth];


(* ::Section:: *)
(*Printing and trees*)


(* ::Text:: *)
(*LambdaTermForm renders a term on one line with the usual conventions: an abstraction body extends as far right as possible, application associates to the left, and parentheses appear only where required. LambdaTermTree gives the expression tree as a Tree object, with \[Lambda] at abstraction nodes, @ at application nodes and de Bruijn indices at the leaves.*)


(* ::Code:: *)
(* context: 0 = top level / abstraction body, 1 = left of an application,
   2 = right of an application *)
showTerm[LambdaVar[i_], _] := IntegerString[i];
showTerm[LambdaAbs[body_], context_] :=
  parenthesizeIf[context > 0, "\[Lambda]" <> showTerm[body, 0]];
showTerm[LambdaApp[fun_, arg_], context_] :=
  parenthesizeIf[context == 2, showTerm[fun, 1] <> " " <> showTerm[arg, 2]];

parenthesizeIf[True, rendered_] := "(" <> rendered <> ")";
parenthesizeIf[False, rendered_] := rendered;

LambdaTermForm[term_] := showTerm[term, 0];

BitStringForm[""] := "\[CurlyEpsilon]";
BitStringForm[bits_String] := bits;

LambdaTermTree[LambdaVar[i_]] := Tree[i, None];
LambdaTermTree[LambdaAbs[body_]] := Tree["\[Lambda]", {LambdaTermTree[body]}];
LambdaTermTree[LambdaApp[fun_, arg_]] :=
  Tree["@", {LambdaTermTree[fun], LambdaTermTree[arg]}];


(* ::Section:: *)
(*The counting table*)


(* ::Text:: *)
(*TermCount[n, m, cap] is the number of terms of size n whose free de Bruijn indices are all at most m and whose indices nowhere exceed cap. It satisfies the recurrence: one variable term (the index n-1, if it is allowed), plus an abstraction for every body of size n-2 counted in the extended context m+1, plus an application for every way of splitting the remaining size between two subterms (a Catalan-style convolution).*)


(* ::Text:: *)
(*Two saturation facts keep the table small and make incremental growth sound. First, a term of size n cannot contain an index above n-1, so both the context and the cap saturate there; memo keys are normalised accordingly, and entries for sizes up to cap+1 are therefore shared between different caps. Second, the recurrence only ever refers to strictly smaller sizes, so growing the table - in either axis - never invalidates an existing entry.*)


(* ::Code:: *)
effectiveContext[n_, m_, cap_] := Min[m, n - 1, cap];

TermCount[n_, m_] := TermCount[n, m, Infinity];
TermCount[n_, m_, cap_] /; n < 2 := 0;
TermCount[n_, m_, cap_] :=
  storedCount[n, effectiveContext[n, m, cap], Min[cap, n - 1]];

ClosedTermCumulative[n_] := ClosedTermCumulative[n, Infinity];
ClosedTermCumulative[n_, cap_] := storedCumulative[Max[n, 3], cap];

ClearLambdaTable[] := (
  Clear[storedCount, storedCumulative];
  storedCount[n_, m_, cap_] := storedCount[n, m, cap] =
    Boole[m == n - 1] +                                  (* the variable n-1 *)
    TermCount[n - 2, m + 1, cap] +                       (* abstractions *)
    Total[TermCount[#, m, cap] TermCount[n - 2 - #, m, cap] & /@
      Range[2, n - 4]];                                  (* applications *)
  storedCumulative[n_, cap_] /; n < 4 := 0;
  storedCumulative[n_, cap_] := storedCumulative[n, cap] =
    storedCumulative[n - 1, cap] + TermCount[n, 0, cap];
);
ClearLambdaTable[];

BuildLambdaTable[sizeLimit_Integer] := BuildLambdaTable[sizeLimit, Infinity];
BuildLambdaTable[sizeLimit_Integer, cap_] := (
  Scan[TermCount[#, 0, cap] &, Range[2, sizeLimit]];
  <|"sizeLimit" -> sizeLimit, "indexCap" -> cap,
    "closedTermsUpToLimit" -> ClosedTermCumulative[sizeLimit, cap]|>);


(* ::Section:: *)
(*Saving and loading the table*)


(* ::Text:: *)
(*The table is the expensive shared artifact, so it can be persisted. SaveLambdaTable harvests every memoised entry (the literal-key down values of the two stores) and writes them with Put as plain Wolfram Language data - portable, human-inspectable, nothing but WL. LoadLambdaTable merges a saved file into the current session: existing entries are kept, loaded ones are added, and subsequent computation continues incrementally from the union.*)


(* ::Code:: *)
(* The harvest patterns are wrapped in HoldPattern because RuleDelayed holds
   only its right-hand side: without the wrapper, storedCount[n_Integer, ...]
   would evaluate while the pattern is being assembled and recurse through
   the general memo rule. *)
storedCountEntries[] :=
  Cases[DownValues[storedCount],
    HoldPattern[RuleDelayed[
      Verbatim[HoldPattern][storedCount[n_Integer, m_Integer, c_Integer]],
      value_Integer]] :> {n, m, c, value}];

storedCumulativeEntries[] :=
  Cases[DownValues[storedCumulative],
    HoldPattern[RuleDelayed[
      Verbatim[HoldPattern][
        storedCumulative[n_Integer, cap : (_Integer | Infinity)]],
      value_Integer]] :> {n, cap, value}];

SaveLambdaTable[path_String] :=
  Module[{counts = storedCountEntries[],
      cumulatives = storedCumulativeEntries[]},
    Put[<|"format" -> "LambdaBinarizationTable",
        "counts" -> counts, "cumulatives" -> cumulatives|>, path];
    <|"path" -> path, "countEntries" -> Length[counts],
      "cumulativeEntries" -> Length[cumulatives]|>];

LoadLambdaTable[path_String] :=
  Module[{data = Get[path]},
    If[!AssociationQ[data] ||
        Lookup[data, "format"] =!= "LambdaBinarizationTable",
      Message[LoadLambdaTable::badfile, path]; Return[$Failed]];
    Scan[Apply[Function[{n, m, c, value}, storedCount[n, m, c] = value]],
      data["counts"]];
    Scan[Apply[Function[{n, cap, value}, storedCumulative[n, cap] = value]],
      data["cumulatives"]];
    <|"path" -> path, "countEntries" -> Length[data["counts"]],
      "cumulativeEntries" -> Length[data["cumulatives"]]|>];


(* ::Section:: *)
(*Encoding and decoding*)


(* ::Text:: *)
(*Both directions traverse the same partition of each size class: first the variable (if one is allowed), then the abstraction block, then the application blocks ordered by left-subterm size. Decoding (number to term) locates the remaining rank inside that partition, splitting it at application nodes with QuotientRemainder; encoding (term to number) reads the very same offsets directly, with no searching, because the term itself says which block it is in. The two functions are exact inverses by construction.*)


(* ::Code:: *)
varCount[n_, m_, cap_] := Boole[n >= 2 && n - 1 <= Min[m, cap]];

(* rank of a term within its size class, in context m *)
termRank[LambdaVar[_], _, _] := 0;
termRank[term : LambdaAbs[body_], m_, cap_] :=
  varCount[LambdaTermSize[term], m, cap] + termRank[body, m + 1, cap];
termRank[term : LambdaApp[fun_, arg_], m_, cap_] :=
  Module[{n = LambdaTermSize[term], leftSize = LambdaTermSize[fun]},
    varCount[n, m, cap] + TermCount[n - 2, m + 1, cap] +
      Total[TermCount[#, m, cap] TermCount[n - 2 - #, m, cap] & /@
        Range[2, leftSize - 1]] +
      termRank[fun, m, cap] TermCount[LambdaTermSize[arg], m, cap] +
      termRank[arg, m, cap]];

(* the rank-th term of size n in context m; total for 0 <= rank < T(n,m) *)
unrankTerm[rank_, n_, m_, cap_] :=
  With[{v = varCount[n, m, cap]},
    If[v == 1 && rank == 0,
      LambdaVar[n - 1],
      With[{rest = rank - v, absBlock = TermCount[n - 2, m + 1, cap]},
        If[rest < absBlock,
          LambdaAbs[unrankTerm[rest, n - 2, m + 1, cap]],
          unrankApplication[rest - absBlock, n, m, cap]]]]];

unrankApplication[rank_, n_, m_, cap_] :=
  Module[{leftSizes = Range[2, n - 4], blocks, partialSums, position, k,
      withinBlock, rightCount, quotientRemainder},
    blocks = TermCount[#, m, cap] TermCount[n - 2 - #, m, cap] & /@ leftSizes;
    partialSums = Accumulate[blocks];
    position = LengthWhile[partialSums, # <= rank &] + 1;
    k = leftSizes[[position]];
    withinBlock = rank - If[position == 1, 0, partialSums[[position - 1]]];
    rightCount = TermCount[n - 2 - k, m, cap];
    quotientRemainder = QuotientRemainder[withinBlock, rightCount];
    LambdaApp[unrankTerm[quotientRemainder[[1]], k, m, cap],
      unrankTerm[quotientRemainder[[2]], n - 2 - k, m, cap]]];

EncodeLambdaTerm[term_] := EncodeLambdaTerm[term, Infinity];
EncodeLambdaTerm[term_, cap_] :=
  Module[{n},
    If[!LambdaTermQ[term],
      Message[EncodeLambdaTerm::notterm, term]; Return[$Failed]];
    If[!closedTermQ[term, 0],
      Message[EncodeLambdaTerm::notclosed, term]; Return[$Failed]];
    If[MaxDeBruijnIndex[term] > cap,
      Message[EncodeLambdaTerm::overcap, cap]; Return[$Failed]];
    n = LambdaTermSize[term];
    StringDrop[
      IntegerString[
        ClosedTermCumulative[n - 1, cap] + termRank[term, 0, cap] + 1, 2],
      1]];

DecodeBitString[bits_String] := DecodeBitString[bits, Infinity];
DecodeBitString[bits_String, cap_] :=
  Module[{number, n},
    If[!StringMatchQ[bits, ("0" | "1") ...],
      Message[DecodeBitString::badbits, bits]; Return[$Failed]];
    number = FromDigits["1" <> bits, 2] - 1;
    n = NestWhile[# + 1 &, 4, ClosedTermCumulative[#, cap] <= number &];
    unrankTerm[number - ClosedTermCumulative[n - 1, cap], n, 0, cap]];


(* ::Section:: *)
(*Self-test*)


(* ::Text:: *)
(*The embedded constants are cross-validated against the Python, C++ and Rust implementations (which produce byte-identical output and are themselves validated against exhaustive enumeration of all bit strings through a binary-lambda-calculus parser). This test checks those constants and additionally exercises round trips in both directions, the capped bijection, agreement of capped and unbounded bijections on small sizes, rejection of non-closed terms by EncodeLambdaTerm, and a save/load cycle through a temporary file. It does not regenerate the 800 cross-language vectors at run time.*)


(* ::Code:: *)
bitsForIndex[number_Integer] := StringDrop[IntegerString[number + 1, 2], 1];

LambdaBijectionSelfTest[] :=
  Module[{expectedCounts, expectedCappedCounts, expectedFirstTerms,
      tempPath, savedValue, results},
    expectedCounts = {1, 0, 1, 1, 2, 1, 6, 5, 13, 14, 37};
    expectedCappedCounts = {1, 0, 1, 1, 2, 1, 5, 5, 12, 13, 30};
    expectedFirstTerms = {
      LambdaAbs[LambdaVar[1]],
      LambdaAbs[LambdaAbs[LambdaVar[1]]],
      LambdaAbs[LambdaAbs[LambdaVar[2]]],
      LambdaAbs[LambdaAbs[LambdaAbs[LambdaVar[1]]]],
      LambdaAbs[LambdaApp[LambdaVar[1], LambdaVar[1]]],
      LambdaAbs[LambdaAbs[LambdaAbs[LambdaVar[2]]]],
      LambdaAbs[LambdaAbs[LambdaAbs[LambdaVar[3]]]],
      LambdaAbs[LambdaAbs[LambdaAbs[LambdaAbs[LambdaVar[1]]]]],
      LambdaAbs[LambdaAbs[LambdaApp[LambdaVar[1], LambdaVar[1]]]],
      LambdaAbs[LambdaApp[LambdaVar[1], LambdaAbs[LambdaVar[1]]]]};
    results = <|
      "countsUnbounded" ->
        (TermCount[#, 0] & /@ Range[4, 14]) === expectedCounts,
      "countsCap2" ->
        (TermCount[#, 0, 2] & /@ Range[4, 14]) === expectedCappedCounts,
      "firstTerms" ->
        (DecodeBitString[bitsForIndex[#]] & /@ Range[0, 9]) ===
          expectedFirstTerms,
      "roundTripUnbounded" ->
        AllTrue[Range[0, 400],
          EncodeLambdaTerm[DecodeBitString[bitsForIndex[#]]] ===
            bitsForIndex[#] &],
      "roundTripCap3" ->
        AllTrue[Range[0, 200],
          EncodeLambdaTerm[DecodeBitString[bitsForIndex[#], 3], 3] ===
            bitsForIndex[#] &],
      "capAgreementOnSmallSizes" ->
        AllTrue[Range[0, ClosedTermCumulative[9] - 1],
          DecodeBitString[bitsForIndex[#]] ===
            DecodeBitString[bitsForIndex[#], 8] &],
      "encodeRejectsNonClosed" ->
        AllTrue[{LambdaVar[1], LambdaAbs[LambdaVar[2]],
            LambdaApp[LambdaVar[1], LambdaVar[1]]},
          Quiet[EncodeLambdaTerm[#]] === $Failed &]|>;
    tempPath = FileNameJoin[{$TemporaryDirectory,
      "lambda-binarization-selftest-table.wl"}];
    savedValue = TermCount[14, 0];
    SaveLambdaTable[tempPath];
    ClearLambdaTable[];
    LoadLambdaTable[tempPath];
    results["saveAndLoad"] = TermCount[14, 0] === savedValue;
    DeleteFile[tempPath];
    Append[results, "allPassed" -> AllTrue[Values[results], TrueQ]]];

End[];

EndPackage[];


(* ::Section:: *)
(*Examples*)


(* ::Text:: *)
(*Evaluate the cells below after loading the package. The first group shows the start of the bijection; the second renders expression trees for a few bit strings; the third demonstrates the capped workflow and table persistence.*)


(* ::Input:: *)
(*LambdaBijectionSelfTest[]*)


(* ::Input:: *)
(*Grid[Table[{BitStringForm[bits], LambdaTermForm[DecodeBitString[bits]]}, {bits, StringDrop[IntegerString[# + 1, 2], 1] & /@ Range[0, 15]}], Alignment -> Left]*)


(* ::Input:: *)
(*LambdaTermTree[DecodeBitString["010001100"]]*)


(* ::Input:: *)
(*Row[LambdaTermTree[DecodeBitString[#]] & /@ {"01", "0110", "10110", "010001100"}, Spacer[20]]*)


(* ::Input:: *)
(*With[{omega = LambdaApp[LambdaAbs[LambdaApp[LambdaVar[1], LambdaVar[1]]], LambdaAbs[LambdaApp[LambdaVar[1], LambdaVar[1]]]]}, {LambdaTermForm[omega], EncodeLambdaTerm[omega], LambdaTermTree[omega]}]*)


(* ::Input:: *)
(*MaxDeBruijnIndex[DecodeBitString["11011000110"]]*)


(* ::Input:: *)
(*BuildLambdaTable[40, 16]*)


(* ::Input:: *)
(*SaveLambdaTable[FileNameJoin[{$HomeDirectory, "lambda-table.wl"}]]*)


(* ::Input:: *)
(*ClearLambdaTable[]; LoadLambdaTable[FileNameJoin[{$HomeDirectory, "lambda-table.wl"}]]*)
