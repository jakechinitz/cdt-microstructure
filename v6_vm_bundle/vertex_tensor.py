# ============================================================================
#  VERTEX TENSOR MODULE: pluggable vertex amplitude for the condensate
#
#  PURPOSE: be the single seam between the (eventual) faithful EPRL vertex from
#  sl2cfoam-next and the condensate Monte Carlo (Link B). Both implementations
#  expose the same interface, so switching between them is a one-line change.
#
#  TWO IMPLEMENTATIONS:
#    * FaithfulVertex  -- wraps the (7,7,7,7,7) tensor produced by
#                         sl2cfoam_driver.jl. Refuses to load if the saved
#                         tensor failed the asymptotic validation gate.
#    * SchematicVertex -- reproduces Link B's existing real-6j-derived weights.
#                         Pipeline-tests the scaffold WITHOUT the real tensor,
#                         but the pipeline runner flags any result built from
#                         it as ARTIFACT/UNTRUSTED.
#
#  INVARIANT THE PIPELINE RELIES ON:
#    a number measured through SchematicVertex MUST NOT be reported as a result
#    of the theory. The Vertex.is_faithful flag carries this and the runner
#    consumes it. Do not flip the flag without a real validated tensor.
# ============================================================================

from __future__ import annotations
import os
from dataclasses import dataclass
from typing import Sequence
import numpy as np

# ----------------------------------------------------------------------------
#  Base interface
# ----------------------------------------------------------------------------
class Vertex:
    """Abstract vertex amplitude. Implementations expose amplitude() and dim."""
    is_faithful: bool = False
    label: str = "abstract"
    dim: int = 0

    def amplitude(self, intws: Sequence[int]) -> float:
        """Return A_v(i1,...,i5). intws must have length 5."""
        raise NotImplementedError

    def dense_tensor(self) -> np.ndarray:
        """Return the full (D,D,D,D,D) amplitude array (for numba interop).
        Default: synthesize by sampling amplitude() over the grid.
        FaithfulVertex overrides to return its stored array directly."""
        D = self.dim
        T = np.empty((D, D, D, D, D), dtype=np.float64)
        for a in range(D):
            for b in range(D):
                for c in range(D):
                    for d in range(D):
                        for e in range(D):
                            T[a, b, c, d, e] = self.amplitude((a, b, c, d, e))
        return T

    def describe(self) -> str:
        return f"<Vertex {self.label} dim={self.dim} faithful={self.is_faithful}>"

# ----------------------------------------------------------------------------
#  Faithful EPRL vertex from sl2cfoam-next
# ----------------------------------------------------------------------------
@dataclass
class FaithfulVertex(Vertex):
    """Wraps the (D,D,D,D,D) tensor written by sl2cfoam_driver.jl.

    Refuses to load unless the saved 'validated' flag is 1. Refuses to load if
    the tensor fails basic sanity (finite, real, expected shape).
    """
    tensor: np.ndarray
    j: int
    immirzi: float
    shells: int
    y_map: str

    is_faithful = True
    label = "faithful (sl2cfoam-next)"

    def __post_init__(self):
        self.dim = self.tensor.shape[0]

    @classmethod
    def load(cls, path: str) -> "FaithfulVertex":
        if not os.path.exists(path):
            raise FileNotFoundError(
                f"vertex tensor not found at {path}. "
                f"Run sl2cfoam_driver.jl first to produce it, or fall back "
                f"to SchematicVertex (pipeline will mark output as artifact)."
            )
        data = np.load(path, allow_pickle=False)
        validated = int(data["validated"]) if "validated" in data.files else 0
        if validated != 1:
            raise ValueError(
                f"vertex tensor at {path} is marked UNVALIDATED (asymptotic "
                f"gate failed in the Julia driver). Refusing to load -- "
                f"investigate vertex_validation.log before using this tensor."
            )
        T = np.asarray(data["tensor"], dtype=np.float64)
        if T.ndim != 5 or len(set(T.shape)) != 1:
            raise ValueError(f"vertex tensor has wrong shape {T.shape}; "
                             f"expected (D,D,D,D,D).")
        if not np.all(np.isfinite(T)):
            raise ValueError("vertex tensor contains non-finite entries; "
                             "library wiring is broken.")
        return cls(
            tensor=T,
            j=int(data["j"]),
            immirzi=float(data["immirzi"]),
            shells=int(data["shells"]),
            y_map=str(data["y_map"]) if "y_map" in data.files else "unknown",
        )

    def amplitude(self, intws: Sequence[int]) -> float:
        if len(intws) != 5:
            raise ValueError(f"expected 5 intertwiner labels, got {len(intws)}")
        D = self.dim
        ii = tuple(min(max(int(i), 0), D-1) for i in intws)
        return float(self.tensor[ii])

    def dense_tensor(self) -> np.ndarray:
        return self.tensor

# ----------------------------------------------------------------------------
#  Schematic fallback: reproduces Link B's current per-tetra-6j behavior
# ----------------------------------------------------------------------------
class SchematicVertex(Vertex):
    """Real per-tetra 6j-derived weights with an APPROXIMATE 5-tetra gluing.

    Reproduces the existing Link B behavior. Lets the pipeline run end-to-end
    without the faithful vertex so the scaffold itself is testable. Pipeline
    runner MUST flag any result built from this as artifact.
    """
    is_faithful = False
    label = "schematic (real 6j, approximate gluing)"

    _WEIGHTS = np.array(
        [0.071429, 0.196429, 0.113095, 0.25, 0.314935, 0.053571, 0.000541]
    )

    def __init__(self, dim: int = 7):
        self.dim = dim
        if dim != len(self._WEIGHTS):
            raise ValueError(f"SchematicVertex hardcoded for dim=7, got {dim}")

    def amplitude(self, intws: Sequence[int]) -> float:
        if len(intws) != 5:
            raise ValueError(f"expected 5 intertwiner labels, got {len(intws)}")
        D = self.dim
        ii = [min(max(int(i), 0), D-1) for i in intws]
        imb = abs(sum(i - (D//2) for i in ii))
        prod = 1.0
        for i in ii:
            prod *= max(self._WEIGHTS[i], 1e-9)
        return prod ** 0.2 * np.exp(-0.15 * imb)

# ----------------------------------------------------------------------------
#  Convenience: auto-pick faithful if present, else schematic with loud notice
# ----------------------------------------------------------------------------
def load_or_fallback(path: str = "vertex_j3.npz", verbose: bool = True) -> Vertex:
    """Load FaithfulVertex if available and validated; otherwise SchematicVertex.

    The pipeline runner uses Vertex.is_faithful to decide whether the output
    is allowed to be reported as a result.
    """
    try:
        v = FaithfulVertex.load(path)
        if verbose:
            print(f"[vertex] loaded FAITHFUL vertex from {path}  ({v.describe()})")
        return v
    except FileNotFoundError as e:
        if verbose:
            print(f"[vertex] {e}")
            print(f"[vertex] FALLBACK: using SchematicVertex. Any d_s reported "
                  f"will be flagged ARTIFACT until the real tensor is plugged in.")
        return SchematicVertex()
    except ValueError as e:
        if verbose:
            print(f"[vertex] ERROR loading {path}: {e}")
            print(f"[vertex] FALLBACK: using SchematicVertex.")
        return SchematicVertex()


def validate_tensor_structure(T: np.ndarray) -> dict:
    """Run structural checks. Returns a report dict; does not raise."""
    report = {
        "shape": T.shape,
        "finite": bool(np.all(np.isfinite(T))),
        "real": bool(np.all(np.isreal(T))),
        "max_abs": float(np.max(np.abs(T))),
        "min_nonzero_abs": float(np.min(np.abs(T[T != 0]))) if np.any(T != 0) else 0.0,
        "nonzero_frac": float(np.mean(T != 0)),
        "rank5": T.ndim == 5,
        "cubic": len(set(T.shape)) == 1,
    }
    if T.ndim == 5 and len(set(T.shape)) == 1:
        Tswap = np.swapaxes(T, 0, 1)
        denom = max(float(np.linalg.norm(T)), 1e-12)
        report["swap01_relerr"] = float(np.linalg.norm(T - Tswap) / denom)
    return report


if __name__ == "__main__":
    print("vertex_tensor.py self-test\n")
    print("1) SchematicVertex round-trip:")
    sv = SchematicVertex()
    print(f"   {sv.describe()}")
    print(f"   amplitude((3,3,3,3,3)) = {sv.amplitude((3,3,3,3,3)):.4e}")
    Tsv = sv.dense_tensor()
    print(f"   dense_tensor shape: {Tsv.shape}, max={Tsv.max():.4e}")
    print()
    print("2) load_or_fallback (will pick faithful if vertex_j3.npz present):")
    v = load_or_fallback("vertex_j3.npz")
    print(f"   chose: {v.describe()}")
    if v.is_faithful:
        T = v.dense_tensor()
        print(f"   faithful dense_tensor: shape={T.shape}, max={T.max():.4e}")
