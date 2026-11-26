// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
//
// WebAssembly bindings (Emscripten / embind) exposing selected entry points
// of polyhedral_common to JavaScript: copositivity, copositive factorization,
// shortest-vector realizability and automorphism group, Gram-matrix
// canonicalization, the packing and covering densities of a polytope,
// minimal covering polytopes, and the covering radius along a Minkowski path.
//
// JS API (after loading the module):
//   const m = await createPolyhedralModule();
//   const r = m.testCopositivity(n, ["1", "0", "0", ...]);
//   r.isCopositive  : bool
//   r.nature        : string
//   r.witness       : string[]   // non-negative vector V with v^T A v < 0,
//                                // empty when isCopositive is true
//
// Matrix entries are passed as strings so the caller can express exact
// rationals like "1/2" or "-7/3" without floating-point loss.

// clang-format off
// WASM-safe header set: cpp_int / cpp_rational from boost.multiprecision are
// header-only and need no native GMP. NumberTheoryBoostGmpInt.h pulls in
// <gmp.h>, which is unavailable in the Emscripten sysroot — keep it out.
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryCommon.h"
#include "MAT_Matrix.h"
#include "Permutation.h"
#include "Group.h"
#include "Copositivity.h"
#include "StrictPositivity.h"
#include "SHORT_Realizability.h"
#include "LatticeStabEquiCan.h"
#include "POLY_DirectTriangulation.h"
#include "PolyNorm_Packing.h"
#include "PolyNorm_Covering.h"
#include "PolyNorm_MinimalCovering.h"
#include "PolyNorm_CoveringPath.h"
// clang-format on

#include <emscripten/bind.h>
#include <emscripten/val.h>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using cpp_rational = boost::multiprecision::cpp_rational;
using cpp_int = boost::multiprecision::cpp_int;

struct CopositivityResultJS {
    bool isCopositive;
    std::string nature;
    std::vector<std::string> witness;
    // Empty on success. Set to a human-readable description on any error
    // (bad input, non-symmetric matrix, etc.) so the JS layer can surface
    // a clean message without going through Emscripten exception plumbing.
    std::string error;
};

// Result of a copositive-factorization test: does A admit a representation
// A = sum_i alpha_i v_i v_i^T with v_i >= 0 (integer) and alpha_i >= 0?
// On success: `coefficients` (alpha_i) and `vectorsFlat` (v_i in row-major order).
// On failure: `certificateFlat` holds an n x n copositive matrix C such that
// <C, A> < 0, which certifies that A is not completely positive.
struct FactorizationResultJS {
    bool hasFactorization;
    int dim;                                   // n
    int nBlocks;                               // length of `coefficients`
    std::vector<std::string> coefficients;     // length nBlocks
    std::vector<std::string> vectorsFlat;      // length nBlocks * dim, row-major
    std::vector<std::string> certificateFlat;  // length dim * dim, row-major
    std::string error;
};

// Result of SHORT_TestRealizabilityShortestFamily: is the family of vectors
// in SHV realizable as the shortest-vector set of some positive-definite
// quadratic form? On success, returns the Gram matrix in `gramFlat`.
struct RealizabilityResultJS {
    bool realizable;
    int dim;                                // n (columns of SHV = ambient dim)
    std::vector<std::string> gramFlat;      // length dim * dim, row-major
    std::string error;
};

// Result of SHORT_GetStabilizer: generators of the automorphism group of
// the shortest-vector configuration as a subgroup of GL(n, Z).
struct AutomorphismResultJS {
    int dim;                                // n
    int nGenerators;
    // Flat list of nGenerators matrices of size dim*dim, each row-major.
    // generator_k(i, j) = generatorsFlat[k*dim*dim + i*dim + j].
    std::vector<std::string> generatorsFlat;
    std::string error;
};

// Result of ComputeCanonicalForm: a canonicalizing integer basis B and the
// canonical Gram matrix B * eMat * B^T.
struct CanonicalFormResultJS {
    int dim;                                  // n
    std::vector<std::string> basisFlat;       // length n*n integers (B)
    std::vector<std::string> canonicalFlat;   // length n*n rationals (B * eMat * B^T)
    std::string error;
};

// Scoped redirection of std::cerr into a caller-supplied buffer. polyhedral_common
// reports fatal conditions by writing a human-readable message to std::cerr and
// then `throw TerminalException{1};`. Without this redirect, the message lands
// on the (invisible) WASM stderr and the JS layer only sees the bare throw. With
// the redirect, the captured text is what we hand back via result.error.
//
// Lifetime: install at the top of every binding call, restore on destruction.
// Not thread-safe — fine here, the WASM main thread runs one binding at a time.
class CerrRedirect {
public:
    explicit CerrRedirect(std::ostringstream &sink)
        : prev_(std::cerr.rdbuf(sink.rdbuf())) {}
    ~CerrRedirect() { std::cerr.rdbuf(prev_); }
    CerrRedirect(CerrRedirect const &) = delete;
    CerrRedirect &operator=(CerrRedirect const &) = delete;
private:
    std::streambuf *prev_;
};

// Build the user-facing error string from whatever polyhedral_common emitted
// before throwing — both via `std::cerr` (redirected here) and via the `os`
// parameter we pass to its routines. Strips trailing newlines so it renders
// cleanly inline in the status bar.
static std::string formatPolyhedralError(std::ostringstream const &sink,
                                         char const *fallback) {
    std::string log = sink.str();
    while (!log.empty() && (log.back() == '\n' || log.back() == '\r')) log.pop_back();
    if (log.empty()) return fallback;
    return std::string("polyhedral_common error: ") + log;
}

// boost::multiprecision::cpp_int's string constructor is lenient on some
// inputs (e.g. parses "Q" as 0 silently on certain builds) so we validate
// the format ourselves before handing the string to cpp_int.
static bool is_valid_int_str(std::string const &s) {
    if (s.empty()) return false;
    size_t start = 0;
    if (s[0] == '-' || s[0] == '+') start = 1;
    if (start == s.size()) return false;
    for (size_t i = start; i < s.size(); ++i) {
        if (s[i] < '0' || s[i] > '9') return false;
    }
    return true;
}

static cpp_rational parse_rational(std::string const &s) {
    auto slash = s.find('/');
    std::string num_str = (slash == std::string::npos) ? s : s.substr(0, slash);
    std::string den_str = (slash == std::string::npos) ? std::string("1")
                                                       : s.substr(slash + 1);
    if (!is_valid_int_str(num_str) || !is_valid_int_str(den_str)) {
        throw std::invalid_argument("not a valid integer or p/q rational");
    }
    cpp_int num(num_str);
    cpp_int den(den_str);
    if (den == 0) {
        throw std::invalid_argument("denominator is zero");
    }
    return cpp_rational(num, den);
}

// Parse + validate the JS-side matrix into M. Returns true on success.
// On failure, leaves M untouched and fills out.error with a clean message.
static bool buildMatrix(int n, emscripten::val const &entries_val,
                        MyMatrix<cpp_rational> &M, CopositivityResultJS &out) {
    if (n <= 0) {
        out.error = "Dimension must be a positive integer";
        return false;
    }
    const unsigned expected = static_cast<unsigned>(n) * static_cast<unsigned>(n);
    const unsigned got = entries_val["length"].as<unsigned>();
    if (got != expected) {
        out.error = "Entries length (" + std::to_string(got) +
                    ") does not match n*n (" + std::to_string(expected) + ")";
        return false;
    }

    M.resize(n, n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            std::string s = entries_val[static_cast<unsigned>(i) * n + j].as<std::string>();
            try {
                M(i, j) = parse_rational(s);
            } catch (std::exception const &) {
                out.error = "Could not parse entry A[" +
                            std::to_string(i + 1) + "," +
                            std::to_string(j + 1) + "] = \"" + s + "\"";
                return false;
            }
        }
    }
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (M(i, j) != M(j, i)) {
                out.error = "Matrix is not symmetric: A[" +
                            std::to_string(i + 1) + "," + std::to_string(j + 1) +
                            "] != A[" + std::to_string(j + 1) + "," +
                            std::to_string(i + 1) + "]";
                return false;
            }
        }
    }
    return true;
}

// Common result-shaping for both copositivity and strict-copositivity tests.
// The boolean field name on the JS side is `isCopositive` for both modes —
// the JS layer knows which test it ran and renders the appropriate sentence.
static void fillResult(CopositivityTestResult<cpp_int> const &r,
                       CopositivityResultJS &out) {
    out.isCopositive = r.test;
    out.nature = r.strNature;
    if (!r.test) {
        out.witness.reserve(static_cast<size_t>(r.eVectResult1.size()));
        for (int i = 0; i < r.eVectResult1.size(); ++i) {
            std::ostringstream oss;
            oss << r.eVectResult1(i);
            out.witness.push_back(oss.str());
        }
    }
}

// Takes the entries as a plain JS Array (via emscripten::val) so callers don't
// have to construct a typed vector on the JS side. All failure modes (bad
// input, non-symmetric matrix, parse errors) are reported via the `error`
// field of the result rather than as exceptions.
CopositivityResultJS testCopositivity(int n, emscripten::val const &entries_val) {
    CopositivityResultJS out;
    std::ostringstream sink;
    CerrRedirect redir(sink);
    try {
        MyMatrix<cpp_rational> M;
        if (!buildMatrix(n, entries_val, M, out)) return out;
        MyMatrix<cpp_int> InitialBasis = IdentityMat<cpp_int>(n);
        fillResult(TestCopositivity<cpp_rational, cpp_int>(M, InitialBasis, sink), out);
    } catch (TerminalException const &) {
        out.error = formatPolyhedralError(sink,
            "polyhedral_common reported a fatal error during copositivity test");
    } catch (std::exception const &e) {
        out.error = std::string("Internal error: ") + e.what();
    } catch (...) {
        out.error = formatPolyhedralError(sink, "Unknown internal error");
    }
    return out;
}

CopositivityResultJS testStrictCopositivity(int n, emscripten::val const &entries_val) {
    CopositivityResultJS out;
    std::ostringstream sink;
    CerrRedirect redir(sink);
    try {
        MyMatrix<cpp_rational> M;
        if (!buildMatrix(n, entries_val, M, out)) return out;
        MyMatrix<cpp_int> InitialBasis = IdentityMat<cpp_int>(n);
        fillResult(TestStrictCopositivity<cpp_rational, cpp_int>(M, InitialBasis, sink), out);
    } catch (TerminalException const &) {
        out.error = formatPolyhedralError(sink,
            "polyhedral_common reported a fatal error during strict copositivity test");
    } catch (std::exception const &e) {
        out.error = std::string("Internal error: ") + e.what();
    } catch (...) {
        out.error = formatPolyhedralError(sink, "Unknown internal error");
    }
    return out;
}

// Test whether A admits a copositive factorization A = sum_i alpha_i v_i v_i^T
// with v_i in Z_{>=0}^n and alpha_i in Q_{>=0}. Wraps the
// TestingAttemptStrictPositivity routine from polyhedral_common/StrictPositivity.h.
FactorizationResultJS testCopositiveFactorization(int n, emscripten::val const &entries_val) {
    FactorizationResultJS out;
    out.dim = 0;
    out.nBlocks = 0;
    std::ostringstream sink;
    CerrRedirect redir(sink);
    try {
        MyMatrix<cpp_rational> M;
        // Reuse buildMatrix's parsing by going through a CopositivityResultJS
        // for the error path only.
        {
            CopositivityResultJS tmp;
            if (!buildMatrix(n, entries_val, M, tmp)) {
                out.error = tmp.error;
                return out;
            }
        }
        out.dim = n;
        MyMatrix<cpp_int> InitialBasis = IdentityMat<cpp_int>(n);
        TestStrictPositivity<cpp_rational, cpp_int> r =
            TestingAttemptStrictPositivity<cpp_rational, cpp_int>(M, InitialBasis, sink);

        out.hasFactorization = r.result;
        if (r.result) {
            const int nb = r.RealizingFamily.rows();
            out.nBlocks = nb;
            out.coefficients.reserve(static_cast<size_t>(nb));
            out.vectorsFlat.reserve(static_cast<size_t>(nb) * static_cast<size_t>(n));
            for (int i = 0; i < nb; ++i) {
                std::ostringstream c;
                c << r.ListCoeff(i);
                out.coefficients.push_back(c.str());
                for (int j = 0; j < n; ++j) {
                    std::ostringstream e;
                    e << r.RealizingFamily(i, j);
                    out.vectorsFlat.push_back(e.str());
                }
            }
        } else {
            const int rows = r.CertificateNonStrictlyPositive.rows();
            const int cols = r.CertificateNonStrictlyPositive.cols();
            out.certificateFlat.reserve(static_cast<size_t>(rows) * static_cast<size_t>(cols));
            for (int i = 0; i < rows; ++i) {
                for (int j = 0; j < cols; ++j) {
                    std::ostringstream e;
                    e << r.CertificateNonStrictlyPositive(i, j);
                    out.certificateFlat.push_back(e.str());
                }
            }
        }
    } catch (TerminalException const &) {
        out.error = formatPolyhedralError(sink,
            "polyhedral_common reported a fatal error during copositive factorization");
    } catch (std::exception const &e) {
        out.error = std::string("Internal error: ") + e.what();
    } catch (...) {
        out.error = formatPolyhedralError(sink, "Unknown internal error");
    }
    return out;
}

// Parse + validate an arbitrary-shape integer matrix from a JS Array of
// strings. Returns true on success; on failure leaves M untouched and sets
// out_error to a clean message. Used by the SHV-based bindings, whose
// inputs are non-square integer matrices.
static bool buildIntMatrix(int nRows, int nCols, emscripten::val const &entries_val,
                           MyMatrix<cpp_int> &M, std::string &out_error) {
    if (nRows <= 0 || nCols <= 0) {
        out_error = "nRows and nCols must be positive";
        return false;
    }
    const unsigned expected = static_cast<unsigned>(nRows) * static_cast<unsigned>(nCols);
    const unsigned got = entries_val["length"].as<unsigned>();
    if (got != expected) {
        out_error = "Entries length (" + std::to_string(got) +
                    ") does not match nRows*nCols (" + std::to_string(expected) + ")";
        return false;
    }

    M.resize(nRows, nCols);
    for (int i = 0; i < nRows; ++i) {
        for (int j = 0; j < nCols; ++j) {
            std::string s = entries_val[static_cast<unsigned>(i) * nCols + j].as<std::string>();
            // Reject rationals here: SHV entries must be integers.
            if (!is_valid_int_str(s)) {
                out_error = "Entry SHV[" + std::to_string(i + 1) + "," +
                            std::to_string(j + 1) + "] = \"" + s + "\" is not an integer";
                return false;
            }
            try {
                M(i, j) = cpp_int(s);
            } catch (std::exception const &) {
                out_error = "Could not parse entry SHV[" + std::to_string(i + 1) + "," +
                            std::to_string(j + 1) + "] = \"" + s + "\"";
                return false;
            }
        }
    }
    return true;
}

// Common permutalib group type used by the SHV bindings. cpp_int replaces
// the mpz_class TintGroup used in the native CP_/SHORT_/LATT_ binaries —
// permutalib is generic over the integer type and Boost cpp_int builds
// cleanly under Emscripten where GMP is unavailable.
using Tidx = uint16_t;
using Telt = permutalib::SingleSidedPerm<Tidx>;
using Tgroup = permutalib::Group<Telt, cpp_int>;

RealizabilityResultJS testShortestVectorsRealizability(
        int nRows, int nCols, emscripten::val const &entries_val) {
    RealizabilityResultJS out;
    out.dim = 0;
    std::ostringstream sink;
    CerrRedirect redir(sink);
    try {
        MyMatrix<cpp_int> SHV;
        if (!buildIntMatrix(nRows, nCols, entries_val, SHV, out.error)) return out;
        out.dim = nCols;
        // The realizability test assumes the configuration spans Q^nCols;
        // a rank-deficient input would have the inner LLL/Gram-extraction
        // path trip on an empty/degenerate matrix.
        int rank = RankMat(SHV);
        if (rank < nCols) {
            out.error = "Vector configuration is rank " + std::to_string(rank) +
                        " in dimension " + std::to_string(nCols) +
                        "; full rank (= nCols) is required.";
            return out;
        }
        ReplyRealizability<cpp_rational, cpp_int> r =
            SHORT_TestRealizabilityShortestFamily<cpp_rational, cpp_int, Tgroup>(SHV, sink);
        out.realizable = r.reply;
        if (r.reply) {
            const int n = r.eMat.rows();
            out.gramFlat.reserve(static_cast<size_t>(n) * static_cast<size_t>(n));
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    std::ostringstream e;
                    e << r.eMat(i, j);
                    out.gramFlat.push_back(e.str());
                }
            }
        }
    } catch (TerminalException const &) {
        out.error = formatPolyhedralError(sink,
            "polyhedral_common reported a fatal error during shortest-vector realizability test");
    } catch (std::exception const &e) {
        out.error = std::string("Internal error: ") + e.what();
    } catch (...) {
        out.error = formatPolyhedralError(sink, "Unknown internal error");
    }
    return out;
}

AutomorphismResultJS shortestVectorsAutomorphismGroup(
        int nRows, int nCols, emscripten::val const &entries_val) {
    AutomorphismResultJS out;
    out.dim = 0;
    out.nGenerators = 0;
    std::ostringstream sink;
    CerrRedirect redir(sink);
    try {
        MyMatrix<cpp_int> SHV;
        if (!buildIntMatrix(nRows, nCols, entries_val, SHV, out.error)) return out;
        out.dim = nCols;
        // Stabilizer computation assumes the configuration spans Q^nCols;
        // a rank-deficient input would otherwise feed a degenerate Gram into
        // the inner shortest-vector routines.
        int rank = RankMat(SHV);
        if (rank < nCols) {
            out.error = "Vector configuration is rank " + std::to_string(rank) +
                        " in dimension " + std::to_string(nCols) +
                        "; full rank (= nCols) is required.";
            return out;
        }
        std::vector<MyMatrix<cpp_int>> generators =
            SHORT_GetStabilizer<cpp_rational, cpp_int, Tgroup>(SHV, sink);
        out.nGenerators = static_cast<int>(generators.size());
        out.generatorsFlat.reserve(
            static_cast<size_t>(out.nGenerators) *
            static_cast<size_t>(nCols) * static_cast<size_t>(nCols));
        for (const auto &G : generators) {
            for (int i = 0; i < nCols; ++i) {
                for (int j = 0; j < nCols; ++j) {
                    std::ostringstream e;
                    e << G(i, j);
                    out.generatorsFlat.push_back(e.str());
                }
            }
        }
    } catch (TerminalException const &) {
        out.error = formatPolyhedralError(sink,
            "polyhedral_common reported a fatal error during shortest-vector automorphism computation");
    } catch (std::exception const &e) {
        out.error = std::string("Internal error: ") + e.what();
    } catch (...) {
        out.error = formatPolyhedralError(sink, "Unknown internal error");
    }
    return out;
}

CanonicalFormResultJS gramCanonicalForm(int n, emscripten::val const &entries_val) {
    CanonicalFormResultJS out;
    out.dim = 0;
    std::ostringstream sink;
    CerrRedirect redir(sink);
    try {
        MyMatrix<cpp_rational> M;
        {
            CopositivityResultJS tmp;
            if (!buildMatrix(n, entries_val, M, tmp)) {
                out.error = tmp.error;
                return out;
            }
        }
        out.dim = n;
        // Canonicalization is only well-defined for positive-definite Gram
        // matrices — the inner shortest-vector enumeration would otherwise
        // return an empty family and downstream Eigen access tramples OOB.
        if (!IsPositiveDefinite(M, sink)) {
            out.error = "Matrix is not positive definite; "
                        "Gram-matrix canonicalization requires a positive-definite form.";
            return out;
        }
        MyMatrix<cpp_int> B = ComputeCanonicalForm<cpp_rational, cpp_int, Tgroup>(M, sink);
        MyMatrix<cpp_rational> B_T = UniversalMatrixConversion<cpp_rational, cpp_int>(B);
        MyMatrix<cpp_rational> canonical = B_T * M * B_T.transpose();

        out.basisFlat.reserve(static_cast<size_t>(n) * static_cast<size_t>(n));
        out.canonicalFlat.reserve(static_cast<size_t>(n) * static_cast<size_t>(n));
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                std::ostringstream eb;
                eb << B(i, j);
                out.basisFlat.push_back(eb.str());
                std::ostringstream ec;
                ec << canonical(i, j);
                out.canonicalFlat.push_back(ec.str());
            }
        }
    } catch (TerminalException const &) {
        out.error = formatPolyhedralError(sink,
            "polyhedral_common reported a fatal error during Gram-matrix canonicalization");
    } catch (std::exception const &e) {
        out.error = std::string("Internal error: ") + e.what();
    } catch (...) {
        out.error = formatPolyhedralError(sink, "Unknown internal error");
    }
    return out;
}

// Result of polytopeNormDensities: the packing and covering of R^n by the
// translates alpha P + Z^n and mu P + Z^n of a polytope P. All values are
// exact rationals, as strings. The vertices are in the input coordinates.
struct PolyNormResultJS {
    int dim;
    int nVertices;
    int nFacets;
    std::vector<std::string> verticesFlat;  // nVertices * dim, row-major
    std::string volume;
    std::string groupOrder;                 // |{A in GL_n(Z) : PA = P}|
    bool hasPacking;
    std::string alpha;                      // packing scalar
    std::string packingDensity;             // alpha^n vol(P)
    int nContacts;
    std::vector<std::string> contactsFlat;  // nContacts * dim, row-major
    bool hasCovering;
    std::string mu;                         // covering radius
    std::string coveringDensity;            // mu^n vol(P)
    std::vector<std::string> lastCoveredPoint;
    int nTight;
    std::vector<std::string> tightFlat;     // nTight * dim: z on the boundary of p - mu P
    std::string error;
};

// Parse a rectangular matrix of rationals (one vector per row).
static bool buildRationalMatrix(int nRows, int nCols, emscripten::val const &entries_val,
                                MyMatrix<cpp_rational> &M, std::string &out_error) {
    if (nRows <= 0 || nCols <= 0) {
        out_error = "The number of vectors and the dimension must be positive";
        return false;
    }
    const unsigned expected = static_cast<unsigned>(nRows) * static_cast<unsigned>(nCols);
    const unsigned got = entries_val["length"].as<unsigned>();
    if (got != expected) {
        out_error = "Entries length (" + std::to_string(got) +
                    ") does not match nRows*nCols (" + std::to_string(expected) + ")";
        return false;
    }
    M.resize(nRows, nCols);
    for (int i = 0; i < nRows; ++i) {
        for (int j = 0; j < nCols; ++j) {
            std::string s = entries_val[static_cast<unsigned>(i) * nCols + j].as<std::string>();
            try {
                M(i, j) = parse_rational(s);
            } catch (std::exception const &) {
                out_error = "Could not parse entry v" + std::to_string(i + 1) + "[" +
                            std::to_string(j + 1) + "] = \"" + s + "\"";
                return false;
            }
        }
    }
    return true;
}

static MyMatrix<cpp_rational> homogenize(MyMatrix<cpp_rational> const &V) {
    MyMatrix<cpp_rational> EXT(V.rows(), V.cols() + 1);
    for (int i = 0; i < V.rows(); ++i) {
        EXT(i, 0) = 1;
        for (int j = 0; j < V.cols(); ++j) EXT(i, j + 1) = V(i, j);
    }
    return EXT;
}

// The vertices of P = {x : v_i . x <= 1}. P is bounded exactly when 0 is in
// the interior of conv(v_i), and then P is the polar of conv(v_i): its
// vertices are the facets a . y <= 1 of conv(v_i).
static bool verticesFromInequalities(MyMatrix<cpp_rational> const &V,
                                     MyMatrix<cpp_rational> &EXT, std::string &out_error,
                                     std::ostream &os) {
    int n = V.cols();
    MyMatrix<cpp_rational> EXTv = homogenize(V);
    if (RankMat(EXTv) != n + 1) {
        out_error = "The normals v_i do not span an n-dimensional polytope: "
                    "the inequalities do not define a bounded polytope";
        return false;
    }
    MyMatrix<cpp_rational> FAC =
        DirectFacetComputationInequalities_General(EXTv, "cdd", os);
    // A facet w_0 + w . y >= 0 of conv(v_i) with w_0 > 0 gives the vertex -w / w_0.
    EXT.resize(FAC.rows(), n + 1);
    for (int i = 0; i < FAC.rows(); ++i) {
        if (FAC(i, 0) <= 0) {
            out_error = "0 is not in the interior of conv(v_1, ..., v_M): "
                        "the inequalities v_i . x <= 1 define an unbounded set";
            return false;
        }
        EXT(i, 0) = 1;
        for (int j = 0; j < n; ++j) EXT(i, j + 1) = -FAC(i, j + 1) / FAC(i, 0);
    }
    return true;
}

// The vertices of the zonotope v_1 [0,1] + ... + v_M [0,1]. The generators
// are added one at a time, and once the point set is full dimensional it is
// reduced to its vertices, so that its size stays that of the zonotope.
static bool verticesFromZonotope(MyMatrix<cpp_rational> const &V,
                                 MyMatrix<cpp_rational> &EXT, std::string &out_error,
                                 std::ostream &os) {
    int n = V.cols();
    std::vector<MyVector<cpp_rational>> pts{ZeroVector<cpp_rational>(n)};
    for (int i = 0; i < V.rows(); ++i) {
        MyVector<cpp_rational> g = GetMatrixRow(V, i);
        if (IsZeroVector(g)) continue;
        std::unordered_set<MyVector<cpp_rational>> seen(pts.begin(), pts.end());
        size_t n_pt = pts.size();
        for (size_t k = 0; k < n_pt; ++k) {
            MyVector<cpp_rational> q = pts[k] + g;
            if (seen.insert(q).second) pts.push_back(q);
        }
        MyMatrix<cpp_rational> EXTp = homogenize(MatrixFromVectorFamily(pts));
        if (RankMat(EXTp) == n + 1) {
            std::vector<MyVector<cpp_rational>> red;
            for (int idx : get_distinct_extreme_rays_indices(EXTp, os)) red.push_back(pts[idx]);
            pts = std::move(red);
        }
    }
    EXT = homogenize(MatrixFromVectorFamily(pts));
    if (RankMat(EXT) != n + 1) {
        out_error = "The generators v_i do not span R^n: the zonotope is not full dimensional";
        return false;
    }
    return true;
}

// The homogeneous vertex matrix (rows (1, v)) of the polytope given by the
// rows v_i of the nRows x nCols matrix `entries`, read according to inputType:
//   "vertices"     P = conv(v_1, ..., v_M)
//   "inequalities" P = {x : v_i . x <= 1}
//   "zonotope"     P = v_1 [0,1] + ... + v_M [0,1]
// For "vertices" the rows that are repeated or not vertices are kept; the
// callers drop them.
static bool buildPolytope(std::string const &inputType, int nRows, int nCols,
                          emscripten::val const &entries_val,
                          MyMatrix<cpp_rational> &EXT, std::string &out_error,
                          std::ostream &os) {
    MyMatrix<cpp_rational> V;
    if (!buildRationalMatrix(nRows, nCols, entries_val, V, out_error)) return false;
    if (inputType == "vertices") {
        EXT = homogenize(V);
        if (RankMat(EXT) != nCols + 1) {
            out_error = "The points do not span an n-dimensional polytope";
            return false;
        }
        return true;
    }
    if (inputType == "inequalities") return verticesFromInequalities(V, EXT, out_error, os);
    if (inputType == "zonotope") return verticesFromZonotope(V, EXT, out_error, os);
    out_error = "Unknown input type \"" + inputType +
                "\", allowed: vertices, inequalities, zonotope";
    return false;
}

// The rows of a homogeneous matrix without their first coordinate, flattened.
static std::vector<std::string> flattenAffine(MyMatrix<cpp_rational> const &EXT) {
    std::vector<std::string> flat;
    for (int i = 0; i < EXT.rows(); ++i) {
        for (int j = 1; j < EXT.cols(); ++j) {
            std::ostringstream oss;
            oss << EXT(i, j);
            flat.push_back(oss.str());
        }
    }
    return flat;
}

static std::vector<std::string> flattenMatrix(MyMatrix<cpp_rational> const &M) {
    std::vector<std::string> flat;
    for (int i = 0; i < M.rows(); ++i) {
        for (int j = 0; j < M.cols(); ++j) {
            std::ostringstream oss;
            oss << M(i, j);
            flat.push_back(oss.str());
        }
    }
    return flat;
}

static std::string toStr(cpp_rational const &x) {
    std::ostringstream oss;
    oss << x;
    return oss.str();
}

static cpp_rational powRational(cpp_rational const &x, int n) {
    cpp_rational r(1);
    for (int i = 0; i < n; ++i) r *= x;
    return r;
}

// Packing and covering of R^n by the translates of the polytope P along Z^n.
//   inputType : "vertices"     P = conv(v_1, ..., v_M)
//               "inequalities" P = {x : v_i . x <= 1}
//               "zonotope"     P = v_1 [0,1] + ... + v_M [0,1]
//   what      : "both", "packing" or "covering"
// The rows of the nRows x nCols matrix `entries` are the v_i.
PolyNormResultJS polytopeNormDensities(std::string const &inputType, int nRows, int nCols,
                                       emscripten::val const &entries_val,
                                       std::string const &what) {
    PolyNormResultJS out;
    out.dim = 0;
    out.nVertices = 0;
    out.nFacets = 0;
    out.hasPacking = false;
    out.nContacts = 0;
    out.hasCovering = false;
    out.nTight = 0;
    std::ostringstream sink;
    CerrRedirect redir(sink);
    try {
        bool doPacking = what == "both" || what == "packing";
        bool doCovering = what == "both" || what == "covering";
        if (!doPacking && !doCovering) {
            out.error = "Unknown computation \"" + what + "\", allowed: both, packing, covering";
            return out;
        }
        int n = nCols;
        MyMatrix<cpp_rational> EXT;
        if (!buildPolytope(inputType, nRows, nCols, entries_val, EXT, out.error, sink))
            return out;
        PolyNormData<cpp_rational, cpp_int, Tgroup> data =
            PolyNorm_BuildData<cpp_rational, cpp_int, Tgroup>(EXT, sink);
        out.dim = n;
        out.nVertices = data.EXT.rows();
        out.nFacets = data.Lmat.rows();
        for (int i = 0; i < data.EXT.rows(); ++i) {
            for (int j = 0; j < n; ++j) {
                out.verticesFlat.push_back(toStr(data.EXT(i, j + 1) + data.center(j)));
            }
        }
        cpp_rational vol =
            DirectVolumePolytope(data.EXT, TriangulationProgram::lrs, sink);
        out.volume = toStr(vol);
        {
            std::ostringstream oss;
            oss << data.GRP.size();
            out.groupOrder = oss.str();
        }
        if (doPacking) {
            PolyNormPacking<cpp_rational, cpp_int> res =
                ComputePolyNormPacking<cpp_rational, cpp_int>(data.LmatDiff, sink);
            out.hasPacking = true;
            out.alpha = toStr(res.alpha);
            out.packingDensity = toStr(powRational(res.alpha, n) * vol);
            out.nContacts = static_cast<int>(res.ListContact.size());
            for (auto const &z : res.ListContact) {
                for (int j = 0; j < n; ++j) {
                    std::ostringstream oss;
                    oss << z(j);
                    out.contactsFlat.push_back(oss.str());
                }
            }
        }
        if (doCovering) {
            PolyNormCovering<cpp_rational, cpp_int> res =
                ComputePolyNormCovering<cpp_rational, cpp_int, Tgroup>(data, sink);
            out.hasCovering = true;
            out.mu = toStr(res.mu);
            out.coveringDensity = toStr(powRational(res.mu, n) * vol);
            for (int j = 0; j < n; ++j) {
                out.lastCoveredPoint.push_back(toStr(res.p(j) + data.center(j)));
            }
            out.nTight = static_cast<int>(res.ListTight.size());
            for (auto const &ePair : res.ListTight) {
                for (int j = 0; j < n; ++j) {
                    std::ostringstream oss;
                    oss << ePair.z(j);
                    out.tightFlat.push_back(oss.str());
                }
            }
        }
    } catch (TerminalException const &) {
        out.error = formatPolyhedralError(sink,
            "polyhedral_common reported a fatal error during the packing / covering computation");
    } catch (std::exception const &e) {
        out.error = std::string("Internal error: ") + e.what();
    } catch (...) {
        out.error = formatPolyhedralError(sink, "Unknown internal error");
    }
    return out;
}

// Result of isMinimalCovering: whether P is covering (P + Z^n = R^n) and
// minimal covering (no convex body strictly contained in P is covering),
// with the vertices of P that can be cut keeping the covering.
struct MinimalCoveringTestResultJS {
    int dim;
    int nVertices;
    int nFacets;
    std::vector<std::string> verticesFlat;   // nVertices * dim, row-major
    std::string volume;
    bool isCovering;
    bool isMinimal;
    std::string mu;                          // covering radius, set when not covering
    int nCuttable;
    std::vector<std::string> cuttableFlat;   // nCuttable * dim, row-major
    std::string error;
};

// Result of minimalCoveringPolytope: a covering polytope Q contained in
// mu(P) P, and whether it is proved minimal covering.
struct MinimalCoveringResultJS {
    int dim;
    std::string muInput;                     // covering radius of P
    bool isCovering;                         // mu(P) <= 1
    bool isMinimal;
    bool isTile;                             // vol(Q) = 1 (grid method)
    std::string description;                 // statistics of the method
    int nVertices;
    std::vector<std::string> verticesFlat;   // nVertices * dim, row-major
    int nFacets;
    std::vector<std::string> facetsFlat;     // nFacets * (dim + 1): (b, a), b + a.x >= 0
    std::string volume;
    std::string error;
};

MinimalCoveringTestResultJS isMinimalCovering(std::string const &inputType, int nRows,
                                              int nCols, emscripten::val const &entries_val) {
    MinimalCoveringTestResultJS out;
    out.dim = 0;
    out.nVertices = 0;
    out.nFacets = 0;
    out.isCovering = false;
    out.isMinimal = false;
    out.nCuttable = 0;
    std::ostringstream sink;
    CerrRedirect redir(sink);
    try {
        MyMatrix<cpp_rational> EXT_input;
        if (!buildPolytope(inputType, nRows, nCols, entries_val, EXT_input, out.error, sink))
            return out;
        MyMatrix<cpp_rational> EXT = SelectRow(
            EXT_input, get_distinct_extreme_rays_indices(EXT_input, sink));
        MyMatrix<cpp_rational> FAC = PolyNorm_MinCovFacets(EXT, sink);
        PolyNormIsMinimalCovering<cpp_rational> res =
            PolyNorm_IsMinimalCovering<cpp_rational, cpp_int, Tgroup>(EXT, FAC, sink);
        out.dim = nCols;
        out.nVertices = EXT.rows();
        out.nFacets = FAC.rows();
        out.verticesFlat = flattenAffine(EXT);
        out.volume = toStr(PolyNorm_MinCovVolume(EXT, sink));
        out.isCovering = res.is_covering;
        out.isMinimal = res.is_minimal;
        if (!res.is_covering) {
            out.mu = toStr(PolyNorm_MinCovCoveringRadius<cpp_rational, cpp_int, Tgroup>(EXT, sink));
        }
        out.nCuttable = static_cast<int>(res.ListCuttable.size());
        for (int iv : res.ListCuttable) {
            for (int j = 0; j < nCols; ++j) out.cuttableFlat.push_back(toStr(EXT(iv, j + 1)));
        }
    } catch (TerminalException const &) {
        out.error = formatPolyhedralError(sink,
            "polyhedral_common reported a fatal error during the minimal covering test");
    } catch (std::exception const &e) {
        out.error = std::string("Internal error: ") + e.what();
    } catch (...) {
        out.error = formatPolyhedralError(sink, "Unknown internal error");
    }
    return out;
}

// A minimal covering polytope inside mu(P) P.
//   method : "grid"      param1 = N >= 1, param2 = n_round >= 1
//            "reduction" param1 = the maximal number of reductions
MinimalCoveringResultJS minimalCoveringPolytope(std::string const &inputType, int nRows,
                                                int nCols, emscripten::val const &entries_val,
                                                std::string const &method, int param1,
                                                int param2) {
    MinimalCoveringResultJS out;
    out.dim = 0;
    out.isCovering = false;
    out.isMinimal = false;
    out.isTile = false;
    out.nVertices = 0;
    out.nFacets = 0;
    std::ostringstream sink;
    CerrRedirect redir(sink);
    try {
        MyMatrix<cpp_rational> EXT;
        if (!buildPolytope(inputType, nRows, nCols, entries_val, EXT, out.error, sink))
            return out;
        out.dim = nCols;
        std::ostringstream desc;
        MyMatrix<cpp_rational> EXTres, FACres;
        if (method == "grid") {
            if (param1 < 1 || param2 < 1) {
                out.error = "The grid method needs N >= 1 and at least one round";
                return out;
            }
            PolyNormMinimalCoveringGrid<cpp_rational> res =
                PolyNorm_MinimalCoveringPolytopeGrid<cpp_rational, cpp_int, Tgroup>(
                    EXT, param1, param2, sink);
            out.muInput = toStr(res.mu_input);
            out.isCovering = res.is_covering;
            if (!res.is_covering) return out;
            out.isMinimal = res.is_minimal;
            out.isTile = res.is_tile;
            desc << "Grid method: " << res.n_round << " round(s), last grid (1/"
                 << res.N_last << ") Z^n with " << res.n_point << " candidate points, "
                 << res.n_drop << " points dropped, " << res.n_covering_test
                 << " covering tests, " << res.n_lp
                 << " linear programs for the local tests of the vertices.";
            EXTres = res.EXT;
            FACres = res.FAC;
        } else if (method == "reduction") {
            if (param1 < 1) {
                out.error = "The reduction method needs at least one reduction step";
                return out;
            }
            PolyNormMinimalCovering<cpp_rational> res =
                PolyNorm_MinimalCoveringPolytope<cpp_rational, cpp_int, Tgroup>(
                    EXT, static_cast<size_t>(param1), sink);
            out.muInput = toStr(res.mu_input);
            out.isCovering = res.is_covering;
            if (!res.is_covering) return out;
            out.isMinimal = res.is_minimal;
            desc << "Reductions after the shrinking by mu: " << res.n_remove
                 << " vertices removed, " << res.n_push << " facets pushed, "
                 << res.n_cut << " vertices cut; " << res.n_node << " nodes, "
                 << res.n_lp << " linear programs.";
            EXTres = res.EXT;
            FACres = res.FAC;
        } else {
            out.error = "Unknown method \"" + method + "\", allowed: grid, reduction";
            return out;
        }
        out.description = desc.str();
        out.nVertices = EXTres.rows();
        out.verticesFlat = flattenAffine(EXTres);
        out.nFacets = FACres.rows();
        out.facetsFlat = flattenMatrix(FACres);
        out.volume = toStr(PolyNorm_MinCovVolume(EXTres, sink));
    } catch (TerminalException const &) {
        out.error = formatPolyhedralError(sink,
            "polyhedral_common reported a fatal error during the minimal covering computation");
    } catch (std::exception const &e) {
        out.error = std::string("Internal error: ") + e.what();
    } catch (...) {
        out.error = formatPolyhedralError(sink, "Unknown internal error");
    }
    return out;
}

// Result of coveringPath: the function f(lambda) = 1 / mu(lambda K1 +
// (1 - lambda) K2) on [0,1], continuous and piecewise affine, as its pieces
// f = s0 + s1 lambda on [alpha, beta], consecutive and each one proved.
struct CoveringPathResultJS {
    int dim;
    int nPieces;
    std::vector<std::string> piecesFlat;     // nPieces * 4: alpha, beta, s0, s1
    std::string description;                 // statistics of the search
    std::string error;
};

// The covering radius along the Minkowski path between K1 and K2, given as
// for polytopeNormDensities (same dimension nCols). PolyNorm_ComputeCoveringPath
// computes 1 / mu((1 - lambda) P1 + lambda P2), so P1 = K2 and P2 = K1.
CoveringPathResultJS coveringPath(std::string const &inputType1, int nRows1,
                                  emscripten::val const &entries1,
                                  std::string const &inputType2, int nRows2,
                                  emscripten::val const &entries2, int nCols) {
    CoveringPathResultJS out;
    out.dim = 0;
    out.nPieces = 0;
    std::ostringstream sink;
    CerrRedirect redir(sink);
    try {
        MyMatrix<cpp_rational> EXT1, EXT2;
        if (!buildPolytope(inputType1, nRows1, nCols, entries1, EXT1, out.error, sink)) {
            out.error = "K1: " + out.error;
            return out;
        }
        if (!buildPolytope(inputType2, nRows2, nCols, entries2, EXT2, out.error, sink)) {
            out.error = "K2: " + out.error;
            return out;
        }
        PolyNormCoveringPath<cpp_rational> res =
            PolyNorm_ComputeCoveringPath<cpp_rational, cpp_int, Tgroup>(EXT2, EXT1, sink);
        out.dim = nCols;
        out.nPieces = static_cast<int>(res.ListPiece.size());
        for (auto const &P : res.ListPiece) {
            out.piecesFlat.push_back(toStr(P.alpha));
            out.piecesFlat.push_back(toStr(P.beta));
            out.piecesFlat.push_back(toStr(P.s0));
            out.piecesFlat.push_back(toStr(P.s1));
        }
        std::ostringstream desc;
        desc << "Search: " << res.n_evaluation << " exact evaluations, "
             << res.n_certificate << " certificates, " << res.n_node << " nodes, "
             << res.n_lp << " linear programs.";
        out.description = desc.str();
    } catch (TerminalException const &) {
        out.error = formatPolyhedralError(sink,
            "polyhedral_common reported a fatal error during the covering path computation");
    } catch (std::exception const &e) {
        out.error = std::string("Internal error: ") + e.what();
    } catch (...) {
        out.error = formatPolyhedralError(sink, "Unknown internal error");
    }
    return out;
}

EMSCRIPTEN_BINDINGS(polyhedral_module) {
    emscripten::value_object<CopositivityResultJS>("CopositivityResult")
        .field("isCopositive", &CopositivityResultJS::isCopositive)
        .field("nature", &CopositivityResultJS::nature)
        .field("witness", &CopositivityResultJS::witness)
        .field("error", &CopositivityResultJS::error);

    emscripten::value_object<FactorizationResultJS>("FactorizationResult")
        .field("hasFactorization", &FactorizationResultJS::hasFactorization)
        .field("dim", &FactorizationResultJS::dim)
        .field("nBlocks", &FactorizationResultJS::nBlocks)
        .field("coefficients", &FactorizationResultJS::coefficients)
        .field("vectorsFlat", &FactorizationResultJS::vectorsFlat)
        .field("certificateFlat", &FactorizationResultJS::certificateFlat)
        .field("error", &FactorizationResultJS::error);

    emscripten::value_object<RealizabilityResultJS>("RealizabilityResult")
        .field("realizable", &RealizabilityResultJS::realizable)
        .field("dim", &RealizabilityResultJS::dim)
        .field("gramFlat", &RealizabilityResultJS::gramFlat)
        .field("error", &RealizabilityResultJS::error);

    emscripten::value_object<AutomorphismResultJS>("AutomorphismResult")
        .field("dim", &AutomorphismResultJS::dim)
        .field("nGenerators", &AutomorphismResultJS::nGenerators)
        .field("generatorsFlat", &AutomorphismResultJS::generatorsFlat)
        .field("error", &AutomorphismResultJS::error);

    emscripten::value_object<CanonicalFormResultJS>("CanonicalFormResult")
        .field("dim", &CanonicalFormResultJS::dim)
        .field("basisFlat", &CanonicalFormResultJS::basisFlat)
        .field("canonicalFlat", &CanonicalFormResultJS::canonicalFlat)
        .field("error", &CanonicalFormResultJS::error);

    emscripten::value_object<PolyNormResultJS>("PolyNormResult")
        .field("dim", &PolyNormResultJS::dim)
        .field("nVertices", &PolyNormResultJS::nVertices)
        .field("nFacets", &PolyNormResultJS::nFacets)
        .field("verticesFlat", &PolyNormResultJS::verticesFlat)
        .field("volume", &PolyNormResultJS::volume)
        .field("groupOrder", &PolyNormResultJS::groupOrder)
        .field("hasPacking", &PolyNormResultJS::hasPacking)
        .field("alpha", &PolyNormResultJS::alpha)
        .field("packingDensity", &PolyNormResultJS::packingDensity)
        .field("nContacts", &PolyNormResultJS::nContacts)
        .field("contactsFlat", &PolyNormResultJS::contactsFlat)
        .field("hasCovering", &PolyNormResultJS::hasCovering)
        .field("mu", &PolyNormResultJS::mu)
        .field("coveringDensity", &PolyNormResultJS::coveringDensity)
        .field("lastCoveredPoint", &PolyNormResultJS::lastCoveredPoint)
        .field("nTight", &PolyNormResultJS::nTight)
        .field("tightFlat", &PolyNormResultJS::tightFlat)
        .field("error", &PolyNormResultJS::error);

    emscripten::value_object<MinimalCoveringTestResultJS>("MinimalCoveringTestResult")
        .field("dim", &MinimalCoveringTestResultJS::dim)
        .field("nVertices", &MinimalCoveringTestResultJS::nVertices)
        .field("nFacets", &MinimalCoveringTestResultJS::nFacets)
        .field("verticesFlat", &MinimalCoveringTestResultJS::verticesFlat)
        .field("volume", &MinimalCoveringTestResultJS::volume)
        .field("isCovering", &MinimalCoveringTestResultJS::isCovering)
        .field("isMinimal", &MinimalCoveringTestResultJS::isMinimal)
        .field("mu", &MinimalCoveringTestResultJS::mu)
        .field("nCuttable", &MinimalCoveringTestResultJS::nCuttable)
        .field("cuttableFlat", &MinimalCoveringTestResultJS::cuttableFlat)
        .field("error", &MinimalCoveringTestResultJS::error);

    emscripten::value_object<MinimalCoveringResultJS>("MinimalCoveringResult")
        .field("dim", &MinimalCoveringResultJS::dim)
        .field("muInput", &MinimalCoveringResultJS::muInput)
        .field("isCovering", &MinimalCoveringResultJS::isCovering)
        .field("isMinimal", &MinimalCoveringResultJS::isMinimal)
        .field("isTile", &MinimalCoveringResultJS::isTile)
        .field("description", &MinimalCoveringResultJS::description)
        .field("nVertices", &MinimalCoveringResultJS::nVertices)
        .field("verticesFlat", &MinimalCoveringResultJS::verticesFlat)
        .field("nFacets", &MinimalCoveringResultJS::nFacets)
        .field("facetsFlat", &MinimalCoveringResultJS::facetsFlat)
        .field("volume", &MinimalCoveringResultJS::volume)
        .field("error", &MinimalCoveringResultJS::error);

    emscripten::value_object<CoveringPathResultJS>("CoveringPathResult")
        .field("dim", &CoveringPathResultJS::dim)
        .field("nPieces", &CoveringPathResultJS::nPieces)
        .field("piecesFlat", &CoveringPathResultJS::piecesFlat)
        .field("description", &CoveringPathResultJS::description)
        .field("error", &CoveringPathResultJS::error);

    emscripten::register_vector<std::string>("VectorString");

    emscripten::function("testCopositivity", &testCopositivity);
    emscripten::function("testStrictCopositivity", &testStrictCopositivity);
    emscripten::function("testCopositiveFactorization", &testCopositiveFactorization);
    emscripten::function("testShortestVectorsRealizability", &testShortestVectorsRealizability);
    emscripten::function("shortestVectorsAutomorphismGroup", &shortestVectorsAutomorphismGroup);
    emscripten::function("gramCanonicalForm", &gramCanonicalForm);
    emscripten::function("polytopeNormDensities", &polytopeNormDensities);
    emscripten::function("isMinimalCovering", &isMinimalCovering);
    emscripten::function("minimalCoveringPolytope", &minimalCoveringPolytope);
    emscripten::function("coveringPath", &coveringPath);
}
