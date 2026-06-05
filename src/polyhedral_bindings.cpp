// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
//
// WebAssembly bindings (Emscripten / embind) exposing selected entry points
// of polyhedral_common to JavaScript: copositivity, copositive factorization,
// shortest-vector realizability and automorphism group, and Gram-matrix
// canonicalization.
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
        MyMatrix<cpp_int> B = ComputeCanonicalForm<cpp_rational, cpp_int>(M, sink);
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

    emscripten::register_vector<std::string>("VectorString");

    emscripten::function("testCopositivity", &testCopositivity);
    emscripten::function("testStrictCopositivity", &testStrictCopositivity);
    emscripten::function("testCopositiveFactorization", &testCopositiveFactorization);
    emscripten::function("testShortestVectorsRealizability", &testShortestVectorsRealizability);
    emscripten::function("shortestVectorsAutomorphismGroup", &shortestVectorsAutomorphismGroup);
    emscripten::function("gramCanonicalForm", &gramCanonicalForm);
}
