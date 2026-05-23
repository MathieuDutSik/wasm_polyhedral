// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
//
// WebAssembly bindings (Emscripten / embind) exposing the
// polyhedral_common copositivity test to JavaScript.
//
// JS API (after loading the module):
//   const m = await createCoposModule();
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
#include "Copositivity.h"
// clang-format on

#include <emscripten/bind.h>
#include <emscripten/val.h>
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

// Takes the entries as a plain JS Array (via emscripten::val) so callers don't
// have to construct a typed vector on the JS side. All failure modes (bad
// input, non-symmetric matrix, parse errors) are reported via the `error`
// field of the result rather than as exceptions — keeps the JS surface
// uniform and avoids depending on the deprecated Emscripten exception
// helpers.
CopositivityResultJS testCopositivity(int n, emscripten::val const &entries_val) {
    CopositivityResultJS out;
    try {
        if (n <= 0) {
            out.error = "Dimension must be a positive integer";
            return out;
        }
        const unsigned expected = static_cast<unsigned>(n) * static_cast<unsigned>(n);
        const unsigned got = entries_val["length"].as<unsigned>();
        if (got != expected) {
            out.error = "Entries length (" + std::to_string(got) +
                        ") does not match n*n (" + std::to_string(expected) + ")";
            return out;
        }

        MyMatrix<cpp_rational> M(n, n);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                std::string s = entries_val[static_cast<unsigned>(i) * n + j].as<std::string>();
                try {
                    M(i, j) = parse_rational(s);
                } catch (std::exception const &e) {
                    out.error = "Could not parse entry A[" +
                                std::to_string(i + 1) + "," +
                                std::to_string(j + 1) + "] = \"" + s + "\"";
                    return out;
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
                    return out;
                }
            }
        }

        MyMatrix<cpp_int> InitialBasis = IdentityMat<cpp_int>(n);
        std::ostringstream sink;
        CopositivityTestResult<cpp_int> r =
            TestCopositivity<cpp_rational, cpp_int>(M, InitialBasis, sink);

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
    } catch (std::exception const &e) {
        out.error = std::string("Internal error: ") + e.what();
    } catch (...) {
        out.error = "Unknown internal error";
    }
    return out;
}

EMSCRIPTEN_BINDINGS(copos_module) {
    emscripten::value_object<CopositivityResultJS>("CopositivityResult")
        .field("isCopositive", &CopositivityResultJS::isCopositive)
        .field("nature", &CopositivityResultJS::nature)
        .field("witness", &CopositivityResultJS::witness)
        .field("error", &CopositivityResultJS::error);

    emscripten::register_vector<std::string>("VectorString");

    emscripten::function("testCopositivity", &testCopositivity);
}
