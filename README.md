# A Faster Software Implementation of SQIsign
## This code is an attachment of the manuscript "A Faster Software Implementation of SQIsign", utilizing several techniques to optimize the implementation of SQIsign.

The code is based on the [code](https://github.com/SQISign/the-sqisign) proposed by Luca De Feo et al. For more details of their implementation see the following original readme.

We only complie and benchmark the code in the case that level=lvl1.

The default compilation option does not adapt the precomputation technique as we proposed in Section 4.3 of our manuscript. This technique may enlarge the cost of key generation, but reduces the signing cost.
To compile the code with the precomputation technique use CMake build option
```
-DPRECOMPUTED=ON
```

**The following is the original readme.**


# SQIsign

This library is a C implementation of SQIsign, short for Short Quaternion and Isogeny Signature (from isogeny graphs of supersingular elliptic curves).

## Requirements

- CMake (version 3.5 or later)
- C99-compatible compiler
- Valgrind (for dynamic testing)
- Clang static analyzer (version 10 or later, for static analysis)
- GMP (version 6.1.2 or later)

## Build

- `mkdir -p build`
- `cd build`
- `cmake -DSQISIGN_BUILD_TYPE=<ref/broadwell> ..`
- `make`

## Build options

CMake build options can be specified with `-D<BUILD_OPTION>=<VALUE>`.

### ENABLE_TESTS

Builds a test harness for the library, the default value is `ON`.

### ENABLE_CT_TESTING

Builds the library with instrumentation for constant-time behavior testing, the default value is `OFF`. Valgrind development files are used for this build option.

### ENABLE_GMP_BUILD

If set to `OFF` (by default), the gmp library on the system is dynamically linked.
If set to `ON`, a custom gmp library is linked, which is built as part of the overall build process. 

In the latter case, the following further options are available:
- `ENABLE_GMP_STATIC`: Does static linking against gmp. The default is `OFF`.
- `GMP_BUILD_CONFIG_ARGS`: Provides additional config arguments for the gmp build (for example `--disable-assembly`). By default, no config arguments are provided.

### ENABLE_DOC_TARGET
If set to `ON`, a doc target is available that builds API documentation. Note that a doxygen installation is required if set to `ON`.

The default is `OFF`.

### SQISIGN_BUILD_TYPE

Specifies the build type for which SQIsign is built. The currently supported flags are:
- `ref`, which builds the plain C reference implementation.
- `broadwell`, which builds an additional implementation with GF assembly optimized code for the Intel Broadwell architecture.

### SQISIGN_TEST_REPS

Specifies and overrides the number of (self-)test repetitions to be run.

### CMAKE_BUILD_TYPE

Can be used to specify special build types. The options are:

- `Release`: Builds with optimizations enabled and assertions disabled.
- `Debug`: Builds with debug symbols.
- `ASAN`: Builds with AddressSanitizer memory error detector.
- `MSAN`: Builds with MemorySanitizer detector for uninitialized reads.
- `LSAN`: Builds with LeakSanitizer for run-time memory leak detection.
- `UBSAN`: Builds with UndefinedBehaviorSanitizer for undefined behavior detection.

The default build type uses the flags `-O3 -Wstrict-prototypes -Wno-error=strict-prototypes -fvisibility=hidden -Wno-error=implicit-function-declaration -Wno-error=attributes`. (Notice that assertions remain enabled in this configuration, which harms performance.)


## Build artifacts

The following libraries are built:

- `libsqisign_common_sys.a`: library with common crypto - AES, Keccak and system random number generator.
- `libsqisign_common_test.a`: library with common crypto for deterministic tests - AES, Keccak and CTR-DRBG PRNG.
- `libsqisign_<level>.a`: library for `SQIsign_<level>`.
- `libsqisign_<level>_test`: library for `SQIsign_<level>`, only for test, using the deterministic CTR-DRBG as backend.
- `libsqisign_<level>_nistapi.a`: library for `SQIsign_<level>` against the NIST API.
- `libsqisign_<level>_nistapi_test.a`: library for `SQIsign_<level>` against the NIST API. Only for test, using the deterministic CTR-DRBG as backend.
- `libsqisign_gf_<level>.a`: gf sub-library, generic or for `<level>`
- `libsqisign_ec_<level>.a`: ec sub-library, generic or for `<level>`
- `libsqisign_klpt_<level>.a`: klpt sub-library, generic or for `<level>`
- `libsqisign_intbig_generic.a`: intbig sub-library, generic
- `libsqisign_quaternion_generic.a`: quaternion sub-library, generic
- `libsqisign_id2iso_<level>.a`: id2iso sub-library, generic or for `<level>`

The following test apps are built:
- `sqisign_bench_<level>`: Benchmarking suites.
- `sqisign_test_kat_<level>`: KAT test suites.
- `sqisign_test_scheme_<level>`: Self-test suites.
- `sqisign_test_prof_<level>`: Profiling suites.

More apps are built in folder `build/apps`:

- `PQCgenKAT_sign_<param>`: App for generating NIST KAT.
- `example_nistapi_<param>`: Example app using the NIST API.

## Test

In the build directory, run: `make test` or `ctest`.

The test harness consists of the following units:

- KAT test: tests against the KAT files in the `KAT` folder - `SQIsign_<level>_KAT`
- Self-tests: runs random self-tests (key-generation, signing and verifying) - `SQIsign_<level>_SELFTEST`
- Sub-library specific unit-tests.

Note that, ctest has a default timeout of 1500s, which is applied to all tests except the KAT tests. To override the default timeout, run `ctest --timeout <seconds>`.

## Known Answer Tests (KAT)

KAT are available in folder `KAT`. They can be generated by running the apps built in the `apps` folder:

- `apps/PQCgenKAT_sign_<level>`

A successful execution will generate the `.req` and `.rsp` files.

A full KAT test is done as part of the test harness (see previous section).

## Benchmarks

A benchmarking suite is built and runs with the following command:

- `test/sqisign_bench_<param> <runs>`, where params specifies the SQIsign parameter set and runs the number of benchmark runs.

The benchmarks profile the `KeyGen`, `Sign` and `Verify` functions. The results are reported in CPU cycles if available on the host platform, and timing in nanoseconds otherwise.

## Examples

Example code that demonstrates how to use SQIsign are available in the `apps` folder:

- `apps/example_nistapi.c`: Example with the NIST API.

## Project Structure

The SQIsign library consists of a number of sub-libraries used to implement the final SQIsign library.

The dependencies are depicted below.
```
    ┌─┬──────┬─┐           ┌─┬────┬─┐            ┌─┬──────┬─┐
    │ ├──────┤ │           │ ├────┤ │            │ ├──────┤ │
    │ │Keygen│ │           │ │Sign│ │            │ │Verify│ │
    │ ├──────┤ │           │ ├────┤ │            │ ├──────┤ │
    └─┴───┬──┴─┘           └─┴─┬──┴─┘            └─┴───┬──┴─┘
          │                    │                       │
          │                    │                       │
          ├────────────────────┼─────────────────┐     │
          │                    │                 │     │
          │                    │                 │     │
      ┌───▼──┐          ┌──────▼────────┐   ┌────▼─────▼───────────┐
      │ PRNG ◄────┬─────┤ Iso <-> Ideal ├───►   Elliptic Curves,   │
      └───▲──┘    │     └──────┬────────┘   │ Pairings & Isogenies │
          │       │            │            └───▲──────┬───────────┘
          │       │            │                │      │
      ┌───┴──┐    │            │                │      │
      │ KLPT ◄────┘            │     ┌──────────┘      │
      └───┬──┘                 │     │                 │
          │                    │     │                 │
┌─────────▼─────────┐          │     │                 │
│ Quaternion orders │          │     │            ┌────▼───┐
│     and ideals    │          │     │            │ GF(p²) │
└─────────┬─────────┘          │     │            └────┬───┘
          │           ┌─┬──────▼─────┴──┬─┐            │
    ┌─────▼─────┐     │ ├───────────────┤ │      ┌─────▼─────┐
    │ MP BigInt │     │ │Precomputations│ │      │ FP BigInt │
    └───────────┘     │ ├───────────────┤ │      └───────────┘
                      └─┴───────────────┴─┘                       
```

There are the following sub-libraries:

- `common`: common code for AES, SHAKE, (P)RNG, memory handling
- `ec`: elliptic curves, isogenies and pairings
- `gf`: GF(p^2) and GF(p) arithmetic, including FP BigInt
- `id2iso`: code for Iso <-> Ideal
- `klpt`: implementation of KLPT
- `quaternion`: quaternion orders and ideals
- `intbig`: multi-precision big integers
- `precomp`: precomputed constants
- `protocols`: protocol implementation


### Folder structure

Folder levels after `src`:
```
SQIsign
└── src
    ├── lib_1
    │   ├── broadwell
    │   │   ├── generic
    │   │   └── lvl1
    │   ├── opt
    │   │   ├── generic
    │   │   └── lvl1
    │   └── ref
    │       └── generic
    ├── lib_2
    │   ├── broadwell
    │   │   └── generic
    │   ├── opt
    │   │   └── generic
    │   └── ref
    │       └── generic
    └── lib_n
        ├── broadwell
        │   └── generic
        ├── opt
        │   └── generic
        └── ref
            └── generic
```

Level 1: Library (e.g. quaternion). A `CMakeLists.txt` file with entry `include(${SELECT_IMPL_TYPE})` takes care of including the implementation Level 2.

Level 2: Implementation type: reference C (ref), optimized C (opt), ASM-optimized (e.g. broadwell, neon, m4). A `CMakeLists.txt` file entry with `include(${SELECT_SQISIGN_VARIANT})` takes care of including the SQIsign variant.

Level 3: SQIsign variant -> generic code or code for a specific parameter set (e.g. lvl1). 


Other folders:
- `apps`: Applications: KAT generation application, examples
- `include`: SQIsign public header files
- `KAT`: Known Answer Test files
- `test`: SQIsign test code

### Sub-library headers

Sub-libraries can define their own headers, which may be different between the implementation types. These header files are used sub-library-internally and by other dependent sub-libraries. The convention is to put the headers in an `include` folder of the sub-library src directory. For example, `src/intbig/ref/generic/include/intbig.h`.

### Sub-library unit tests

Sub-libraries can implement their own, self-contained unit tests. The convention is to put the unit tests in a `test` folder of the sub-library `src` directory. For example, `src/intbig/ref/generic/test/test_intbig.c`.

### Shared implementation types

It is possible to share implementations between implementation types. For example, the broadwell optimized implementation might use the same code as the reference implementation except in the GF module.

## License

SQIsign is licensed under Apache-2.0. See [LICENSE](LICENSE) and [NOTICE](NOTICE).

Third party code is used in some test and common code files:

- `src/common/aes_c.c`; MIT: "Copyright (c) 2016 Thomas Pornin <pornin@bolet.org>"
- `src/common/fips202.c`: Public Domain
- `src/common/randombytes_system.c`: MIT: Copyright (c) 2017 Daan Sprenkels <hello@dsprenkels.com>
- `apps/PQCgenKAT_sign.c`, `common/randombytes_ctrdrbg.c`, `test/test_kat.c`: by NIST (Public Domain)
