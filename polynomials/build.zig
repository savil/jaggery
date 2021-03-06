const std = @import("std");

pub fn build(b: *std.build.Builder) void {
    // Standard release options allow the person running `zig build` to select
    // between Debug, ReleaseSafe, ReleaseFast, and ReleaseSmall.
    const mode = b.standardReleaseOptions();

    const lib = b.addStaticLibrary("polynomials", "src/main.zig");
    lib.addPackagePath("finite-fields", "../finite-fields/src/main.zig");
    lib.setBuildMode(mode);
    lib.install();

    const main_tests = b.addTest("src/main.zig");
    main_tests.setBuildMode(mode);

    const polynomial_tests = b.addTest("src/polynomial.zig");
    polynomial_tests.addPackagePath("finite-fields", "../finite-fields/src/main.zig");
    polynomial_tests.setBuildMode(mode);
    
    const test_step = b.step("test", "Run library tests");
    test_step.dependOn(&main_tests.step);
    test_step.dependOn(&polynomial_tests.step);
}
