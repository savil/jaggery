const std = @import("std");

const FieldElement = struct {

};

/// xgcd is the extended euclidean algorithm for computing multiplicative inverses:
/// for input a and b, returns their greatest common divisor g, 
/// along with Bezout coefficients a,b such that ax + by = g.
/// https://aszepieniec.github.io/stark-anatomy/basic-tools#finite-field-arithmetic
///
/// The multiplicative inverse of A is A^-1 such that:
/// - A * A^-1 = 1 mod C or equivalently (A * A^-1) mod C = 1
/// https://www.khanacademy.org/computing/computer-science/cryptography/modarithmetic/a/modular-inverses
///
const xgcdResult = struct {a: i64, b: i64, g: i64};
fn xgcd(x_in: i64, y_in: i64) !xgcdResult {
    var x = x_in;
    var y = y_in;
    var a: i64 = 1;
    var new_a: i64 = 0;
    var b: i64 = 0;
    var new_b: i64 = 1;

    while (y != 0) {
        var quotient = try std.math.divFloor(@TypeOf(x), x, y);

        var old_x = x;
        x = y;
        y = old_x - quotient * y; // this is the remainder

        var old_a = a;
        a = new_a;
        new_a = old_a - quotient * new_a;

        var old_b = b;
        b = new_b;
        new_b = old_b - quotient * new_b;
    }

    return xgcdResult{.a = a, .b = b, .g = x};
}

test "xgcd" {
    var result = try xgcd(1, 1);
    try std.testing.expectEqual(@as(i64, 1), result.g);

    result = try xgcd(1, 2);
    try std.testing.expectEqual(@as(i64, 1), result.g);

    result = try xgcd(10, 4);
    try std.testing.expectEqual(@as(i64, 2), result.g);

    // example from:
    // https://www.khanacademy.org/computing/computer-science/cryptography/modarithmetic/a/the-euclidean-algorithm
    result = try xgcd(270, 192);
    try std.testing.expectEqual(@as(i64, 6), result.g);
}