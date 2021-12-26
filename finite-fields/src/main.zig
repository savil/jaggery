const std = @import("std");
const testing = std.testing;

export fn add(a: i32, b: i32) i32 {
    return a + b;
}

test "basic add functionality" {
    try testing.expect(add(3, 7) == 10);
}

/// this field must have a sufficiently large
/// subgroup of power-of-two order
fn field() Field {
    p = 1 + 407 * ( 1 << 119 ); // 1 + 11 * 37 * 2^119;
    return Field.init(p);
}